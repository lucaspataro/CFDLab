#include <stdio.h>

#include <chrono>
#include <iostream>
#include <ctime>
#include <Eigen/Dense>

#include "boundary_val.hpp"
#include "helper.hpp"
#include "init.hpp"
#include "uvp.hpp"
#include "sor.hpp"
#include "multigrid.hpp"
#include "visual.hpp"
#include "temperature.hpp"
#include "parallel.hpp"
#include "enums.hpp"
#include "solver.hpp"

#include "mpi.h"

#ifdef gpp9
// gcc Version >= 9
#include "filesystem"
namespace filesystem = std::filesystem;
#else
// gcc Version < 9
#include "experimental/filesystem"
namespace filesystem = std::experimental::filesystem;
#endif // DEBUG

using namespace Eigen;

/**
 * The main operation reads the configuration file, initializes the scenario and
 * contains the main loop. So here are the individual steps of the algorithm:
 *
 * - read the program configuration file using read_parameters()
 * - set up the matrices (arrays) needed. Use the predefined matrix<typename>
 * type and give initial values in the constructor.
 * - perform the main loop
 * - at the end: destroy any memory allocated and // print some useful statistics
 *
 * The layout of the grid is decribed by the first figure below, the enumeration
 * of the whole grid is given by the second figure. All the unknowns corresond
 * to a two-dimensional degree of freedom layout, so they are not stored inGrid::
 * arrays, but in a matrix.
 *
 * @image html grid.jpg
 *
 * @image html whole-grid.jpg
 *
 * Within the main loop, the following steps are required (for some of the
 * operations, a definition is defined already within uvp.h):
 *
 * - calculate_dt() Determine the maximal time step size.
 * - boundaryvalues() Set the boundary values for the next time step.
 * - calculate_fg() Determine the values of F and G (diffusion and confection).
 *   This is the right hand side of the pressure equation and used later on for
 *   the time step transition.
 * - calculate_rs()
 * - Iterate the pressure poisson equation until the residual becomes smaller
 *   than eps or the maximal number of iterations is performed. Within the
 *   iteration loop the operation sor() is used.
 * - calculate_uv() Calculate the velocity at the next time step.
 */
void run_simulation(filesystem::path input, filesystem::path output) {
  int rank, num_proc;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_proc);


  auto total = std::chrono::high_resolution_clock::now();
  auto start = std::chrono::high_resolution_clock::now();

  int hitIterMax = 0;
  int setup_time = 0;
  int calculation_time = 0;
  int output_time = 0;

  std::string problem;
  std::string geometry;

  Config config;

  config.num_proc = num_proc;
  config.rank = rank;
  config.boundary_size = 1;

  // Read the problem parameters
  read_parameters(input, config);

  if (config.num_proc > 1 && config.solver != solver_type::SOR) {
    printf("\u001B[31m Only solver SOR can be used for parallel calculation!");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  init_parallel(config);

  // abort if number of processes does not match the 
  // desired domain grid
  if (config.rank == 0 && (config.num_proc != config.jproc * config.iproc)) {
    printf("\u001B[31m Number of processes %d insufficient for domain decomposition iproc %d jproc %d. \n", config.num_proc, config.iproc, config.jproc);
    printf("Please set the number of processes to %d. \u001B[0m \n", config.iproc * config.jproc);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  printf("[%d/%d] omg_i %d, omg_j %d rankl %d, rankr %d, rankb %d, rankt %d \n", 
          config.rank, 
          config.num_proc, 
          config.omg_i, 
          config.omg_j, 
          config.rank_l, 
          config.rank_r, 
          config.rank_b, 
          config.rank_t);

  config.l_imax = config.ir - config.il + 1;
  config.l_jmax = config.jt - config.jb + 1;

  matrix<int> GeoArray = read_pgm(config);
  size_t cols = (config.l_imax + 2 * config.boundary_size);
  size_t rows = (config.l_jmax + 2 * config.boundary_size);

  config.imax = cols-2;
  config.jmax = rows-2;

  MatrixXd F(cols, rows);
  MatrixXd G(cols, rows);
  MatrixXd RS(cols, rows);
  MatrixXd U(cols, rows);
  MatrixXd V(cols, rows);
  MatrixXd P(cols, rows);
  MatrixXd T(cols, rows);

  matrix<cell_type> type;
  type.resize(cols, std::vector<cell_type>(rows, cell_type::FLUID));

  Grid grid(config, GeoArray, U, V, P, T, type);
  Solver* solver = Solver::create(config, P, RS, type);

  auto stop = std::chrono::high_resolution_clock::now();
  auto milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
          .count();

  setup_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
          .count();


  double t = 0.0, res = 1.0, global_res = 1.0, last = -config.dt_value, averageIterations = 0;
  int n = 0, it = 0, vtkFileCount = 0;

  // Initial output
  boundaryvalues(grid, config);

  if (config.calcTemp) {
    boundaryvalues_t(grid, config, T);

    boundary_comm(config, T, cols, rows, 2);
  }
  
  uv_comm(config, U, V, cols, rows);

  write_vtkFile(output.c_str(), vtkFileCount, config, U, V, P, T, GeoArray);
  vtkFileCount++;

  // sync before simulation start
  MPI_Barrier(MPI_COMM_WORLD);

  while (t <= config.t_end) {

    // // print progress update every 200 timesteps
    if ( n % 200 == 0 && rank == 0) {
      stop = std::chrono::high_resolution_clock::now();
      milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(stop - total)
          .count();
      printf("Current simulation time %f, elapsed real time %s, last iterations: %d, last dt: %f \n", t, format_duration(milliseconds).c_str(), it, config.dt);
    }

    start = std::chrono::high_resolution_clock::now();
    
    // Select dt according to (13)
    calculate_dt(grid, config);

    if (config.calcTemp != 0) {
        calculate_t(grid, config, T, U, V, type);
        
        boundaryvalues_t(grid, config, T);

        boundary_comm(config, T, cols, rows, 2);
    }

    // Compute Fn and Gn according to (9),(10),(17)
    calculate_fg(grid, config, F, G, U, V, T, type);

    FG_comm(config, F, G, cols, rows);
    
    // Compute the right-hand side rhs of the pressure equation (11)
    calculate_rs(config, F, G, RS);

    solver -> solve(it, res);

    averageIterations = averageIterations + it;
    if(it >= config.itermax && config.rank == 0){
        hitIterMax++;
        std::cout << "\u001B[31m";
        std::cout << "WARNING: Iteration ended due to reaching itermax at " << t;
        std::cout << " and a res of " << res;
        std::cout << "\u001B[0m" << std::endl;
    }

    // Compute un+1 and vn+1 according to (7),(8)
    calculate_uv(grid, config, F, G);
    
    uv_comm(config, U, V, cols, rows);
    
    // Set boundary values for u and v according to (14),(15)
    boundaryvalues(grid, config);

    stop = std::chrono::high_resolution_clock::now();
    calculation_time +=
        std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
            .count();

    // only write vtk file when `dt_value` simulation time has passed
    if ((t - last) >= config.dt_value ) {
      start = std::chrono::high_resolution_clock::now();

      // Output of u; v; p; t values for visualization, if necessary
      write_vtkFile(output.c_str(), vtkFileCount, config, U, V, P, T, GeoArray);

      stop = std::chrono::high_resolution_clock::now();
    
      output_time +=
        std::chrono::duration_cast<std::chrono::milliseconds>(stop - start)
            .count();

      milliseconds =
        std::chrono::duration_cast<std::chrono::milliseconds>(stop - total)
            .count();

      vtkFileCount++;

      last = t;
    }

    t += config.dt;
    n++;
  }

  // devide the sum of the iterations in the sor by the number of timesteps
  averageIterations /= n;


  stop = std::chrono::high_resolution_clock::now();
  milliseconds =
      std::chrono::duration_cast<std::chrono::milliseconds>(stop - total)
          .count();

  if (config.rank == 0) {
    printf("The average iterations in the SOR were: %f\n", averageIterations);
    printf("Final simulation time: %f \n", t);
    printf("Runtime: %s \n", format_duration(milliseconds).c_str());
    printf("Setup time: %s \n", format_duration(setup_time).c_str());
    printf("Calculation time: %s \n", format_duration(calculation_time).c_str());
    printf("Output time: %s \n", format_duration(output_time).c_str());

    if (hitIterMax > 0) {
        std::cout << "\u001B[31m";
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        std::cout<<"WARNING: The SOR exited " << hitIterMax<<" time(s) due to the maximum iteration condition!"<<std::endl;
        std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
        std::cout << "\u001B[0m" << std::endl;
    }
  }
}