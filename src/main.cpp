#include <stdio.h>

#include <iostream>
#include <ctime>

#include <mpi.h>

#include "helper.hpp"
#include "simulation.hpp"

#ifdef gpp9
// gcc Version >= 9
#include "filesystem"
namespace filesystem = std::filesystem;
#else
// gcc Version < 9
#include "experimental/filesystem"
namespace filesystem = std::experimental::filesystem;
#endif // DEBUG

inline void print_scenarios(filesystem::path folder) {
  std::string file;
  for (auto entry : filesystem::directory_iterator(folder)) {
    file = std::string(entry.path().filename());
    std::size_t dat = file.find(std::string(".dat"));
    if (dat != std::string::npos) {
      printf("- %s \n", file.substr(0, dat).c_str());
    }
  }
}

int main(int argn, char** args) {

  int rank, input_size, output_size;

  MPI_Init(&argn, &args);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  filesystem::path input, output;

  if (rank == 0) {
    std::string input_path = get_cmd_argument(argn, args, "-in");
    std::string output_path = get_cmd_argument(argn, args, "-out");
    
    if(input_path.empty()) {
      input_path = "../scenarios";
    }

    if(output_path.empty()) {
      output_path = "../vtk-files";
    }

    input = filesystem::path(input_path);
    
    if(!filesystem::exists(input)) {
      printf("\u001B[31m Input path %s not found!\u001B[0m \n", input.c_str());
			MPI_Abort(MPI_COMM_WORLD, -1);
    }

    if (argn == 1) {
      // no problem name given, exit program with error message
      printf("Please give the scenario name that should be simulated!\n");
      print_scenarios(input);
      printf("\n example: ./sim Cavity100 \n");
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // take second argument since first is the program name
    std::string scenario = std::string(args[1]);

    input /= (scenario + ".dat");

    // check if .dat file for scenario exists
    if(!filesystem::exists(input)) {
      // no .dat file found
      printf("\u001B[31m No .dat file found for scenario '%s' in input folder '%s'! \u001B[0m \n", scenario.c_str(), input.c_str());
      MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // create output path
    output = filesystem::path(output_path);
    output /= (scenario + "-" + std::to_string(time(0)));
    filesystem::create_directories(output);

    // add scenario name so vtk files have the correct name
    output /= scenario; // TODO: maybe this should be moved inside run_simulation and use the 'problem' parameter

    input_size = input.string().size();
    output_size = output.string().size();
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(&input_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&output_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  std::string in = std::string(input_size, ' ');
  std::string out = std::string(output_size, ' ');
  
  if(rank == 0) {
    in = input.string();
    out = output.string();
  }

  MPI_Bcast(const_cast<char*>(in.data()), input_size, MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(const_cast<char*>(out.data()), output_size, MPI_CHAR, 0, MPI_COMM_WORLD);

  // append rank to output path
  sprintf(out.data(), "%s_p%d", out.c_str(), rank);

  input = filesystem::path(in.c_str());
  output = filesystem::path(out.c_str());

  MPI_Barrier(MPI_COMM_WORLD);
  
  run_simulation(input, output);

  MPI_Finalize();

  return 0;
}
