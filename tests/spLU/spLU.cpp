// Definition of sample test
#include "../catch_config.hpp"
#include "../catch.hpp"
#include "helper.hpp"
#include "simulation.hpp"
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

void simulate(std::string input, std::string output){
    int argn = 0; 
    char **args = new char*[0];

    MPI_Init(&argn, &args);
    run_simulation(filesystem::path(input), filesystem::path(output));
    MPI_Finalize();
    delete [] args;
}

TEST_CASE( "Sparse Lu direct solver", "[spLU spLU_test]" ) {
    // create output folder
    filesystem::create_directories(filesystem::path("spLU/output/"));
    double epsilon = 0.00001;
    // run simulation
    simulate("spLU/NaturalConvection.dat","spLU/output/NaturalConvection");

    // compare results
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.0.vtk","spLU/output/NaturalConvection.0.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.1.vtk","spLU/output/NaturalConvection.1.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.2.vtk","spLU/output/NaturalConvection.2.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.3.vtk","spLU/output/NaturalConvection.3.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.4.vtk","spLU/output/NaturalConvection.4.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.5.vtk","spLU/output/NaturalConvection.5.vtk", epsilon));
}
