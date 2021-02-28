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


TEST_CASE( "Check output for multi", "[multi multi]" ) {
    // create output folder
    filesystem::create_directories(filesystem::path("multi/output/"));

    // run simulation
    simulate("multi/NaturalConvection.dat","multi/output/NaturalConvection");

    double epsilon = 0.00001;

    // compare results
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.0.vtk","multi/output/NaturalConvection.0.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.1.vtk","multi/output/NaturalConvection.1.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.2.vtk","multi/output/NaturalConvection.2.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.3.vtk","multi/output/NaturalConvection.3.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.4.vtk","multi/output/NaturalConvection.4.vtk", epsilon));
    REQUIRE( compareFilesMinimum("solverCompareVTK/NaturalConvection.5.vtk","multi/output/NaturalConvection.5.vtk", epsilon));
}
