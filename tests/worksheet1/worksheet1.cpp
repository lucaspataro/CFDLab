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

void cavity100(std::string input, std::string output){
    int argn = 0; 
    char **args = new char*[0];

    MPI_Init(&argn, &args);
    run_simulation(filesystem::path(input), filesystem::path(output));
    MPI_Finalize();
    delete [] args;
}


TEST_CASE( "Check output for worksheet 1", "[ws1 ws1_test]" ) {
    // create output folder
    filesystem::create_directories(filesystem::path("worksheet1/output/"));

    // run simulation
    cavity100("worksheet1/Cavity100.dat","worksheet1/output/Cavity100");

    // compare results
    REQUIRE( compareFiles("worksheet1/compare/Cavity100.0.vtk","worksheet1/output/Cavity100.0.vtk"));
    REQUIRE( compareFiles("worksheet1/compare/Cavity100.1.vtk","worksheet1/output/Cavity100.1.vtk"));
    REQUIRE( compareFiles("worksheet1/compare/Cavity100.2.vtk","worksheet1/output/Cavity100.2.vtk"));
    REQUIRE( compareFiles("worksheet1/compare/Cavity100.3.vtk","worksheet1/output/Cavity100.3.vtk"));
    REQUIRE( compareFiles("worksheet1/compare/Cavity100.4.vtk","worksheet1/output/Cavity100.4.vtk"));
    REQUIRE( compareFiles("worksheet1/compare/Cavity100.5.vtk","worksheet1/output/Cavity100.5.vtk"));
}
