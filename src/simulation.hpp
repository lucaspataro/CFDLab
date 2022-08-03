#ifndef __SIMULATION_H_
#define __SIMULATION_H_

#ifdef gpp9
// gcc Version >= 9
#include "filesystem"
namespace filesystem = std::filesystem;
#else
// gcc Version < 9
#include "experimental/filesystem"
namespace filesystem = std::experimental::filesystem;
#endif // DEBUG

void run_simulation(filesystem::path input, filesystem::path output);

#endif