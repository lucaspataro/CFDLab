#ifndef __INIT_H_
#define __INIT_H_

#include <string>
#include "datastructures.hpp"

#ifdef gpp9
// gcc Version >= 9
#include "filesystem"
namespace filesystem = std::filesystem;
#else
// gcc Version < 9
#include "experimental/filesystem"
namespace filesystem = std::experimental::filesystem;
#endif // DEBUG

/**
 * Maximum length of input lines
 */
#define MAX_LINE_LENGTH 1024

/**
 * This operation initializes all the local variables reading a configuration
 * file. For every variable a macro like READ_INT() is called passing it the
 * input filename and the variable to be written to. This macro calls
 * an operation read_int() augmenting the parameter set with the name of the
 * variable to be read. The read_int() operation parses the input file, extracts
 * the value of the variable, sets the variable and finally prints some debug
 * information. This is possible as the macro adds the name of the variable to
 * be set. All the helper operations can be found within helper.hpp and helper.c.
 *
 * @param szFileName char pointer to the filename
 * @param config Config struct containing all configuration options
 */
int read_parameters( 
  filesystem::path szFileName,
  Config& conf
);


#endif

