# Definition of the minimum required cmake Version
cmake_minimum_required(VERSION 3.23)
# Definition of the Project
# Later you can access the project variable like ${CFDLAB_SOURCE_DIR}
project(CFDLAB VERSION 1.0)

# Definition of the C++ Standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(MPI REQUIRED)
if (NOT MPI_FOUND)
  message("MPI not found")
  return ()
endif()
include_directories(SYSTEM ${MPI_INCLUDE_PATH})


find_package (Eigen3 3.3 REQUIRED)


# VTK Library
find_package(VTK REQUIRED)
if (NOT VTK_FOUND)
  message("VTK not found")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")


# Filesystem
# set(gpp9 False)
set(gpp9 True)
if(gpp9)
    message("g++ Version is 9 or higher")
    add_compile_definitions(gpp9)
else()
    message("g++ Version is lower than Version 9")
    # Link Filesystem Libraries
    link_libraries(stdc++fs)
endif()

# create a library that can be used form catch
#add_library(
#        CFD
#        main.cpp sor.cpp sor.hpp helper.cpp helper.hpp init.cpp init.hpp visual.cpp visual.hpp
#        boundary_val.cpp boundary_val.hpp uvp.cpp uvp.hpp grid.cpp grid.hpp cell.cpp cell.hpp
#        simulation.cpp simulation.hpp temperature.cpp temperature.hpp parallel.cpp parallel.hpp)
#
# make it public
# target_include_directories(CFD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} "/usr/include/" )

# if you use external libraries you have to link them like 
# target_link_libraries(CFD 
#  PUBLIC
#  ${VTK_LIBRARIES}
#  MPI::MPI_CXX
#  Eigen3::Eigen)



# Creating the executable of our project and the required dependencies
# the executable is called sim

aux_source_directory(src SOURCES)


add_executable(sim ${SOURCES})
  
target_include_directories(sim PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# if you use external libraries you have to link them like 
target_link_libraries(sim 
  PUBLIC
  ${VTK_LIBRARIES}
  MPI::MPI_CXX
  Eigen3::Eigen)

# If you write tests, you can include your subdirectory (in this case tests) as done here
# Testing
# enable_testing ()
# add_subdirectory(tests)
