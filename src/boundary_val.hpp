#ifndef __RANDWERTE_HPP__
#define __RANDWERTE_HPP__

#include "cstring"
#include "helper.hpp"
#include "datastructures.hpp"
#include "grid.hpp"
#include <Eigen/Dense>

using namespace Eigen;

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(Grid& grid, Config& config);

void boundaryvalues_t(Grid &grid, Config& config, MatrixXd& T);
#endif
