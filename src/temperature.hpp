#ifndef __TEMPERATUR_H_
#define __TEMPERATUR_H_
#include "grid.hpp"

void calculate_t(Grid& grid, Config& config, MatrixXd& T, MatrixXd& U, MatrixXd& V, matrix<cell_type>& type);

#endif
