//
// Created by felix on 02.07.20.
//

#include "spLU.hpp"
#include <Eigen/Sparse>
#include <iostream>
#include "helper.hpp"

using namespace Eigen;


spLU::spLU(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types):
        Solver(config, types) {
    int imax = p.rows();
    int jmax = p.cols();
    double dx = config.dx;
    double dy = config.dy;

    _x = &p;
    _rhs = &rhs;


    // build system matrix for coarsest grid
    imax = imax - 2;
    jmax = jmax - 2;

    double h_xx_inv = 1 / (dx*dx);
    double h_yy_inv = 1 / (dy*dy);
    int nx = 2;
    int ny = 2;
    int node = 0;

    _A = SparseMatrix<double>(imax*jmax, imax*jmax);

    for(auto i = 0; i < imax; ++i) {
        for (auto j = 0; j < jmax; ++j) {
            nx = 2;
            ny = 2;
            node = j*imax+i;
            if(j != 0) {
                // lower node
                _A.insert(node, node - imax)= -h_yy_inv;
            } else {
                if(types[i+1][j] != cell_type::OUTLET) {
                    ny = ny-1;
                }
            }

            if(j != jmax - 1) {
                // upper node
                _A.insert(node, node + imax)= -h_yy_inv;
            } else {
                if(types[i+1][j+2] != cell_type::OUTLET) {
                    ny = ny-1;
                }
            }

            if(i != 0) {
                // left node
                _A.insert(node, node - 1)= -h_xx_inv;
            } else {
                if(types[i][j+1] != cell_type::OUTLET) {
                    nx = nx-1;
                }
            }

            if(i != imax - 1) {
                // right node
                _A.insert(node, node + 1)= -h_xx_inv;
            } else {
                if(types[i+2][j+1] != cell_type::OUTLET) {
                    nx = nx-1;
                }
            }

            _A.insert(node, node)= nx*h_xx_inv + ny*h_yy_inv;
        }
    }
    _solver.compute(_A);
}

void spLU::solve(int &it, double &res) {
    int rows = (*_x).rows()-2;
    int cols = (*_x).cols()-2;
    MatrixXd x_block = (*_x).block(1, 1, rows, cols);
    MatrixXd b_block = (*_rhs).block(1, 1, rows, cols);

    x_block.resize(rows * cols, 1);
    b_block.resize(rows * cols, 1);

    x_block = _solver.solve(b_block);
    x_block.resize(rows, cols);

    (*_x).block(1, 1, rows, cols) = -x_block;
    boundary_p(*_x, _types);
}