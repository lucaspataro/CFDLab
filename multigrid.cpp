#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <Eigen/Cholesky>

#include <iostream>
#include <vector>
#include <algorithm> 

#include "sor.hpp"
#include "enums.hpp"
#include "multigrid.hpp"
#include "helper.hpp"

#define SMOOTHING_ITERATIONS 3

using namespace Eigen;

template<typename T>
void slice(matrix<T>& org, matrix<T>& sliced, int nth, int imax, int jmax) {
    
    sliced.resize(imax, std::vector<T>(jmax, cell_type::FLUID));
    for(auto i=0, x=0; i < org.size(); i+=nth, x++) {
        for(auto j=0, y=0; j < org[0].size(); j+=nth, y++) {
            sliced[x][y] = org[i][j];
        }
    }
}

Level::Level(int level_, int imax_, int jmax_, double dx_, double dy_, matrix<cell_type>& types_):
    level(level_), imax(imax_), jmax(jmax_), 
    dx(dx_), dy(dy_) {
    
    _types = &types_;
    _x = new MatrixXd(imax_, jmax_);
    _b = new MatrixXd(imax_, jmax_);
    _e = new MatrixXd(imax_, jmax_);
    _res = new MatrixXd(imax_, jmax_);
}

Level::Level(int level_, int imax_, int jmax_, double dx_, double dy_, matrix<cell_type>& types_, MatrixXd& p, MatrixXd& rhs):
    level(level_), imax(imax_), jmax(jmax_), 
    dx(dx_), dy(dy_) {
    
    _types = &types_;
    _x = &p;
    _b = &rhs;

    _e = new MatrixXd(imax_, jmax_);
    _res = new MatrixXd(imax_, jmax_);
}

MatrixXd& Level::x() {
    return *_x;
}

MatrixXd& Level::b() {
    return *_b;
}

MatrixXd& Level::e() {
    return *_e;
}

MatrixXd& Level::res() {
    return *_res;
}

matrix<cell_type>& Level::types() {
    return *_types;
}

Multigrid::Multigrid(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types):
    Solver(config, types) {
    
    int imax = p.rows();
    int jmax = p.cols();
    double dx = config.dx;
    double dy = config.dy;

    _residual = MatrixXd::Zero(imax, jmax);

    for (auto i=0; i < config.levels; ++i) {
        
        if(i == 0) {
            _levels.push_back(new Level(i, imax, jmax, dx, dy, types, p, rhs));
        } else {
            matrix<cell_type>* S = new matrix<cell_type>(imax, std::vector<cell_type>(jmax, cell_type::FLUID));
            
            slice(types, *S, pow(2, i), imax, jmax);
            _levels.push_back(new Level(i, imax, jmax, dx, dy, *S));
        }

        if(i < config.levels - 1) {
            imax = (imax - 1) / 2 + 1;
            jmax = (jmax - 1) / 2 + 1;

            dx = 1.0 / (imax - 1.0);
            dy = 1.0 / (jmax - 1.0);
        }
    }

    // build system matrix for coarsest grid
    imax = imax - 2;
    jmax = jmax - 2;
    
    double h_xx_inv = 1 / (dx*dx);
    double h_yy_inv = 1 / (dy*dy);
    int nx = 2;
    int ny = 2;
    int node = 0;

    _A = MatrixXd::Zero(imax*jmax, imax*jmax);

    for(auto i = 0; i < imax; ++i) {
        for (auto j = 0; j < jmax; ++j) {
            nx = 2;
            ny = 2;
            node = j*imax+i;
            if(j != 0) {
                // lower node
                _A(node, node - imax)= -h_yy_inv;
            } else {
                if(_levels[config.levels-1]->types()[i+1][j] != cell_type::OUTLET) {
                    ny = ny-1;
                }
            }

            if(j != jmax - 1) {
                // upper node
                _A(node, node + imax)= -h_yy_inv;
            } else {
                if(_levels[config.levels-1]->types()[i+1][j+2] != cell_type::OUTLET) {
                    ny = ny-1;
                }
            }

            if(i != 0) {
                // left node
                _A(node, node - 1)= -h_xx_inv;
            } else {
                if(_levels[config.levels-1]->types()[i][j+1] != cell_type::OUTLET) {
                    nx = nx-1;
                }
            }

            if(i != imax - 1) {
                // right node
                _A(node, node + 1)= -h_xx_inv;
            } else {
                if(_levels[config.levels-1]->types()[i+2][j+1] != cell_type::OUTLET) {
                    nx = nx-1;
                }
            }

            _A(node, node)= nx*h_xx_inv + ny*h_yy_inv;
        }
    }
    //_A.makeCompressed();
    _solver = FullPivLU<MatrixXd>(_A);
}

Level& Multigrid::levels(int l) {
    return *_levels[l];
}

void Multigrid::solve(int& it, double& res) {
    res = 0.0;
    it = 0;
    Level fine = levels(0);
    do {
        v_cycle(0);

        // apply pressure boundary
        boundary_p(fine.x(), fine.types());

        // retrieve the residual res
        comp_residual(fine.x(), fine.b(), _residual, _config.dx, _config.dy);

        compute_l2Norm(&res, _residual);
        it += 1;
    } while (it < _config.itermax && res > _config.eps);
}

void Multigrid::v_cycle(int depth) {
    Level level = levels(depth);
    
    // print_matrix(level.types(), level.imax-2, level.jmax-2);

    if (depth < _config.levels - 1) {
        Level next = levels(depth+1);
        
        // presmoothing
        smoothing(level.x(), level.b(), level.types(), level.dx, level.dy, depth == 0);

        // compute residual
        level.res().setZero();
        comp_residual(level.x(), level.b(), level.res(), level.dx, level.dy);

        // restrict residual
        next.b().setZero();
        restriction_fullweight(level.res(), next.b());


        // go one level deeper
        next.x().setZero();
        v_cycle(depth + 1);

        // prolongate residual
        level.e().setZero();
        prolongate(level.e(), next.x());

        // add 
        level.x() = level.x() + level.e();
        // Post-Smoothing
        smoothing(level.x(), level.b(), level.types(), level.dx, level.dy, depth == 0);

    } else {
        // use gaussSeidel until it works afterwards direct solve
        /**
        double res = 0.0;
        
        for (auto i = 0; i < _config.itermax; ++i) {
            gaussSeidel(level.x(), level.b(), level.types(), level.dx, level.dy);
            boundary_p(level.x(), _types);

            comp_residual(level.x(), level.b(), level.res(), level.dx, level.dy);

            compute_l2Norm(&res, level.res());

            if(res <= _config.eps) {
                break;
            }
        }**/

        // direkt solve
        // level.x().setZero();

        int rows = level.imax - 2;
        int cols = level.jmax - 2;

        MatrixXd x_block = level.x().block(1, 1, rows, cols);
        MatrixXd b_block = level.b().block(1, 1, rows, cols);

        x_block.resize(rows * cols, 1);
        b_block.resize(rows * cols, 1);

        x_block = _solver.solve(b_block);

        x_block.resize(rows, cols);

        level.x().block(1, 1, rows, cols) = -x_block;
        boundary_p(level.x(), level.types());
    }
}

void Multigrid::smoothing(
        MatrixXd &P,
        MatrixXd &RS,
        matrix<cell_type>& types,
        double dx,
        double dy,
        bool boundary
) {
    // do 1 - 3 steps of gauss-seidel
    for (auto i = 0; i < SMOOTHING_ITERATIONS; ++i) {
        gaussSeidel(P, RS, types, dx, dy);
        boundary_p(P, types);
    }
}

/**
void restriction_semi(MatrixXd &m, MatrixXd &b) {
    checkCoarseningAllowed(m.cols(), m.rows());
    auto imax = (m.cols()-1) / 2;
    auto jmax = (m.rows()-1) / 2;
    assert(imax  == b.cols());

    auto x = 0;
    auto y = 0;
// Question: doesnt write in first element...
    for (auto j = 0; j < jmax; j++) {
        for (auto i = 0; i < imax; i++) {
            x = 2 * i;
            y = 2 * j;
            b(i, j) = 0.25 * (m.coeffRef(x+1,y)+m.coeffRef(x,y+1)+m.coeffRef(x+2,y+1)+m.coeffRef(x+1,y+2)) +
                      0.5 * m.coeffRef(x+1,y+1);
        }
    }
}**/


void Multigrid::restriction_fullweight(MatrixXd &m, MatrixXd &b) {
    auto imax = (m.cols()-1) / 2;
    auto jmax = (m.rows()-1) / 2;

    auto x = 0;
    auto y = 0;

    for (auto j = 1; j < jmax; j++) {
        for (auto i = 1; i < imax; i++) {
            x = 2 * i;
            y = 2 * j;
            b.coeffRef(i, j) = 1/16 * (m.coeffRef(x-1, y-1) + m.coeffRef(x-1, y+1) + m.coeffRef(x+1, y-1) + m.coeffRef(x + 1, y + 1)) +
                    1/8 * (m.coeffRef(x+1,y)+m.coeffRef(x,y+1)+m.coeffRef(x-1,y)+m.coeffRef(x,y-1)) +
                    0.25 * m.coeffRef(x,y);
        }
    }

    // left and right boundary
    for (auto j = 0; j< jmax+1; j++) {
        // left
        b.coeffRef(0,j) = m.coeffRef(0,2*j);
        //right
        b.coeffRef(imax,j) = m.coeffRef(m.rows()-1,2*j);
    }

    // bottom and top boundary
    for (auto i = 0; i<imax+1; i++) {
        // bottom
        b.coeffRef(i,0) = m.coeffRef(2*i,0);
        // top
        b.coeffRef(i,jmax) = m.coeffRef(2*i,m.rows()-1);
    }
}

void Multigrid::prolongate(MatrixXd &fine, MatrixXd &coarse) {
    auto imax_coarse = coarse.rows()-1;
    auto jmax_coarse = coarse.cols()-1;
    auto imax_fine = fine.rows()-1;
    auto jmax_fine = fine.cols()-1;

    auto x = 0;
    auto y = 0;
    for (auto j = 0; j < jmax_coarse ; j++) {
        y = 2 * j;
        for (auto i = 0; i < imax_coarse; i++) {
            x = 2 * i;
            fine.coeffRef(x, y) = coarse.coeffRef(i, j);
            fine.coeffRef(x + 1, y) = 0.5 * (coarse.coeffRef(i, j) + coarse.coeffRef(i+1,j));
            fine.coeffRef(x, y+1) = 0.5 * (coarse.coeffRef(i, j) + coarse.coeffRef(i, j+1));
            fine.coeffRef(x+1, y + 1) = 0.25 * (coarse.coeffRef(i, j) + coarse.coeffRef(i+1,j) + coarse.coeffRef(i, j+1) + coarse.coeffRef(i+1, j+1));
        }
    }

    for (auto j = 0; j< jmax_coarse; j++) {
        // right boundary
        fine.coeffRef(imax_fine,2*j) = coarse.coeffRef(imax_coarse,j); //center
        fine.coeffRef(imax_fine,2*j+1) = 0.5 * (coarse.coeffRef(imax_coarse,j) + coarse.coeffRef(imax_coarse,j+1)); //top
    }

    for (auto i = 0; i< imax_coarse; i++) {
        // top
        fine.coeffRef(2*i,jmax_fine) = coarse.coeffRef(i,jmax_coarse); //center
        fine.coeffRef(2*i+1,jmax_fine) = 0.5 * (coarse.coeffRef(i,jmax_coarse) + coarse.coeffRef(i + 1,jmax_coarse)); //right
    }

    // upper right corner
    fine.coeffRef(imax_fine,jmax_fine) = coarse.coeffRef(imax_coarse,jmax_coarse); //center
}

void Multigrid::gaussSeidel(MatrixXd &P, MatrixXd &RS, matrix<cell_type>& _types, double dx, double dy) {
    double rloc;
    double dxdx_inv = 1.0/(dx*dx);
    double dydy_inv = 1.0/(dy*dy);
    double coeff = 1/(2.0*(dxdx_inv + dydy_inv));
    int imax = P.rows()-2;
    int jmax = P.cols()-2;
    /* SOR iteration */
    for(auto i = 1; i <= imax; i++) {
        for(auto j = 1; j<= jmax; j++) {
            if (_types[i][j] == cell_type::FLUID) {
                P(i,j) = coeff*(( P(i+1,j)+P(i-1,j))* dxdx_inv + ( P(i,j+1)+P(i,j-1))*dydy_inv - RS(i,j));
            }
        }
    }
}