#include <Eigen/Dense>
#include <iostream>

#include "solver.hpp"
#include "multigrid.hpp"
#include "sor.hpp"
#include "helper.hpp"
#include "CG.hpp"
#include "spLU.hpp"

using namespace Eigen;

Solver* Solver::create(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types) {
    switch(config.solver) {
        case solver_type::MULTIGRID:
            return new Multigrid(config, p, rhs, types);
            break;

        case solver_type::SOR:
            return new SOR(config, p, rhs, types);
            break;
        case solver_type::CG:
            return new CG(config, p, rhs, types);
            break;
        case solver_type::spLU:
            return new spLU(config, p, rhs, types);
            break;
    }

    return nullptr;
}

void Solver::comp_residual(
        MatrixXd &x,
        MatrixXd &b,
        MatrixXd &res,
        double dx,
        double dy
) {
    //std::cout << "-----------------------------------" << std::endl;
    // print_matrix(res, res.rows() - 2, res.cols() - 2);
    /* compute the residual */
    double dxdx_inv = 1.0 / (dx * dx);
    double dydy_inv = 1.0 / (dy * dy);
    double xip1j, xij, xi_1j, xijp1, xij_1;
    double rs = 0;
    for (auto i = 1; i < x.rows()-1; i++) {
        for (auto j = 1; j < x.cols()-1; j++) {
            xip1j = x.coeffRef(i+1,j);
            xij = x.coeffRef(i,j);
            xi_1j = x.coeffRef(i-1,j);
            xijp1 = x.coeffRef(i,j+1);
            xij_1 = x.coeffRef(i,j-1);
            rs = b.coeffRef(i,j);
            res.coeffRef(i,j) = rs - (xip1j-2.0*xij+xi_1j)*dxdx_inv - ( xijp1-2.0*xij+xij_1)*dydy_inv;
        }
    }
}

void Solver::compute_l2Norm(double* residual, MatrixXd& residual_) {
    double res = 0.0;
    for (int j = 0; j < residual_.cols(); j++) {
        for (int i = 0; i < residual_.rows(); i++) {
            res += residual_(i, j) * residual_(i, j);
        }
    }

    res /= (residual_.rows()-2)*(residual_.cols()-2);
    res = sqrt(res);

    *residual = res;
}

void Solver::compute_average(double* residual, MatrixXd& residual_) {
    double res = 0.0;
    for (int j = 0; j < residual_.cols(); j++) {
        for (int i = 0; i < residual_.rows(); i++) {
            res += abs(residual_(i, j));
        }
    }

    res /= (residual_.rows()-2)*(residual_.cols()-2);

    *residual = res;
}

solver_type Solver::string2type(std::string& t) {
    if (t == "multigrid") {
        return solver_type::MULTIGRID;
    } else if (t == "cg") {
        return solver_type::CG;
    } else if ( t == "spLU") {
        return solver_type::spLU;
    } else {
        return solver_type::SOR;
    }
}


Solver::Solver(Config &config, matrix<cell_type>& types): _config(config), _types(types) {}

void Solver::boundary_p(
        MatrixXd &P,
        matrix<cell_type>& type
) {
    bool isBoundary = false;
    int imax = P.rows()-2;
    int jmax = P.cols()-2;

    /* set boundary values */
    /*
     * this sets the pressure for the ghost cells of each process if they are a boundary
     * this is meaningless for edges between prcessors in the global domain as it they will be overwritten
     * from the processor next to you. However, you need it for the boundaries of the global domain.
    */
    for(auto i = 1; i <= imax; i++) {
        isBoundary = (
                type[i][0] != cell_type::FLUID &&
                type[i][0] != cell_type::OUTLET
        );
        if (isBoundary && type[i][1] == cell_type::FLUID) {
            P.coeffRef(i,0) = P.coeffRef(i,1);
        }
        isBoundary = (
                type[i][jmax+1] != cell_type::FLUID &&
                type[i][jmax+1] != cell_type::OUTLET
        );
        if (isBoundary && type[i][jmax] == cell_type::FLUID) {
            P.coeffRef(i,jmax+1) = P.coeffRef(i,jmax);
        }
    }

    for(auto j = 1; j <= jmax; j++) {
        isBoundary = (
                type[0][j] != cell_type::FLUID &&
                type[0][j] != cell_type::OUTLET
        );
        if (isBoundary && type[1][j] == cell_type::FLUID) {
            P.coeffRef(0,j) = P.coeffRef(1,j);
        }
        isBoundary = (
                type[imax+1][j] != cell_type::FLUID &&
                type[imax+1][j] != cell_type::OUTLET
        );
        if (isBoundary && type[imax][j] == cell_type::FLUID) {
            P.coeffRef(imax+1,j) = P.coeffRef(imax,j);
        }
    }

    for(auto j = 1; j <= jmax; j++) {
        for (auto i = 1; i <= imax; i++) {
            if (type[i][j] == cell_type::NOSLIP) {
                P.coeffRef(i,j) = 0;
                int count = 0;
                if (type[i + 1][j] == cell_type::FLUID) {
                    P.coeffRef(i,j) += P.coeffRef(i + 1,j);
                    count++;
                }
                if (type[i-1][j] == cell_type::FLUID) {
                    P.coeffRef(i,j) += P.coeffRef(i - 1,j);
                    count++;
                }
                if (type[i][j+1] == cell_type::FLUID) {
                    P.coeffRef(i,j) += P.coeffRef(i,j+1);
                    count++;
                }
                if (type[i][j-1] == cell_type::FLUID) {
                    P.coeffRef(i,j) += P.coeffRef(i,j-1);
                    count++;
                }
                if (count != 0) {
                    P.coeffRef(i,j) /= count;
                }
            }
        }
    }
}

