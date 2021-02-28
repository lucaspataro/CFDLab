#include "sor.hpp"
#include <cmath>
#include "mpi.h"
#include "parallel.hpp"
#include <Eigen/Dense>
#include "helper.hpp"
using namespace Eigen;

SOR::SOR(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types): Solver(config, types), _x(p), _rhs(rhs) {
}

void SOR::solve(int& it, double& res) {
    it = 0;
    double r = 0.0, global_res;
    do {
        it++;
        // Perform a SOR iteration according to (18) using the provided function
        sor();

        // apply pressure boundary
        boundary_p(_x, _types);

        // retrieve the residual res
        calculate_res(r);

        // communicate pressure 
        boundary_comm(_config, _x, _config.imax + 2, _config.jmax + 2, 1);
       
        // communicate residual
        MPI_Allreduce(&r, &global_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        global_res /= _config.total_n_fluid;
        global_res = sqrt(global_res);

    } while (it < _config.itermax && global_res > _config.eps);
    res = global_res;
}

void SOR::sor() {
    MatrixXd& P = _x;
    MatrixXd& RS = _rhs;
    double dxdx_inv = 1.0/(_config.dx*_config.dx);
    double dydy_inv = 1.0/(_config.dy*_config.dy);
    double coeff = _config.omg/(2.0*(dxdx_inv + dydy_inv));
    /* SOR iteration */
    for(auto i = 1; i <= _config.imax; i++) {
        for(auto j = 1; j<= _config.jmax; j++) {
            if (_types[i][j] == cell_type::FLUID) {
                P.coeffRef(i,j) = (1.0-_config.omg)*P.coeffRef(i,j)
                                + coeff*(( P.coeffRef(i+1,j)+P.coeffRef(i-1,j))* dxdx_inv + ( P.coeffRef(i,j+1)+P.coeffRef(i,j-1))*dydy_inv - RS.coeffRef(i,j));
            }
        }
    }
}

void SOR::calculate_res(double& res) {
    MatrixXd& P = _x;
    MatrixXd& RS = _rhs;
    /* compute the residual */
    double dxdx_inv = 1.0/(_config.dx*_config.dx);
    double dydy_inv = 1.0/(_config.dy*_config.dy);
    double Pip1j,Pij,Pi_1j,Pijp1,Pij_1,rs;

    double rloc = 0;
    for(auto i = 1; i <= _config.imax; i++) {
        for(auto j = 1; j <= _config.jmax; j++) {
            if (_types[i][j] == cell_type::FLUID) {
                Pip1j = P.coeffRef(i+1,j);
                Pij = P.coeffRef(i,j);
                Pi_1j = P.coeffRef(i-1,j);
                Pijp1 = P.coeffRef(i,j+1);
                Pij_1 = P.coeffRef(i,j-1);
                rs = RS.coeffRef(i,j);
                rloc += ( (Pip1j-2.0*Pij+Pi_1j)*dxdx_inv + ( Pijp1-2.0*Pij+Pij_1)*dydy_inv - rs) *
                        ( (Pip1j-2.0*Pij+Pi_1j)*dxdx_inv + ( Pijp1-2.0*Pij+Pij_1)*dydy_inv - rs);
            }
        }
    }
    
    /* set residual */
    res = rloc;
}


