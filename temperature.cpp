//
// Created by felix on 23.05.20.
//
#include "temperature.hpp"
#include "grid.hpp"
#include <iostream>
using namespace std;

void calculate_t(Grid& grid, Config& config, MatrixXd& T, MatrixXd& U, MatrixXd& V, matrix<cell_type>& types) {
    double dTdt, duTdx, dvTdy, d2Tdx2, d2Tdy2;
    double T_ij, T_ip1j, T_im1j, T_ijm1, T_ijp1;
    double u_ij, u_im1j;
    double v_ij, v_ijm1;
    
    MatrixXd T_hat = MatrixXd::Zero(T.rows(), T.cols());

    for (auto i = 1; i < grid.imaxb()-1; i++) {
        for (auto j = 1; j < grid.jmaxb()-1; j++) {
            if(types[i][j] == cell_type::FLUID) {
                T_ij = T.coeffRef(i,j);
                T_ip1j= T.coeffRef(i + 1,j);
                T_im1j = T.coeffRef(i - 1,j);
                T_ijm1 = T.coeffRef(i,j-1);
                T_ijp1 = T.coeffRef(i,j+1);
                u_ij = U.coeffRef(i,j);
                u_im1j = U.coeffRef(i-1,j);
                v_ij = V.coeffRef(i,j);
                v_ijm1 = V.coeffRef(i,j-1);

                duTdx = (u_ij * (T_ij + T_ip1j) - u_im1j * (T_im1j + T_ij) + config.alpha * (abs(u_ij) * (T_ij - T_ip1j) - abs(u_im1j) * (T_im1j - T_ij))) * 0.5 / config.dx;

                dvTdy = (v_ij * (T_ij + T_ijp1) - v_ijm1 * (T_ijm1 + T_ij) + config.alpha * (abs(v_ij) * (T_ij - T_ijp1) - abs(v_ijm1) * (T_ijm1 - T_ij))) * 0.5 / config.dy;

                d2Tdx2 = (T_ip1j - 2 * T_ij + T_im1j) / (config.dx * config.dx);

                d2Tdy2 = (T_ijp1 - 2 * T_ij + T_ijm1) / (config.dy * config.dy);

                dTdt = ((d2Tdx2 + d2Tdy2)/(config.PR * config.Re)) - duTdx - dvTdy;

                T_hat(i,j) = T.coeffRef(i,j) + config.dt * dTdt;
            }
        }
    }

    grid.set_temperature(T_hat);
}
