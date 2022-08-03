//
// Created by felix on 27.06.20.
//

#include <iostream>
#include "CG.hpp"
#include "sor.hpp"
#include "helper.hpp"

void CG::solve(int& it, double& residual) {
    it = 0;
    residual = 0;
    double alpha = 0;
    double beta = 0;
    double r_2 = 0;
    double r_2_old = 0;
    double dAd = 0;
    MatrixXd& d = *_d;
    MatrixXd& res = *_res;
    MatrixXd& Ad = *_Ad;
    // boundary has to be applied only once in the beginning as you work only on the residual
    boundary_p(*_x, _types);

    comp_residual(*_x, *_rhs, res, _config.dx, _config.dy);
    d = res;

    //double r_2_test = 0;


    do{
        it++;
        calc_Ad();
        r_2 = res.squaredNorm();

        dAd = 0;
        for (auto j = 1; j < d.cols()-1; j++) {
            for (auto i = 1; i < d.rows()-1; i++) {
                dAd = dAd + d.coeffRef(i,j) * Ad.coeffRef(i,j);
            }
        }


        alpha = r_2 / (dAd);
        // update x
        *_x = *_x + alpha * d;

        r_2_old = r_2;
        // update residual
        res = res - alpha * Ad;
        r_2 = res.squaredNorm();
        beta = r_2/(r_2_old);
        // update conjugate d
        d = res + beta * d;

        compute_l2Norm(&residual, res);
    } while (it < _config.itermax && residual > _config.eps);
}

CG::CG(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types): Solver(config, types) {
    int rows = p.rows();
    int cols = p.cols();
    _d = new MatrixXd(rows, cols);
    _Ad = new MatrixXd(rows, cols);
    _res = new MatrixXd(rows, cols);
    _x = &p;
    _rhs = &rhs;
    _dyy_inv = 1 / (_config.dy*_config.dy);
    _dxx_inv = 1 / (_config.dx*_config.dx);
}

void CG::calc_Ad() {
    MatrixXd& d = *_d;
    MatrixXd& Ad = *_Ad;
    for (int j = 1; j < d.cols()-1; j++) {
        for (int i = 1; i < d.rows()-1; i++) {
            Ad.coeffRef(i,j) = (d.coeffRef(i+1,j) - 2 * d.coeffRef(i,j) + d.coeffRef(i-1,j))* _dxx_inv +
                    (d.coeffRef(i,j+1) - 2 * d.coeffRef(i,j) + d.coeffRef(i,j-1))* _dyy_inv;
        }
    }
}
