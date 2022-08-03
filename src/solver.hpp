#ifndef __SOLVER_H_
#define __SOLVER_H_

#include <Eigen/Dense>

#include "enums.hpp"
#include "datastructures.hpp"

using namespace Eigen;

class Solver {
protected:
    Config& _config;
    matrix<cell_type>& _types;
public:
    virtual void solve(int& it, double& res) = 0;
    Solver(Config& config, matrix<cell_type>& types);
    void comp_residual(MatrixXd &x, MatrixXd &b, MatrixXd &res, double dx, double dy);
    void boundary_p(MatrixXd &P, matrix<cell_type>& type);
    void compute_l2Norm(double* residual, MatrixXd& residual_);
    void compute_average(double* residual, MatrixXd& residual_);

    static Solver* create(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types);
    static solver_type string2type(std::string& t);

};

#endif