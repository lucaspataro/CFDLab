#ifndef __SOR_H_
#define __SOR_H_

#include <Eigen/Dense>

#include "datastructures.hpp"
#include "grid.hpp"
#include "solver.hpp"

using namespace Eigen;

class SOR: public Solver {
private:
    MatrixXd& _x;
    MatrixXd& _rhs;
    public:
        SOR(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types);
        void solve(int& it, double& res);
        void sor();
        void calculate_res(double& res);
};
#endif
