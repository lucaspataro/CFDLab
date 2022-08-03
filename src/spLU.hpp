//
// Created by felix on 02.07.20.
//

#ifndef CFDLAB_SPLU_H
#define CFDLAB_SPLU_H


#include "solver.hpp"
#include <Eigen/Sparse>

class spLU: public Solver {
public:
    spLU(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types);
    void solve(int& it, double& res);

private:
    MatrixXd* _x = nullptr;
    MatrixXd* _rhs = nullptr;
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   _solver;
    SparseMatrix<double> _A;
};


#endif //CFDLAB_SPLU_H
