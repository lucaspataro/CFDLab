//
// Created by felix on 27.06.20.
//

#ifndef CFDLAB_CG_H
#define CFDLAB_CG_H


#include "solver.hpp"

class CG: public Solver {
public:
    CG(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types);
    void solve(int& it, double& res);
private:
    MatrixXd* _d = nullptr;
    MatrixXd* _Ad = nullptr;
    MatrixXd* _res = nullptr;
    MatrixXd* _x = nullptr;
    MatrixXd* _rhs = nullptr;
    double _dyy_inv;
    double _dxx_inv;
    void calc_Ad();
};


#endif //CFDLAB_CG_H
