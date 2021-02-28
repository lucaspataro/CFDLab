#include "vector"
#include <string>
#include "enums.hpp"

#ifndef DATA_STRUCTURES_H
#define DATA_STRUCTURES_H

template<typename T>
using matrix = std::vector<std::vector<T>>;

struct Config {
    int xlength;
    int ylength;
    int Re;
    int itermax;
    int imax;
    int jmax;
    int l_imax;
    int l_jmax;
    int calcTemp = 0;
    int iproc = 1;
    int jproc = 1;
    int levels;
    int n_fluid;
    int total_n_fluid;
    int boundary_size;

    int num_proc;
    int rank;
    int rank_l;
    int rank_r;
    int rank_t;
    int rank_b;

    int il;
    int ir;
    int jb;
    int jt;

    int omg_i;
    int omg_j;

    double t_end;
    double dt;
    double omg;
    double eps;
    double tau;
    double alpha;
    double dt_value;
    double UI;
    double VI;
    double GX;
    double GY;
    double PI;
    double PR;
    double TI;
    double T_h;
    double T_c;
    double beta;
    double dx;
    double dy;

    solver_type solver;
    std::string problem;
    std::string geometry;
};

#endif //DATA_STRUCTURES_H
