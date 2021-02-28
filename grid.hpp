

#ifndef CFDLAB_GRID_H
#define CFDLAB_GRID_H
#include <vector>
#include "cell.hpp"
#include "datastructures.hpp"
#include <Eigen/Dense>
using namespace Eigen;
/*
 * CamelCase for function that returns or sets Objects
 * snake_case for functions that return a specific attribute like imax or jmax
 */

class Grid {
public:
    Grid(Config& config, matrix<int> GeoArray, MatrixXd& U, MatrixXd& V, MatrixXd& P, MatrixXd& T, matrix<cell_type>& type);

    double max_u();
    double max_v();

    // set Temperatur
    void set_temperature(MatrixXd& vec);
    
    // specific row and column
    Cell& cell(int i, int j);


    // Get Dimensions
    int imax() const;
    int jmax() const;

    // i/jmax with borders
    int imaxb() const;
    int jmaxb() const;

    // Print matrices
    void print_pressure();
    void print_temperatur();

private:
    cell_type lookupCellType(int type_int);
    matrix<Cell> _cells;

    const int _imax;
    const int _jmax;
    const int _imax_b;
    const int _jmax_b;
    const int _boundary_size;

};

#endif //CFDLAB_GRID_H
