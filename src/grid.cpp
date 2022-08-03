#include "grid.hpp"
#include "datastructures.hpp"
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>

Grid::Grid(Config& config, matrix<int> GeoArray, MatrixXd& U, MatrixXd& V, MatrixXd& P, MatrixXd& T, matrix<cell_type>& types):
    _boundary_size(config.boundary_size),
    _imax_b(config.imax + 2 * config.boundary_size), // +2 for the total size of the vector / matrix
    _jmax_b(config.jmax + 2 * config.boundary_size),
    _imax(config.imax),
    _jmax(config.jmax) {

    // Resizing the grid cells
    // the boundary size is given as a single value for all borders
    _cells.resize(_imax_b, std::vector<Cell>(_jmax_b, Cell()));

    // first set all the cell types
    for(auto j = 0; j< _jmax_b; j++) {
        for (auto i = 0; i < _imax_b; i++) {
            // set cell type
            cell_type type;
            type = lookupCellType(GeoArray[i][j]);
            types[i][j] = type;

            // set velocities
            U.coeffRef(i,j) = config.UI;
            V.coeffRef(i,j) = config.VI;

            // set Pressure
            P.coeffRef(i,j) = config.PI;

            // set temmperature
            T.coeffRef(i,j) = config.TI;

            _cells[i][j].set_type_addr(types[i][j]);
            _cells[i][j].set_velocity_addr_u(U.coeffRef(i,j));
            _cells[i][j].set_velocity_addr_v(V.coeffRef(i,j));
            _cells[i][j].set_pressure_addr(P.coeffRef(i,j));
            _cells[i][j].set_temperature_addr(T.coeffRef(i,j));

            // set init values to 0 if not fluid cell
            if (type != cell_type::FLUID) {
                double zeroVal = 0;
                _cells[i][j].set_pressure(zeroVal);
                _cells[i][j].set_temperature(zeroVal);
                _cells[i][j].set_velocity_u(zeroVal);
                _cells[i][j].set_velocity_v(zeroVal);
            }
        }
    }
    // then add all the neighbours
    for(auto j = 0; j < _jmax_b; j++) {
        for (auto i = 0; i < _imax_b; i++) {
            Cell* top = nullptr;
            Cell* bot = nullptr;
            Cell* left = nullptr;
            Cell* right = nullptr;
            if (i != 0) {
                left = &_cells[i-1][j];
            }
            if (i != _imax_b-1) {
                right = &_cells[i+1][j];
            }
            if (j != 0) {
                bot = &_cells[i][j-1];
            }
            if (j != _jmax_b-1) {
                top = &_cells[i][j+1];
            }

            _cells[i][j].set_neighbours(top, bot, left, right);
        }
    }
};

cell_type Grid::lookupCellType(int type_int) {
    cell_type type;
    switch(type_int) {
        case 1:
            type = cell_type::FLUID;
            break;
        case 2:
            type = cell_type::NOSLIP;
            break;
        case 3:
            type = cell_type::INLET;
            break;
        case 4:
            type = cell_type::OUTLET;
            break;
        case 5:
            type = cell_type::FREESLIP;
            break;
        case 6:
            type = cell_type::LID;
            break;
    }
    return type;
}

int Grid::jmaxb() const {
    return _jmax_b;
};

int Grid::imaxb() const {
    return _imax_b;
};

int Grid::imax() const {
    return _imax;
};

int Grid::jmax() const {
    return _jmax;
};

Cell& Grid::cell(int i, int j) {
    return _cells[i][j];
};

double Grid::max_u() {
    double max = 0.0;
    for(int x=0; x< Grid::imaxb();x++) {
        for (int y=0; y < Grid::jmaxb(); y++) {
            // Accessing velocity
            max = std::max(max, _cells[x][y].velocity_u());
        }
    }
    return max;
}

double Grid::max_v() {
    double max = 0.0;
    for(int x=0; x< Grid::imaxb();x++) {
        for (int y=0; y < Grid::jmaxb(); y++) {
            // Accessing velocity
            max = std::max(max, _cells[x][y].velocity_v());
        }
    }
    return max;
}

void Grid::set_temperature(MatrixXd& vec) {
    // Iterate over cells
    for(int y=0; y < Grid::jmaxb(); y++) {
        for (int x=0; x < Grid::imaxb();x++) {
            // Accessing velocity
            _cells[x][y].set_temperature(vec(x, y));
        }
    }
}

void Grid::print_pressure() {
    for(int y=Grid::jmaxb()-1;y >= 0;y--) {
        for(int x=0;x<Grid::imaxb();x++){
                // Accessing pressure
                std::cout<< Grid::_cells[x][y].pressure() << " ";
            }
            // Print new line
            std::cout << std::endl;
        }
}

void Grid::print_temperatur() {
    for(int y=Grid::jmaxb()-1;y >= 0;y--) {
        for(int x=0;x<Grid::imaxb();x++){
            // Accessing temperature
            std::cout<< Grid::_cells[x][y].temperature() << " ";
        }
        // Print new line
        std::cout << std::endl;
    }
}


