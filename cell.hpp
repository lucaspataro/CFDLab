#ifndef CFDLAB_CELL_HPP
#define CFDLAB_CELL_HPP
#include <array>
#include "enums.hpp"


class Cell {
public:
    // Constructors
    Cell();
    Cell(double& PI, double& UI, double& VI, double& TI);

    // Get + Set pressure
    double& pressure();
    void set_pressure(double& value);
    void set_pressure_addr(double& addr);

    // Get + Set temperature
    double& temperature();
    void set_temperature(double& value);
    void set_temperature_addr(double& addr);

    // Get + Set velocity
    double& velocity(velocity_type type);
    double& velocity_u();
    double& velocity_v();
    void set_velocity(double& value, velocity_type type);
    void set_velocity_u(double& value);
    void set_velocity_v(double& value);
    void set_velocity_addr_u(double& addr);
    void set_velocity_addr_v(double& addr);

    // Get + Set border
    bool& border(border_position position);
    void set_border(border_position position);

    // Get + Set neighbours
    Cell* neighbour(neighbour_position position);
    Cell* neighbour_top();
    Cell* neighbour_bottom();
    Cell* neighbour_left();
    Cell* neighbour_right();
    void set_neighbours(Cell *top, Cell *bottom, Cell *left, Cell *right);

    // Get + Set cell type
    cell_type& type();
    void set_type(cell_type type);
    bool is_type(cell_type type);
    bool is_fluid();
    void set_type_addr(cell_type& addr);

    // Getter boundary
    boundary_type& boundary();


private:
    // one pressure value per call
    double* _pressure = nullptr;

    // one Temperatur value per cell
    double* _temperatur = nullptr;

    // Fixed size velocity
    double* _vel_u = nullptr;
    double* _vel_v = nullptr;
    
    // Fixed number of borders
    std::array<bool, 4> _border = {false};

    // the type of cell: FLUID, OBSTACLE, ...
    cell_type* _type = nullptr;

    boundary_type _boundary = boundary_type::UNSET;

    // Neighbours
    Cell* _top;
    Cell* _bottom;
    Cell* _left;
    Cell* _right;

    bool is_fluid(Cell* cell);
    bool is_inlet(Cell* cell);
    bool is_outlet(Cell* cell);
};

#endif //CFDLAB_CELL_HPP
