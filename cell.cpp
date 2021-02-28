#include "cell.hpp"

Cell::Cell() {};

// Pressure Get and Set
double& Cell::pressure() {
    return *_pressure;
}

void Cell::set_pressure(double& value) {
    *_pressure = value;
}

void Cell::set_pressure_addr(double& addr){
    _pressure = &addr;
};

// Temperatur Get and Set
double& Cell::temperature() {
    return *_temperatur;
}

void Cell::set_temperature(double &value) {
    *_temperatur = value;
}

void Cell::set_temperature_addr(double& addr){
    _temperatur = &addr;
};

// Velocity Get and Set
double& Cell::velocity(velocity_type type) {
    double* vel;
    switch (type) {
        case velocity_type::U:
            vel = _vel_u;
            break;
        case velocity_type::V:
            vel = _vel_v;
            break;
    }
    return *vel;
};

// Velocity Get and Set
double& Cell::velocity_u() {
    return *_vel_u;
};

// Velocity Get and Set
double& Cell::velocity_v() {
    return *_vel_v;
};

void Cell::set_velocity(double& value, velocity_type type) {
    switch (type) {
        case velocity_type::U:
            *_vel_u = value;
            break;
        case velocity_type::V:
            *_vel_v = value;
            break;
    }
};

void Cell::set_velocity_u(double& value) {
    *_vel_u = value;
};

void Cell::set_velocity_addr_u(double& addr){
    _vel_u = &addr;
};

void Cell::set_velocity_addr_v(double& addr){
    _vel_v = &addr;
};


void Cell::set_velocity_v(double& value) {
    *_vel_v = value;
};

// borders Get and Set
bool& Cell::border(border_position position) {
    return _border[static_cast<int>(position)];
};

// Set border
void Cell::set_border(border_position position) {
    _border[static_cast<int>(position)] = true;
};

Cell* Cell::neighbour(neighbour_position position) {
    Cell* ret = nullptr;
    switch(position) {
        case neighbour_position::TOP:
            ret = _top;

        case neighbour_position::BOTTOM:
            ret = _bottom;

        case neighbour_position::LEFT:
            ret = _left;

        case neighbour_position::RIGHT:
            ret = _right;
    }

    return ret;
}

Cell* Cell::neighbour_top() {
    return _top;
}

Cell* Cell::neighbour_bottom() {
    return _bottom;
}

Cell* Cell::neighbour_left() {
    return _left;
}

Cell* Cell::neighbour_right() {
    return _right;
}

void Cell::set_neighbours(Cell *top, Cell *bottom, Cell *left, Cell *right) {
    _top = top;
    _bottom = bottom;
    _left = left;
    _right = right;
    if (*_type == cell_type::NOSLIP ) {
        int count = 0;
        if (_top != nullptr) count += _top->type() == cell_type::FLUID;
        if (_bottom != nullptr) count += _bottom->type() == cell_type::FLUID;
        if (_left != nullptr) count += _left->type() == cell_type::FLUID;
        if (_right != nullptr) count += _right->type() == cell_type::FLUID;
        if (count > 2) {
            throw "Obstical cell has more than 2 fluid neighbours!";
        }

        if((is_fluid(top) || is_inlet(top) || is_outlet(top)) && !is_fluid(left) && !is_fluid(right) && !is_fluid(bottom)) {
            _boundary = boundary_type::B_N;
        }

        if(!is_fluid(top) && !is_fluid(left) && is_fluid(right) && !is_fluid(bottom)) {
            _boundary = boundary_type::B_E;
        }

        if(!is_fluid(top) && is_fluid(left) && !is_fluid(right) && !is_fluid(bottom)) {
            _boundary = boundary_type::B_W;
        }

        if(!is_fluid(top) && !is_fluid(left) && !is_fluid(right) && (is_fluid(bottom) || is_inlet(bottom) || is_outlet(bottom))) {
            _boundary = boundary_type::B_S;
        }

        if(is_fluid(top) && !is_fluid(left) && is_fluid(right) && !is_fluid(bottom)) {
            _boundary = boundary_type::B_NE;
        }

        if(is_fluid(top) && is_fluid(left) && !is_fluid(right) && !is_fluid(bottom)) {
            _boundary = boundary_type::B_NW;
        }

        if(!is_fluid(top) && !is_fluid(left) && is_fluid(right) && is_fluid(bottom)) {
            _boundary = boundary_type::B_SE;
        }

        if(!is_fluid(top) && is_fluid(left) && !is_fluid(right) && is_fluid(bottom)) {
            _boundary = boundary_type::B_SW;
        }
    }
}

cell_type& Cell::type() {
    return *_type;
}

void Cell::set_type(cell_type type) {
    *_type = type;
}

bool Cell::is_type(cell_type type) {
    return *_type == type;
}

bool Cell::is_fluid() {
    return *_type == cell_type::FLUID ;
}

bool Cell::is_fluid(Cell *cell) {
    return cell != nullptr && cell->is_fluid();
}

bool Cell::is_inlet(Cell *cell) {
    return cell != nullptr && cell->is_type(cell_type::INLET);
}

bool Cell::is_outlet(Cell *cell) {
    return cell != nullptr && cell->is_type(cell_type::OUTLET);
}

boundary_type& Cell::boundary() {
    return _boundary;
}

void Cell::set_type_addr(cell_type& addr){
    _type = &addr;
};