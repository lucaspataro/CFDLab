#include <stdio.h>
#include <iostream>

#include "boundary_val.hpp"

#include "datastructures.hpp"
#include "grid.hpp"
#include <Eigen/Dense>
using namespace Eigen;

inline void no_slip(Cell* cell, Cell* top, Cell* bottom, Cell* left, Cell* right) {
    double zero = 0.0, u = 0.0, v = 0.0;

    switch(cell -> boundary()) {
        case boundary_type::B_N:
            cell -> set_velocity_v(zero);
            u = - top -> velocity_u();
            cell -> set_velocity_u(u);
            break;

        case boundary_type::B_S:
            bottom -> set_velocity_v(zero);
            u = - bottom -> velocity_u();
            cell -> set_velocity_u(u);
            break;

        case boundary_type::B_E:
            cell -> set_velocity_u(zero);
            v = - right -> velocity_v();
            cell -> set_velocity_v(v);
            break;

        case boundary_type::B_W:
            left -> set_velocity_u(zero);
            v = - left -> velocity_v();
            cell -> set_velocity_v(v);
            break;

        case boundary_type::B_NE:
            cell -> set_velocity_u(zero);
            cell -> set_velocity_v(zero);
            break;

        case boundary_type::B_NW:
            left -> set_velocity_u(zero);
            cell -> set_velocity_v(zero);
            u = - top -> velocity_u();
            cell -> set_velocity_u(u);
            break;
            
        case boundary_type::B_SE:
            cell -> set_velocity_u(zero);
            bottom -> set_velocity_v(zero);

            //u = - (left -> neighbour_bottom() -> velocity_u());
            v = - right -> velocity_v();

            //bottom -> set_velocity_u(u);
            cell -> set_velocity_v(v);
            break;

        case boundary_type::B_SW:
            left -> set_velocity_u(zero);
            bottom -> set_velocity_v(zero);

            u = - bottom -> velocity_u();
            v = - left -> velocity_v();

            cell -> set_velocity_u(u);
            cell -> set_velocity_v(v);
            break;
    }
}

inline void free_slip(Cell* cell, Cell* top, Cell* bottom, Cell* left, Cell* right) {
    double zero = 0.0, u = 0.0, v = 0.0;

    switch(cell -> boundary()) {
        case boundary_type::B_N:
            cell -> set_velocity_v(zero);

            u = top -> velocity_u();
            cell -> set_velocity_u(u);

            u = left -> neighbour_top()->velocity_u();
            left -> set_velocity_u(u);
            break;

        case boundary_type::B_S:
            bottom -> set_velocity_v(zero);

            u = left -> neighbour_bottom()->velocity_u();
            left -> set_velocity_u(u);

            u = bottom -> velocity_u();
            cell -> set_velocity_u(u);
            break;

        case boundary_type::B_E:
            cell -> set_velocity_u(zero);
            
            v = right -> velocity_v();
            cell -> set_velocity_v(v);

            v = right -> neighbour_bottom()->velocity_v();
            bottom -> set_velocity_v(v);

            break;

        case boundary_type::B_W:
            left -> set_velocity_u(zero);

            v = left -> neighbour_bottom()->velocity_v();
            bottom -> set_velocity_v(v);

            v = left -> velocity_v();
            cell -> set_velocity_v(v);
            break;

        case boundary_type::B_NE:
            cell -> set_velocity_u(zero);
            cell -> set_velocity_v(zero);

            u = left -> neighbour_top()->velocity_u();
            v = right -> neighbour_bottom()->velocity_v();
            left -> set_velocity_u(u);
            bottom -> set_velocity_v(v);
            break;

        case boundary_type::B_NW:
            cell -> set_velocity_v(zero);
            left -> set_velocity_u(zero);
            
            u = top -> velocity_u();
            v = left -> neighbour_bottom()->velocity_v();

            cell -> set_velocity_u(u);
            bottom -> set_velocity_v(v);
            break;
            
        case boundary_type::B_SE:
            cell -> set_velocity_u(zero);
            bottom -> set_velocity_v(zero);

            u = left -> neighbour_bottom()->velocity_u();
            v = right -> velocity_v();

            left -> set_velocity_u(u);
            cell -> set_velocity_v(v);
            break;

        case boundary_type::B_SW:
            left -> set_velocity_u(zero);
            bottom -> set_velocity_v(zero);

            u = bottom -> velocity_u();
            v = left -> velocity_v();

            cell -> set_velocity_u(u);
            cell -> set_velocity_v(v);
            break;
    }
}

inline void boundary(Cell& cell) {
  if (!cell.is_fluid()) {
      cell_type type = cell.type();
      Cell* top = cell.neighbour_top();
      Cell* bottom = cell.neighbour_bottom();
      Cell* left = cell.neighbour_left();
      Cell* right = cell.neighbour_right();
      
      double zero = 0.0;
      double one = 1.0;
      double vel = 0.0;

      switch(type) {
          case cell_type::FREESLIP:
              free_slip(&cell, top, bottom, left, right);
              break;

          case cell_type::NOSLIP:
              no_slip(&cell, top, bottom, left, right);
              break;

          case cell_type::INLET:
              // for inlet set initial values
              cell.set_velocity_u(one);
              cell.set_velocity_v(zero);
              break;

          case cell_type::OUTLET:
              // for outlet use value of neighbour
              if(top != nullptr && top->is_fluid()) {
                  cell.set_velocity_v(top->velocity_v());
              }
              if(bottom != nullptr && bottom->is_fluid()) {
                  cell.set_velocity_v(bottom->velocity_v());
              }

              if(left != nullptr && left->is_fluid()) {
                  cell.set_velocity_u(left -> velocity_u());
              } else if ( right != nullptr && right->is_fluid()) {
                  cell.set_velocity_u(right->velocity_u());
              }
              break;

          case cell_type::LID:
              // special lid boundary for cavity scenario
              if(bottom != nullptr) {
                  cell.set_velocity_v(zero);
                  vel = 2.0 - (bottom->velocity_u());
                  cell.set_velocity_u(vel);
              }
              break;
      }
  }
}

void boundaryvalues(Grid& grid, Config& config) {
  for (auto j = 0; j < config.jmax + 2; j++) {
      for (auto i = 0; i < config.imax + 2; i++) {
          boundary(grid.cell(i, j));
      }
  }
}

void boundaryvalues_t(Grid &grid, Config& config, MatrixXd& T) {
    
    for(int i = 0; i<config.imax+2; ++i){
        for(int j = 0; j<config.jmax+2; ++j){
            boundary_type type = grid.cell(i,j).boundary();
            switch(type) {
                case boundary_type::B_E:
                    T(i,j) = T(i+1,j);
                    break;
                case boundary_type::B_W:
                    T(i,j) = T(i-1,j);
                    break;
                case boundary_type::B_N:
                    T(i,j) = T(i,j+1);
                    break;
                case boundary_type::B_S:
                    T(i,j) = T(i,j-1);
                    break;
                case boundary_type::B_NE:
                    T(i,j) = (T(i,j+1) + T(i+1,j))/2;
                    break;
                case boundary_type::B_NW:
                    T(i,j) = (T(i,j+1) + T(i-1,j))/2;
                    break;
                case boundary_type::B_SE:
                    T(i,j) = (T(i,j-1) + T(i+1,j))/2;
                    break;
                case boundary_type::B_SW:
                    T(i,j) = (T(i,j-1) + T(i-1,j))/2;
                    break;
            }
        }
    }

    if(config.problem == "NaturalConvection" || config.problem == "FluidTrap") {
        if (config.omg_i == 1) {
            for (int j=0; j<config.jmax+2; j++) {
                T(0,j) = 2*config.T_h - T(1,j);
            }
        }
        if (config.omg_i == config.iproc) {
            for (int j=0; j<config.jmax+2; j++) {
                T(config.imax+1,j) = 2*config.T_c - T(config.imax,j);
            }
        }
    }


    if(config.problem == "FluidTrapReversed") {
        if (config.omg_i == 1) {
            for (int j=0; j<config.jmax+2; j++) {
                T(0,j) = 2*config.T_c - T(1,j);
            }
        }
        if (config.omg_i == config.iproc) {
            for (int j=0; j<config.jmax+2; j++) {
                T(config.imax+1, j) = 2*config.T_h - T(config.imax,j);
            }
        }
    }

    if(config.problem == "RayleighBenardConvection") {
        if (config.omg_j == 1) {
            for(int i=0; i<config.imax+2; i++) {
                T(i,0) = 2*config.T_h - T(i, 1);
            }
        }
        if (config.omg_j == config.jproc) {
            for(int i=0; i<config.imax+2; i++) {
                T(i, config.jmax+1) = 2*config.T_c - T(i,config.jmax);
            }
        }
    }
}
