#include "uvp.hpp"

#include <bits/stdc++.h>
#include "mpi.h"

#include <cmath>
#include <cstdlib>
#include <vector>

#include "datastructures.hpp"
#include "helper.hpp"
#include <Eigen/Dense>
using namespace Eigen;

using namespace std;

static bool abs_compare(double a, double b) {
  return (abs(a) < abs(b));
}

// Determines the value of F and G
void calculate_fg(Grid &grid,
                  Config& config,
                  MatrixXd &F, MatrixXd &G,
                  MatrixXd &U, MatrixXd &V,
                  MatrixXd &T, matrix<cell_type> &type) {

  double temp, temp2;
  double d2udx2, d2udy2, du2dx, duvdy;
  double d2vdx2, d2vdy2, dv2dy, duvdx;
  double u_ij, u_ip1j, u_im1j, u_ijp1, u_ijm1;
  double v_ij, v_ip1j, v_im1j, v_ijp1, v_ijm1;
  double v_ip1jm1, u_im1jp1;
  double T_ij, T_ip1j, T_ijp1;

  for (auto j = 1; j < config.jmax + 1; j++) {
    for (auto i = 1; i < config.imax + 1; i++) {
        // check if its a fluid if it is not it has to be a noslip type as it is inside the domain
        if (type[i][j] == cell_type::FLUID && type[i+1][j] == cell_type::FLUID) {
            u_ij = U.coeffRef(i,j);
            u_ip1j = U.coeffRef(i + 1,j);
            u_im1j = U.coeffRef(i - 1,j);
            u_ijp1 = U.coeffRef(i,j + 1);
            u_ijm1 = U.coeffRef(i,j - 1);

            v_ij = V.coeffRef(i,j);
            v_ip1j = V.coeffRef(i + 1,j);
            v_im1j = V.coeffRef(i - 1,j);
            v_ijp1 = V.coeffRef(i,j + 1);
            v_ijm1 = V.coeffRef(i,j - 1);

            u_im1jp1 = U.coeffRef(i - 1,j + 1);
            v_ip1jm1 = V.coeffRef(i + 1,j - 1);

            T_ij = T(i,j);
            T_ip1j = T(i+1,j);
            T_ijp1 = T(i,j+1);

            // calculate second derivative of u after x
            d2udx2 = (u_ip1j - 2 * u_ij + u_im1j) / (config.dx * config.dx);

            // calculate second derivative of u after y
            d2udy2 = (u_ijp1 - 2 * u_ij + u_ijm1) / (config.dy * config.dy);

            // calculate first derivative of u² after x
            du2dx= pow(u_ij + u_ip1j, 2)-pow(u_im1j + u_ij, 2);

            temp = abs(u_ij + u_ip1j) * (u_ij - u_ip1j);
            temp2 = abs(u_im1j + u_ij) * (u_im1j - u_ij);

            du2dx += (config.alpha * (temp - temp2));
            du2dx /= (4 * config.dx);

            // calculate first derivative of uv after y
            temp = (v_ij + v_ip1j) * (u_ij + u_ijp1);
            temp2 = (v_ijm1 + v_ip1jm1) * (u_ijm1 + u_ij);
            duvdy = temp - temp2;

            temp = abs(v_ij + v_ip1j) * (u_ij - u_ijp1);
            temp2 = abs(v_ijm1 + v_ip1jm1) * (u_ijm1 - u_ij);

            duvdy += (config.alpha * (temp - temp2));
            duvdy /= (4 * config.dy);

            // calcute F
            F.coeffRef(i,j) = u_ij + config.dt * ((d2udx2 + d2udy2) / config.Re - du2dx - duvdy);
            if (config.calcTemp) {
                F.coeffRef(i,j) -= config.beta * 0.5 * config.dt * (T_ij + T_ip1j) * config.GX;
            }
        }
    }
  }

  for (auto i = 1; i < config.imax + 1; i++) {
    for (auto j = 1; j < config.jmax + 1; j++) {
        if (type[i][j] == cell_type::FLUID && type[i][j+1] == cell_type::FLUID) {
            u_ij = U.coeffRef(i,j);
            u_im1j = U.coeffRef(i - 1,j);
            u_ijp1 = U.coeffRef(i,j + 1);

            v_ij = V.coeffRef(i,j);
            v_ip1j = V.coeffRef(i + 1,j);
            v_im1j = V.coeffRef(i - 1,j);
            v_ijp1 = V.coeffRef(i,j + 1);
            v_ijm1 = V.coeffRef(i,j - 1);

            u_im1jp1 = U.coeffRef(i - 1,j + 1);
            v_ip1jm1 = V.coeffRef(i + 1,j - 1);

            T_ij = T(i,j);
            T_ip1j = T(i+1,j);
            T_ijp1 = T(i,j+1);

            // calculate second derivative of v after x
            d2vdx2 = (v_ip1j - 2 * v_ij + v_im1j) / (config.dx * config.dx);

            // calculate second derivative of v after y
            d2vdy2 = (v_ijp1 - 2 * v_ij + v_ijm1) / (config.dy * config.dy);

            // calculate first derivative of v² after y
            temp = pow(v_ij + v_ijp1, 2);
            temp2 = pow(v_ijm1 + v_ij, 2);
            dv2dy = temp - temp2;

            temp = abs(v_ij + v_ijp1) * (v_ij - v_ijp1);
            temp2 = abs(v_ijm1 + v_ij) * (v_ijm1 - v_ij);

            dv2dy += (config.alpha * (temp - temp2));
            dv2dy /= (4 * config.dy);

            // calculate first derivative of uv after x
            temp = (u_ij + u_ijp1) * (v_ij + v_ip1j);
            temp2 = (u_im1j + u_im1jp1) * (v_im1j + v_ij);
            duvdx = temp - temp2;

            temp = abs(u_ij + u_ijp1) * (v_ij - v_ip1j);
            temp2 = abs(u_im1j + u_im1jp1) * (v_im1j - v_ij);

            duvdx += (config.alpha * (temp - temp2));
            duvdx /= (4 * config.dx);

            // calcute G
            G.coeffRef(i,j) = v_ij + config.dt * ((d2vdx2 + d2vdy2) / config.Re - dv2dy - duvdx);
            if (config.calcTemp) {
                G.coeffRef(i,j) -= config.beta * 0.5 * config.dt * (T_ij + T_ijp1) * config.GY;
            }
        }
    }
  }

  for (auto i = 0; i < config.imax+2; i++) {
      for (auto j = 0; j < config.jmax+2; j++) {
          cell_type type_ij = type[i][j];
          if ( type_ij != cell_type::FLUID) {
              Cell* top = nullptr;
              Cell* bottom = nullptr;
              Cell* left = nullptr;
              Cell* right = nullptr;
              cell_type temp_type;
              if(grid.cell(i,j).neighbour_top()) {
                  temp_type = grid.cell(i,j).neighbour_top()->type();
                  if (temp_type == cell_type::FLUID || temp_type == cell_type::INLET || temp_type == cell_type::OUTLET) {
                      top = grid.cell(i,j).neighbour_top();
                  }
              }
              if(grid.cell(i,j).neighbour_bottom()) {
                  temp_type = grid.cell(i,j).neighbour_bottom()->type();
                  if (temp_type == cell_type::FLUID || temp_type == cell_type::INLET || temp_type == cell_type::OUTLET) {
                      bottom = grid.cell(i,j).neighbour_bottom();
                  }
              }
              if(grid.cell(i,j).neighbour_left()) {
                  temp_type = grid.cell(i,j).neighbour_left()->type();
                  if (temp_type == cell_type::FLUID || temp_type == cell_type::INLET || temp_type == cell_type::OUTLET) {
                      left = grid.cell(i,j).neighbour_left();
                  }
              }
              if(grid.cell(i,j).neighbour_right()) {
                  temp_type = grid.cell(i,j).neighbour_right()->type();
                  if (temp_type == cell_type::FLUID || temp_type == cell_type::INLET || temp_type == cell_type::OUTLET) {
                      right = grid.cell(i,j).neighbour_right();
                  }
              }
              if (type_ij == cell_type::FREESLIP || type_ij == cell_type::NOSLIP || type_ij == cell_type::LID) {
                  if (top != nullptr) {
                      G.coeffRef(i,j) = V.coeffRef(i,j);
                  }
                  if (bottom!= nullptr) {
                      G.coeffRef(i,j-1) = V.coeffRef(i,j-1);
                  }
                  if (left != nullptr) {
                      F.coeffRef(i-1,j) = U.coeffRef(i-1,j);
                  }
                  if (right != nullptr) {
                      F.coeffRef(i,j) = U.coeffRef(i,j);
                  }
              }
              if (type_ij == cell_type::INLET || type_ij == cell_type::OUTLET) {
                  if (top != nullptr && top->is_fluid()) {
                      G.coeffRef(i,j) = V.coeffRef(i,j);
                  }
                  if (bottom!= nullptr && bottom->is_fluid()) {
                      G.coeffRef(i,j-1) = V.coeffRef(i,j-1);
                  }
                  if (left != nullptr && left->is_fluid()) {
                      F.coeffRef(i-1,j) = U.coeffRef(i-1,j);
                  }
                  if (right != nullptr && right->is_fluid()) {
                      F.coeffRef(i,j) = U.coeffRef(i,j);
                  }
              }
          }
      }
  }
}

// This operation computes the right hand side of the pressure poisson equation.
void calculate_rs(Config& config, MatrixXd &F, MatrixXd &G, MatrixXd &RS) {
  for (auto i = 1; i < config.imax + 1; i++) {
    for (auto j = 1; j < config.jmax + 1; j++) {
      RS(i,j) = (((F.coeffRef(i,j) - F.coeffRef(i - 1,j)) / config.dx) + ((G.coeffRef(i,j) - G.coeffRef(i,j - 1)) / config.dy)) * (1/config.dt);
    }
  }
}

// Determines the maximal time step size
void calculate_dt(Grid &grid, Config& config) {

  if (config.tau <= 0) {
    return;
  }

  // get max velocity
  double uMax = grid.max_u();
  double vMax = grid.max_v();
  double globuMax, globvMax;

  MPI_Allreduce(&uMax, &globuMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&vMax, &globvMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

  double dx2 = config.dx * config.dx;
  double dy2 = config.dy * config.dy;
  double dt1 = (config.Re / 2) / (1 / (dx2) + 1 / (dy2));
  double dt2 = config.dx / globuMax;
  double dt3 = config.dy / globvMax;
  double dt4 = std::numeric_limits<double>::max();

  if(config.calcTemp){
      dt4 = ( (config.Re * config.PR) / 2.0 ) * ( (dx2 * dy2) / (dx2 + dy2) );
  }

  config.dt = config.tau * min(min(dt1, dt2), min(dt3, dt4));
}

void calculate_uv(Grid &grid, Config& config, MatrixXd &F, MatrixXd &G) {
  
  double u = 0, v = 0;
  for (auto i = 1; i < config.imax + 1; i++) {
    for (auto j = 1; j < config.jmax + 1; j++) {
        if (grid.cell(i,j).is_fluid()) {
            if (grid.cell(i+1,j).is_fluid() || grid.cell(i+1,j).is_type(cell_type::OUTLET)) {
                u = F.coeffRef(i,j) - (config.dt / config.dx) *
                                    (grid.cell(i + 1, j).pressure() - grid.cell(i, j).pressure());
                grid.cell(i, j).set_velocity_u(u);
            }
        }
    }
  }
  
  // calc v for imax
  for (auto i = 1; i < config.imax + 1; i++) {
      for (auto j = 1; j < config.jmax + 1; j++) {
          if (grid.cell(i,j).is_fluid()) {
              if (grid.cell(i, j + 1).is_fluid()) {
                  v = G.coeffRef(i,j) - (config.dt / config.dy) *
                                      (grid.cell(i, j + 1).pressure() - grid.cell(i, j).pressure());
                  grid.cell(i, j).set_velocity_v(v);
              }
          }
      }
  }
}