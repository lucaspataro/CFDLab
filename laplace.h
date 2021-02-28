/*
 * ---------------------------------------------------------------------
 * Copyright (C) 2012, 2014 Tino Kluge (ttk448 at gmail.com)
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ---------------------------------------------------------------------
 */

#ifndef _laplace_h
#define _laplace_h

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <boost/multi_array.hpp>
#include <Eigen/Dense>
using namespace Eigen;

// 2d-array utility functions
// ---------------------------
namespace arr
{

// L^2 difference between two 2d-arrays
// allow_shift: removes a constant shift between 2 matrices,
//    i.e.  comparing X[][] and Y[][] + (X[1][1]-Y[1][1])
// ignore_corner: ignores all 4 corners, X[0][0], ...
double diff(MatrixXd& X,
            MatrixXd& Y,
            bool allow_shift, bool ignore_corner);

// L^2 difference between two vectors
double diff(MatrixXd& X,
            MatrixXd& Y,
            bool allow_shift=false);

// L^2 difference between vector and a constant
double diff(MatrixXd& X, double c);

} // namespace arr




// PDE functions (Laplace/Poisson)
// -------------------------------
namespace pde
{

namespace types
{
// boundary type
enum boundary {
   Dirichlet = 0,        // value fixed
   Neumann   = 1,        // first derivative fixed
};
} // namespace types


// sets the number of threads for the fftw solver
int fftw_threads(int n);
// cleans up memory used by fftw for threads
void fftw_clean();


// PDE functions to solve the discrete Poisson's Equation
// ------------------------------------------------------

// discrete Gradient operator on a uniform grid
// given u(x,y), it returns DX = a_1 du/dx, DY = a_2 du/dy,
// using one-sided finite differences
// same as above, but the input U is assumed to contain only inner grid points
// and the boundary grid points are inferred from the boundary condition
// U[i][j] = u(x_i, y_j), i=0...n-1
// DX[i][j] = a1*(U[i][j]-U[i-1][j])/h, i=0...n+1

// simplified version where only one boundary value can be specified
void grad(MatrixXd& DX, MatrixXd& DY,
          MatrixXd& U,
          double a1, double a2, double h1, double h2,
          double bdvalue, types::boundary bdcond);

// discrete Divergence operator on a uniform grid
// given u(x,y), v(x,y), it returns f = a_1 du/dx + a_2 dv/dx,
// using one-sided finite differences
void div(MatrixXd& F,
         MatrixXd& U,
         MatrixXd& V,
         double a1, double a2, double h1, double h2);

// discrete Laplace operator on a uniform grid
// the input U is assumed to contain only inner grid points and the
// boundary grid points are inferred from the boundary condition
// U[i][j] = u(x_i, y_j)
// F[i][j] = f(x_i, y_i)

// simplified version where only one boundary value can be specified
void laplace(MatrixXd& F,
             MatrixXd& U
             );

// Neumann boundary condition with a general right hand side F will not
// have a solution, however by looking at F we can say what the minimum
// L2-norm error will be for 0-Neumann boundary conditions:
// return = min ||Laplace U - F||
// min is taken over all U consistent with 0-Neumann condition
double neumann_error(const boost::multi_array<double,2>& F);

// given a right hand side F, we can find a constant Neumann-boundary
// value which will have a solution U, so that Laplace U = F,
// returns the boundary value
double neumann_compat(const boost::multi_array<double,2>& F,
                      double a1, double a2, double h1, double h2);




// solves the 2D Poisson equation: a1 u_xx + a2 u_yy = f
//
// discretised with uniform grid, x_{i+1} = x_i + h1, y_{i+1} = y_i + h2
// input: rhs       F[i][j] = f(x_i, y_j)
// output: solution  U[i][j] = u(x_i, y_j), ie only inner points by default

// Poisson solver, assuming uniform boundary
double poisolve(MatrixXd& U,
                MatrixXd& F
                );

} // namespace pde


#endif /* _laplace_h */
