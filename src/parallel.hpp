#ifndef __PARALLEL_H_
#define __PARALLEL_H_

#include "mpi.h"
#include <vector>
#include "datastructures.hpp"
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

/**
 * master process splits up domain in as uniform as
 * possible splitted domains of [il-1,ir+1]×[jb-1,jt+1]
 * also make each proc aware of its neighbour, if there is no neighbour set to MPI_PROC_NULL
 * @param config configuration struct containing
 */
void init_parallel (Config& config);

/**
 * called within sor to communicate pressure
 * order :
 * send to left - receive from right
 * send to right - receive from left
 * send to the top - receive from the bottom
 * send to the bottom - receive from the top
 * @param P
 * @param il
 * @param ir
 * @param jb
 * @param jt
 * @param rank_l
 * @param rank_r
 * @param rank_b
 * @param rank_t
 * @param bufSend
 * @param bufRecv
 * @param status
 * @param chunk
 */

void boundary_comm(Config &config, MatrixXd &M,int height, int width, int tag);

/**
 * exchange the velocity values U and V between the processes: called in the end of calculate_uv
 * @param U
 * @param V
 * @param il
 * @param ir
 * @param jb
 * @param jt
 * @param rank_l
 * @param rank_r
 * @param rank_b
 * @param rank_t
 * @param bufSend
 * @param bufRecv
 * @param status
 * @param chunk
 */
void uv_comm (Config &config, MatrixXd  &U, MatrixXd &V, int height, int width);

void FG_comm (Config &config, MatrixXd  &F, MatrixXd &G, int height, int width);
#endif