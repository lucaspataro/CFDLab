#include "parallel.hpp"
#include "datastructures.hpp"
#include "enums.hpp"
#include <vector>

using namespace std;

void init_parallel (Config& config) {
  
  if(config.iproc == 1 && config.jproc == 1){
    int ndim = 2;
    int dimensions[ndim];

    dimensions[0] = dimensions[1] = 0;
    MPI_Dims_create(config.num_proc, ndim, dimensions);

    config.iproc = dimensions[0];
    config.jproc = dimensions[1];
  }

  if (config.rank == 0) {
      printf("Number of Processors: %d: Grid Dimensions = [%d x %d] \n",
             config.num_proc, config.iproc, config.jproc);
  }

  int rows = config.jmax / config.jproc;
  int cols = config.imax / config.iproc;

// calculate domain indices
  config.omg_i = (config.rank % config.iproc) + 1;
  config.omg_j = (config.rank / config.iproc) + 1;

  // total leftover
  int leftover_it = config.imax % config.iproc;
  int leftover_jt = config.jmax % config.jproc;

  // leftovers applied to each process
  int leftover_i = leftover_it >= config.omg_i ? 1 : 0;
  int leftover_j = leftover_jt >= config.omg_j ? 1 : 0;


  // calculate domain bounds
  config.il = (config.omg_i - 1) * cols;
  config.il += leftover_it > (config.omg_i-1) ? (config.omg_i-1): leftover_it;
  config.il += 1;

  config.ir = config.omg_i * cols;
  config.ir += leftover_i == 1 ? config.omg_i: leftover_it;

  config.jb = (config.omg_j - 1) * rows;
  config.jb += leftover_jt > (config.omg_j-1) ? (config.omg_j-1): leftover_jt;
  config.jb += 1;

  config.jt = config.omg_j * rows;
  config.jt += leftover_j == 1 ? config.omg_j: leftover_jt;

  // set processes on right and top to end of domain
  if (config.omg_i == config.iproc) {
    config.ir = config.imax;
  }

  if (config.omg_j == config.jproc) {
    config.jt = config.jmax;
  }

  // assign neighbour process ranks
  config.rank_t = MPI_PROC_NULL;
  config.rank_l = MPI_PROC_NULL;
  config.rank_r = MPI_PROC_NULL;
  config.rank_b = MPI_PROC_NULL;

  if (config.il > 1) {
    config.rank_l = config.rank - 1;
  }

  if (config.ir < config.imax) {
      config.rank_r = config.rank + 1;
  }

  if (config.jb > 1) {
      config.rank_b = config.rank - config.iproc;
  }

  if (config.jt < config.jmax) {
      config.rank_t = config.rank + config.iproc;
  }

}

void store_in_matrix(MatrixXd &M, double* vals, int width, int height, neighbour_position pos) {
    switch (pos) {
        case neighbour_position::LEFT:
            for (auto i = 0; i < height; i++) {
                M.coeffRef(0,i) = vals[i];
            }
            break;
        case neighbour_position::RIGHT:
            for (auto i = 0; i < height; i++) {
                M.coeffRef(width-1,i) = vals[i];
            }
            break;
        case neighbour_position::TOP:
            for (auto i = 0; i < width; i++) {
                M.coeffRef(i,height-1) = vals[i];
            }
            break;
        case neighbour_position::BOTTOM:
            for (auto i = 0; i < width; i++) {
                M.coeffRef(i,0) = vals[i];
            }
            break;
    }

}

void get_from_matrix(MatrixXd &M, double* vals, int width, int height, neighbour_position pos) {
    switch (pos) {
        case neighbour_position::LEFT:
            for (auto i = 0; i < height; i++) {
                vals[i] = M.coeffRef(1,i);
            }
            break;
        case neighbour_position::RIGHT:
            for (auto i = 0; i < height; i++) {
                vals[i] = M.coeffRef(width-2,i);
            }
            break;
        case neighbour_position::TOP:
            for (auto i = 0; i < width; i++) {
                vals[i] = M.coeffRef(i,height-2);
            }
            break;
        case neighbour_position::BOTTOM:
            for (auto i = 0; i < width; i++) {
                vals[i] = M.coeffRef(i,1);
            }
            break;
    }

}

void boundary_comm(Config& config, MatrixXd &M, int width, int height, int tag) {
    double* verticalVecSend = new double[height];
    double* verticalVecRec = new double[height];
    double* horVecSend = new double[width];
    double* horVecRec = new double[width];

    // communicate in i-direction (x-axis)
    if (config.omg_i % 2) {
        // omg_i = 1, 3, 5, 7, ..., (2k + 1 <= iproc)
        // receive - send
        if (config.omg_i != 1) {
            // omg_i = 3, 5, 7, ..., (2k + 1 <= iproc)
            get_from_matrix(M, verticalVecSend, width, height, neighbour_position::LEFT);
            MPI_Recv(verticalVecRec, height, MPI_DOUBLE, config.rank_l, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(verticalVecSend, height, MPI_DOUBLE, config.rank_l, tag, MPI_COMM_WORLD);
            store_in_matrix(M, verticalVecRec, width, height, neighbour_position::LEFT);
        }
        if (config.omg_i != config.iproc) {
            // omg_i = 1, 3, 5, 7, ..., (2k + 1 < iproc)
            get_from_matrix(M, verticalVecSend, width, height, neighbour_position::RIGHT);
            MPI_Recv(verticalVecRec, height, MPI_DOUBLE, config.rank_r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(verticalVecSend, height, MPI_DOUBLE, config.rank_r, tag, MPI_COMM_WORLD);
            store_in_matrix(M, verticalVecRec, width, height, neighbour_position::RIGHT);
        }
    } else {
        // omg_i = 2, 4, 6, ..., (2k <= iproc)
        // receive - send
        if (config.omg_i != config.iproc) {
            // omg_i = 2, 4, 6, ..., (2k < iproc)
            get_from_matrix(M, verticalVecSend, width, height, neighbour_position::RIGHT);
            MPI_Send(verticalVecSend, height, MPI_DOUBLE, config.rank_r, tag, MPI_COMM_WORLD);
            MPI_Recv(verticalVecRec, height, MPI_DOUBLE, config.rank_r, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            store_in_matrix(M, verticalVecRec, width, height, neighbour_position::RIGHT);
        }

        // omg_i = 2, 4, 6, ..., (2k <= iproc)
        get_from_matrix(M, verticalVecSend, width, height, neighbour_position::LEFT);
        MPI_Send(verticalVecSend, height, MPI_DOUBLE, config.rank_l, tag, MPI_COMM_WORLD);
        MPI_Recv(verticalVecRec, height, MPI_DOUBLE, config.rank_l, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        store_in_matrix(M, verticalVecRec, width, height, neighbour_position::LEFT);
    }

    // communicate in j-direction (y-axis)
    if (config.omg_j % 2) {
        // omg_j = 1, 3, 5, 7, ..., (2k + 1 <= jproc)
        if (config.omg_j != 1) {
            // omg_j = 3, 5, 7, ..., (2k + 1 <= jproc)
            get_from_matrix(M, horVecSend, width, height, neighbour_position::BOTTOM);
            MPI_Recv(horVecRec, width, MPI_DOUBLE, config.rank_b, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(horVecSend, width, MPI_DOUBLE, config.rank_b, tag, MPI_COMM_WORLD);
            store_in_matrix(M, horVecRec, width, height, neighbour_position::BOTTOM);
        }
        if (config.omg_j != config.jproc) {
            // omg_j = 1, 3, 5, 7, ..., (2k + 1 < jproc)
            get_from_matrix(M, horVecSend, width, height, neighbour_position::TOP);
            MPI_Recv(horVecRec, width, MPI_DOUBLE, config.rank_t, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(horVecSend, width, MPI_DOUBLE, config.rank_t, tag,MPI_COMM_WORLD);
            store_in_matrix(M, horVecRec, width, height, neighbour_position::TOP);
        }
    } else {
        // omg_j = 2, 4, 6, ..., (2k <= jproc)
        if (config.omg_j != config.jproc) {
            // omg_j = 2, 4, 6, ..., (2k < jproc)
            get_from_matrix(M, horVecSend, width, height, neighbour_position::TOP);
            MPI_Send(horVecSend, width, MPI_DOUBLE, config.rank_t, tag,MPI_COMM_WORLD);
            MPI_Recv(horVecRec, width, MPI_DOUBLE, config.rank_t, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            store_in_matrix(M, horVecRec, width, height, neighbour_position::TOP);
        }

        // omg_j = 2, 4, 6, ..., (2k <= jproc)
        get_from_matrix(M, horVecSend, width, height, neighbour_position::BOTTOM);
        MPI_Send(horVecSend, width, MPI_DOUBLE, config.rank_b, tag,MPI_COMM_WORLD);
        MPI_Recv(horVecRec, width, MPI_DOUBLE, config.rank_b, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        store_in_matrix(M, horVecRec, width, height, neighbour_position::BOTTOM);
    }


    delete[] verticalVecSend;
    delete[] verticalVecRec;
    delete[] horVecSend;
    delete[] horVecRec;
}

void store_vel_in_matrix(MatrixXd &U, MatrixXd &V, double* vals, int width, int height, neighbour_position pos) {
    switch (pos) {
        case neighbour_position::LEFT:
            for (auto i = 0; i < height; i++) {
                U.coeffRef(0,i) = vals[2*i];
                V.coeffRef(0,i) = vals[2*i + 1];
            }
            break;
        case neighbour_position::RIGHT:
            for (auto i = 0; i < height; i++) {
                U.coeffRef(width-1,i) = vals[2*i];
                V.coeffRef(width-1,i) = vals[2*i+1];
            }
            break;
        case neighbour_position::TOP:
            for (auto i = 0; i < width; i++) {
                U.coeffRef(i,height-1) = vals[2*i];
                V.coeffRef(i,height-1) = vals[2*i+1];
            }
            break;
        case neighbour_position::BOTTOM:
            for (auto i = 0; i < width; i++) {
                U.coeffRef(i,0) = vals[2*i];
                V.coeffRef(i,0) = vals[2*i+1];
            }
            break;
    }

}

void get_vel_from_matrix(MatrixXd &U, MatrixXd &V, double* vals, int width, int height, neighbour_position pos) {
    switch (pos) {
        case neighbour_position::LEFT:
            for (auto i = 0; i < height; i++) {
                vals[2*i] = U.coeffRef(1,i);
                vals[2*i+1] = V.coeffRef(1,i);
            }
            break;
        case neighbour_position::RIGHT:
            for (auto i = 0; i < height; i++) {
                vals[2*i] = U.coeffRef(width-2,i);
                vals[2*i+1] = V.coeffRef(width-2,i);
            }
            break;
        case neighbour_position::TOP:
            for (auto i = 0; i < width; i++) {
                vals[2*i] = U.coeffRef(i,height-2);
                vals[2*i+1] = V.coeffRef(i,height-2);
            }
            break;
        case neighbour_position::BOTTOM:
            for (auto i = 0; i < width; i++) {
                vals[2*i] = U.coeffRef(i,1);
                vals[2*i+1] = V.coeffRef(i,1);
            }
            break;
    }

}


void uv_comm (Config &config, MatrixXd  &U, MatrixXd &V, int width, int height) {
    double* verticalVecSend = new double[2*height];
    double* verticalVecRec = new double[2*height];
    double* horVecSend = new double[2*width];
    double* horVecRec = new double[2*width];

    // communicate in i-direction (x-axis)
    if (config.omg_i % 2) {
        // omg_i = 1, 3, 5, 7, ..., (2k + 1 <= iproc)
        // receive - send
        if (config.omg_i != 1) {
            // omg_i = 3, 5, 7, ..., (2k + 1 <= iproc)
            get_vel_from_matrix(U, V, verticalVecSend, width, height, neighbour_position::LEFT);
            MPI_Recv(verticalVecRec, 2*height, MPI_DOUBLE, config.rank_l, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(verticalVecSend, 2*height, MPI_DOUBLE, config.rank_l, 0,MPI_COMM_WORLD);
            store_vel_in_matrix(U, V, verticalVecRec, width, height, neighbour_position::LEFT);
        }
        if (config.omg_i != config.iproc) {
            // omg_i = 1, 3, 5, 7, ..., (2k + 1 < iproc)
            get_vel_from_matrix(U,V,verticalVecSend, width, height, neighbour_position::RIGHT);
            MPI_Recv(verticalVecRec, 2*height, MPI_DOUBLE, config.rank_r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(verticalVecSend ,2*height, MPI_DOUBLE, config.rank_r, 0,MPI_COMM_WORLD);
            store_vel_in_matrix(U,V,verticalVecRec, width, height, neighbour_position::RIGHT);
        }
    } else {
        // omg_i = 2, 4, 6, ..., (2k <= iproc)
        // receive - send
        if (config.omg_i != config.iproc) {
            // omg_i = 2, 4, 6, ..., (2k < iproc)
            get_vel_from_matrix(U,V,verticalVecSend, width, height, neighbour_position::RIGHT);
            MPI_Send(verticalVecSend, 2*height, MPI_DOUBLE, config.rank_r, 0, MPI_COMM_WORLD);
            MPI_Recv(verticalVecRec, 2*height, MPI_DOUBLE, config.rank_r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            store_vel_in_matrix(U,V,verticalVecRec, width, height, neighbour_position::RIGHT);
        }

        // omg_i = 2, 4, 6, ..., (2k <= iproc)
        get_vel_from_matrix(U,V, verticalVecSend, width, height, neighbour_position::LEFT);
        MPI_Send(verticalVecSend, 2*height, MPI_DOUBLE, config.rank_l, 0, MPI_COMM_WORLD);
        MPI_Recv(verticalVecRec, 2*height, MPI_DOUBLE, config.rank_l, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        store_vel_in_matrix(U,V,verticalVecRec, width, height, neighbour_position::LEFT);
    }
    
    // communicate in j-direction (y-axis)
    if (config.omg_j % 2) {
        // omg_j = 1, 3, 5, 7, ..., (2k + 1 <= jproc)
        if (config.omg_j != 1) {
            // omg_j = 3, 5, 7, ..., (2k + 1 <= jproc)
            get_vel_from_matrix(U,V,horVecSend, width, height, neighbour_position::BOTTOM);
            MPI_Recv(horVecRec, 2*width, MPI_DOUBLE, config.rank_b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(horVecSend ,2*width, MPI_DOUBLE, config.rank_b, 0,MPI_COMM_WORLD);
            store_vel_in_matrix(U,V,horVecRec, width, height, neighbour_position::BOTTOM);
        }
        if (config.omg_j != config.jproc) {
            // omg_j = 1, 3, 5, 7, ..., (2k + 1 < jproc)
            get_vel_from_matrix(U,V,horVecSend, width, height, neighbour_position::TOP);
            MPI_Recv(horVecRec, 2*width, MPI_DOUBLE, config.rank_t, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Send(horVecSend, 2*width, MPI_DOUBLE, config.rank_t, 0, MPI_COMM_WORLD);
            store_vel_in_matrix(U,V,horVecRec, width, height, neighbour_position::TOP);
        }
    } else {
        // omg_j = 2, 4, 6, ..., (2k <= jproc)
        if (config.omg_j != config.jproc) {
            // omg_j = 2, 4, 6, ..., (2k < jproc)
            get_vel_from_matrix(U,V,horVecSend, width, height, neighbour_position::TOP);
            MPI_Send(horVecSend, 2*width, MPI_DOUBLE, config.rank_t, 0,MPI_COMM_WORLD);
            MPI_Recv(horVecRec, 2*width, MPI_DOUBLE, config.rank_t, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            store_vel_in_matrix(U,V,horVecRec, width, height, neighbour_position::TOP);
        }

        // omg_j = 2, 4, 6, ..., (2k <= jproc)
        get_vel_from_matrix(U,V,horVecSend, width, height, neighbour_position::BOTTOM);
        MPI_Send(horVecSend, 2*width, MPI_DOUBLE, config.rank_b, 0,MPI_COMM_WORLD);
        MPI_Recv(horVecRec, 2*width, MPI_DOUBLE, config.rank_b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        store_vel_in_matrix(U,V,horVecRec, width, height, neighbour_position::BOTTOM);
    }

    delete[] verticalVecSend;
    delete[] verticalVecRec;
    delete[] horVecSend;
    delete[] horVecRec;
}


void FG_comm (Config &config, MatrixXd  &F, MatrixXd &G, int width, int height) {
    
    double* vertVec = new double[height];
    double* horVec = new double[width];

    // communicate in i-direction (x-axis)
    if (config.omg_i % 2) {
        // omg_i = 1, 3, 5, 7, ..., (2k + 1 <= iproc)
        // receive - send
        if (config.omg_i != 1) {
            // omg_i = 3, 5, 7, ..., (2k + 1 <= iproc)
            MPI_Recv(vertVec, height, MPI_DOUBLE, config.rank_l, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            store_in_matrix(F,vertVec, width, height, neighbour_position::LEFT);
        }
        if (config.omg_i != config.iproc) {
            // omg_i = 1, 3, 5, 7, ..., (2k + 1 < iproc)
            get_from_matrix(F,vertVec, width, height, neighbour_position::RIGHT);
            MPI_Send(vertVec, height, MPI_DOUBLE, config.rank_r, 0, MPI_COMM_WORLD);
        }
    } else {
        // omg_i = 2, 4, 6, ..., (2k <= iproc)
        // receive - send
        if (config.omg_i != config.iproc) {
            // omg_i = 2, 4, 6, ..., (2k < iproc)
            get_from_matrix(F,vertVec, width, height, neighbour_position::RIGHT);
            MPI_Send(vertVec, height, MPI_DOUBLE, config.rank_r, 0, MPI_COMM_WORLD);
        }
        
        // omg_i = 2, 4, 6, ..., (2k <= iproc)
        MPI_Recv(vertVec, height, MPI_DOUBLE, config.rank_l, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        store_in_matrix(F,vertVec, width, height, neighbour_position::LEFT);
    }

    // communicate in j-direction (y-axis)
    if (config.omg_j % 2) {
        // omg_j = 1, 3, 5, 7, ..., (2k + 1 <= jproc)
        // receive - send
        if (config.omg_j != 1) {
            // omg_j = 3, 5, 7, ..., (2k + 1 <= jproc)
            MPI_Recv(horVec, width, MPI_DOUBLE, config.rank_b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            store_in_matrix(G,horVec, width, height, neighbour_position::BOTTOM);
        }
        if (config.omg_j != config.jproc) {
            // omg_j = 1, 3, 5, 7, ..., (2k + 1 < jproc)
            get_from_matrix(G,horVec, width, height, neighbour_position::TOP);
            MPI_Send(horVec, width, MPI_DOUBLE, config.rank_t, 0, MPI_COMM_WORLD);
        }
    } else {
        // omg_j = 2, 4, 6, ..., (2k <= jproc)
        if (config.omg_j != config.jproc) {
            // omg_j = 2, 4, 6, ..., (2k < jproc)
            get_from_matrix(G,horVec, width, height, neighbour_position::TOP);
            MPI_Send(horVec, width, MPI_DOUBLE, config.rank_t, 0, MPI_COMM_WORLD);
        }

        // omg_j = 2, 4, 6, ..., (2k <= jproc)
        MPI_Recv(horVec, width, MPI_DOUBLE, config.rank_b, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        store_in_matrix(G,horVec, width, height, neighbour_position::BOTTOM);
    }

    delete[] vertVec;
    delete[] horVec;
}

