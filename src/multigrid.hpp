#ifndef __MULTIGRID_H_
#define __MULTIGRID_H_

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <vector>
#include "enums.hpp"
#include "solver.hpp"

using namespace Eigen;

class Level {
    public:
        Level(int level, int imax, int jmax, double dx, double dy, matrix<cell_type>& types);
		Level(int level, int imax, int jmax, double dx, double dy, matrix<cell_type>& types_, MatrixXd& p, MatrixXd& rhs);

		MatrixXd& x();
		MatrixXd& b();
		MatrixXd& e();
		MatrixXd& res();
		matrix<cell_type>& types();

	private:
		friend class Multigrid;

        int level;
        int imax;
        int jmax;

        double dx;
        double dy;

		matrix<cell_type>* _types;
        MatrixXd* _x = nullptr;
        MatrixXd* _b = nullptr;
        MatrixXd* _e = nullptr;
		MatrixXd* _res = nullptr;
};

class Multigrid: public Solver {
    public:
        Multigrid(Config& config, MatrixXd& p, MatrixXd& rhs, matrix<cell_type>& types);
        void solve(int& it, double& res);
		Level& levels(int l);

	private:
		MatrixXd _residual;
        FullPivLU<MatrixXd>_solver;
        MatrixXd _A;
        
		std::vector<Level*> _levels;

		void v_cycle(int level);
		void smoothing(MatrixXd& x, MatrixXd& b, matrix<cell_type>& types, double dx, double dy, bool boundary);
		void restriction_fullweight(MatrixXd& fine, MatrixXd& coarse);
		void prolongate(MatrixXd& fine, MatrixXd& coarse);
		void gaussSeidel(MatrixXd &P, MatrixXd &RS, matrix<cell_type>& _types, double dx, double dy);

};

#endif