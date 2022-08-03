#include "init.hpp"
#include "datastructures.hpp"
#include "solver.hpp"
#include <fstream>
#ifdef gpp9
// gcc Version >= 9
#include "filesystem"
namespace filesystem = std::filesystem;
#else
// gcc Version < 9
#include "experimental/filesystem"
namespace filesystem = std::experimental::filesystem;
#endif // DEBUG

int read_parameters(filesystem::path filepath, Config& conf) {
	std::ifstream file(filepath);
	if (!file.is_open()) return -1;
	std::string var;
	std::string solver;
	while (!file.eof() && file.good()) {
		file >> var;
		if (var[0] == '#') {     /* ignore comment line*/
			file.ignore(MAX_LINE_LENGTH, '\n');
		}
		else {
			if (var == "xlength")  file >> conf.xlength;
			if (var == "ylength")  file >> conf.ylength;
			if (var == "Re")       file >> conf.Re;
			if (var == "t_end")    file >> conf.t_end;
			if (var == "dt")       file >> conf.dt;
			if (var == "omg")      file >> conf.omg;
			if (var == "eps")      file >> conf.eps;
			if (var == "tau")      file >> conf.tau;
			if (var == "alpha")    file >> conf.alpha;
			if (var == "dt_value") file >> conf.dt_value;
			if (var == "UI")       file >> conf.UI;
			if (var == "VI")       file >> conf.VI;
			if (var == "GX")       file >> conf.GX;
			if (var == "GY")       file >> conf.GY;
			if (var == "PI")       file >> conf.PI;
			if (var == "itermax")  file >> conf.itermax;
			if (var == "imax")     file >> conf.imax;
			if (var == "jmax")     file >> conf.jmax;
			if (var == "problem")  file >> conf.problem;
			if (var == "geometry") file >> conf.geometry;
            if (var == "PR")       file >> conf.PR;
            if (var == "TI")       file >> conf.TI;
			if (var == "T_h")      file >> conf.T_h;
			if (var == "T_c")      file >> conf.T_c;
            if (var == "beta")     file >> conf.beta;
            if (var == "calcTemp" ) file >> conf.calcTemp;
            if (var == "iproc") 	file >> conf.iproc;
            if (var == "jproc") 	file >> conf.jproc;
			if (var == "levels") 	file >> conf.levels;
			if (var == "solver") 	file >> solver;
		}
	}

	if (!file.good() && !file.eof()) return -1;

	conf.solver = Solver::string2type(solver);
	conf.dx = conf.xlength / (double)(conf.imax);
	conf.dy = conf.ylength / (double)(conf.jmax);

	return 1;
}
