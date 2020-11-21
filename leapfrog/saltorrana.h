#include <iostream>
#include <cmath>
#include <vector>

using std::vector;

namespace leap{

	/* 1-D verlet velocity */
	vector<double> leap_step(double (*)(double), double, double, double);

	/* n-D verlet velocity */
	vector<double> upd_x_n(vector<double>, vector<double> , double);
	vector<double> upd_v_nHalf(vector<double> (*)(vector<double>), vector<double>, vector<double>, double);
	vector<double> upd_v_n(vector<double> (*)(vector<double>), vector<double>, vector<double>, double);

	vector<vector<double>> leap_step_n(vector<double> (*) (vector<double>), vector<double>, vector<double>, double);

}

