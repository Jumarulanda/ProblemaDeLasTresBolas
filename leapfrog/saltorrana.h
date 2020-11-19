#include <iostream>
#include <cmath>
#include <vector>

using std::vector;

namespace leap{

	vector<double> leap_step(double (*)(double), double, double, double, bool);
		
	void leap_rangeInt(double (*)(double), double, double, double, double);
}
