#include "saltorrana.h"

using std::sin;
using std::pow;
using std::cout;
using std::endl;

double omeg_p(double);
void print_vector(const vector<double>);

double PI = 3.14159265359;
double os_frq = 2; 

int main(){

	double th_0 = 1., omeg_0 = 0.;
	double h = 2./os_frq/25.;

	vector<double> sols = leap :: leap_step(*omeg_p, th_0, omeg_0, h);

	for (int i=0; i < 100; i++) {
		sols = leap :: leap_step(*omeg_p, sols[0], sols[1], h);
		print_vector(sols);
	}

	/* leap :: leap_rangeInt(*omeg_p, th_0, omeg_0, h, 50.*h); */

	return 0;
}

double omeg_p(double theta){
	double thpp = -1*pow(os_frq,2)*theta;

	return thpp;
}

void print_vector(const vector<double> pvec){

	int v_size = pvec.size();

	for (int i = 0; i < v_size - 1; i++) {
		cout << pvec[i] << ",";
	}

	cout << pvec[v_size-1] << endl;
}
