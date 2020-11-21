#include <vector>
#include <iostream>
#include <cmath>
#include "RK4_symplectic.h"

using namespace std;

vector <vector<double>> f(vector <double>,vector <double>);

int main()
{
  vector <double> init_q {0};
  vector <double> init_p {1};
  rk4_s integrate(init_q,init_p,f);
  integrate.rk4(0.1,100);
  return 0;
}

vector <vector<double>> f(vector <double> q,vector <double> p)
{
  double k = 1;
  double m = 1;

  double xpoint = p[0]/m;
  double ppoint = -k*q[0];
  vector <vector<double>> variables {{xpoint},{ppoint}};
  return variables;
  
}
