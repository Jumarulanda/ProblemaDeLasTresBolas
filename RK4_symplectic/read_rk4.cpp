#include <vector>
#include <iostream>
#include <cmath>
#include "RK4_symplectic.h"

using namespace std;

vector <vector<double>> f(double,vector <double>,vector <double>);
vector <vector<double>> g(double,vector <double>,vector <double>);

int main()
{
  vector <double> init_q {0}; //initial values of q
  vector <double> init_p {1}; //initial values of p
  rk4_s integrate(init_q,init_p,f); //Initializing RK4_S class
  integrate.rk4(0.1,100); //solving equation
  return 0;
}

//Hamilton equations 
vector <vector<double>> f(double t,vector <double> q,vector <double> p)
{
  double k = 1;
  double m = 1;

  double xpoint = p[0]/m; 
  double ppoint = -k*q[0];
  vector <vector<double>> variables {{xpoint},{ppoint}};
  return variables;
}




