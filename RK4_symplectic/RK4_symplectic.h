#include <vector>
#include <iostream>
#include <cmath>

using namespace std;

class rk4_s
{
 public:
  rk4_s(vector  <double>,vector  <double>,vector <vector<double>> (*)(double,vector <double>,vector<double>));
  vector <vector <double>> rk4(const double, const int);
 private:
  vector <vector<double>> (*f)(double, vector<double>,vector<double>);
  vector <double> vector_state_q;
  vector <double> vector_state_p;
  vector <double> next(const double,double);
  void v_q(double,double,double);
  void v_p(double,double,double);
  void print(vector <double> const &);
  double t;
};
