#include "RK4_symplectic.h"


rk4_s :: rk4_s (vector  <double> init_q,vector  <double> init_p,vector <vector<double>> (*func)(vector <double>,vector<double>))
{
  vector_state_q  = init_q;
  vector_state_p  = init_p;
  f = func;
}


vector <vector <double>> rk4_s :: rk4(const double h,const int N)
{
  vector <vector <double>> results; 
  for (int i = 0; i<N; i++)
    {
      vector <double> vec =  next(h);
      print(vec);
      cout<< "\n";
      results.push_back(vec);
    }
  return results;
}
void rk4_s :: print(vector<double> const &input)
{
  for (int i = 0; i < input.size(); i++)
    {
      cout << input.at(i) << ' ';
    }
}

vector <double> rk4_s :: next(const double h)
{
  double c1 = 1./(2.*(2.-pow(2.,1./3.)));
  double c2 = (1.-pow(2.,1./3.))/(2*(2.-pow(2.,1./3.)));
  double c3 = c2;
  double c4 = c1;
  vector <double> c {c1,c2,c3,c4};

  double d1  = 1./(2.-pow(2.,1./3.));
  double d2  = -pow(2.,1./3.)/(2.-pow(2.,1./3.));
  double d3 = d1;
  double d4 = 0;
  vector <double> d {d1,d2,d3,d4};
  
  for (int i = 0; i<4;i++)
    {
      v_q(h,c[i]);
      v_p(h,d[i]);	      
    }
  vector <double> result;
  result = vector_state_q;
  result.insert(result.end(),vector_state_p.begin(),vector_state_p.end());
  return result;
}
void rk4_s :: v_q(double h,double c)
{
  int l = vector_state_q.size();
  vector <double> q1_vec;
  for (int i=0; i<l;i++)
    {
      vector <vector<double>> fun = f(vector_state_q,vector_state_p);
      double q1 = vector_state_q[i] + c*h*fun[0][i];
      q1_vec.push_back(q1);
    }
  vector_state_q = q1_vec;
}

void rk4_s :: v_p(double h,double d)
{
  int l = vector_state_p.size();
  vector <double> p_vec;
  for (int i=0; i<l;i++)
    {
      double p1 = vector_state_p[i] + d*h*f(vector_state_q,vector_state_p)[1][0];
      p_vec.push_back(p1);
    }
  vector_state_p = p_vec; 
}