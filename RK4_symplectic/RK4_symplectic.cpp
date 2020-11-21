#include "RK4_symplectic.h"


rk4_s :: rk4_s (vector  <double> init_q,vector  <double> init_p,vector <vector<double>> (*func)(double, vector <double>,vector<double>))
{
  vector_state_q  = init_q; //Initializing q states
  vector_state_p  = init_p;  //Initializing p states
  f = func; 
}

//Solution for the system of equation
vector <vector <double>> rk4_s :: rk4(const double h,const int N)
{
  t = 0; //initializing the time
  vector <vector <double>> results; //vector with solutions  
  for (int i = 0; i<N; i++)
    {
      t+=h;
      vector <double> vec =  next(h,t); //finding the nexts q and p
      print(vec);
      cout << t << endl;;
      results.push_back(vec); //adding q_{i+1} and p_{i+1} to the solution in every time
    }
  return results;
}

//printing the results 
void rk4_s :: print(vector<double> const &input)
{
  for (int i = 0; i < input.size(); i++)
    {
      cout << input.at(i) << ' ';
    }
  
}

//finding  next steps for q and p
vector <double> rk4_s :: next(const double h,double t)
{
  //Definitions of constants that define the symplectic method
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

  //updating states of q and p with the fourth order runge kuta
  for (int i = 0; i<4;i++)
    {
      v_q(h,t,c[i]);
      v_p(h,t,d[i]);	      
    }
  vector <double> result;
  result = vector_state_q;
  result.insert(result.end(),vector_state_p.begin(),vector_state_p.end());
  return result;
}

//update states of q for one order 
void rk4_s :: v_q(double h,double t, double c)
{
  int l = vector_state_q.size();
  vector <double> q1_vec;
  for (int i=0; i<l;i++)
    {
      vector <vector<double>> fun = f(t,vector_state_q,vector_state_p);
      double q1 = vector_state_q[i] + c*h*fun[0][i];
      q1_vec.push_back(q1);
    }
  vector_state_q = q1_vec;
}

//update states of p for one order
void rk4_s :: v_p(double h,double t,double d)
{
  int l = vector_state_p.size();
  vector <double> p_vec;
  for (int i=0; i<l;i++)
    {
      double p1 = vector_state_p[i] + d*h*f(t,vector_state_q,vector_state_p)[1][0];
      p_vec.push_back(p1);
    }
  vector_state_p = p_vec; 
}
