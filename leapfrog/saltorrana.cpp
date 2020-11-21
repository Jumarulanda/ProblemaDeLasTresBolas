#include "saltorrana.h"

/* Leapfrog integration method. It is a generalization of the Euler method, and 
 * of a higher order. This method is also simplectic in the sense that it 
 * conserves the area of the parallelogram (q,p) of the phase space when an 
 * integration step is performed */

/* This method is applicable to second order differential equations of the form 
 *
 * x'' = F(x)
 *
 * where there is no first order term in x. The system of two coupled first order
 * differential equations may be constructed
 *
 * x' = v
 * v' = F(x)
 *
 * Notice that F(x) does not depend on the velocities, but only on the positions.
 *
 * This integration method is designed on the fact that the slope between the 
 * points (x0,v0) and (x1,v1) is a better approximation to the slope at (x 1/2, v 1/2) 
 * than to the extremals themselves. Therefore, this method performs an Euler 
 * integration of x_n using the next half-step of v, that is
 *
 * x_n+1 = x_n + h * v_n+1/2
 *
 * The next v half step is calculated as
 *
 * v_n+3/2 = v_n+1/2 + h * F(x_n+1)
 *
 * Where, as t_n+1 is at a distance 1/2 from t_n+3/2, v_n+3/2 is also computed by 
 * a frog's leap.
 *
 * However, the initial conditions are at (x0,v0), not (x0,v_1/2). The velovity at 
 * v_1/2 is computed by a half Euler step. Even though Euler is of order 1, as only this
 * computation is made, the method as a whole is still of order 2. Then, v_1/2 is computed as
 *
 * v_1/2 = v_0 + h/2 * F(x_0)
 *
 * Likewise, the full step on v is calculated by a Euler half-step at v_n+1/2
 *
 * v_n+1 = v_n+1/2 + h/2 * F(x_n+1)
 *
 * */

/* ------------------------------------------------------------------------------------------------ */

void print_vector(const vector<double>, double);

/* Integration method */

/* Verlet algorithm for 1 degree of freedom */

vector<double> leap :: leap_step(double (*F)(double), double x_n, double v_n, double h){
	/* Leapfrog integration step.
	 * 
	 * Inputs:
	 * 		- F(x): differential equation function for the velocity
	 * 		- x_n: position initial condition of present step 
	 * 		- v_n: velocity at present step.
	 *
	 * Outputs:
	 * 		- A double vector containing the next full step for x and v, and the next half
	 * 		  step for v, like this
	 *
	 * 		  {x_n+1, v_n+1, v_n+1/2}
	 *
	 * */

	double v_nHalf = v_n + h*F(x_n)*0.5;	
	double x_np1 = x_n + h*v_nHalf;
	double v_np1 = v_nHalf + h*0.5*F(x_np1);

	vector<double> n_state = {x_np1, v_np1};

	return n_state;
}


/* Verlet algorithm for n degrees of freedom */

vector<double> leap :: upd_x_n(vector<double> x_n, vector<double> v_nHalf, double h) {
	/* Function that updates the n generalized positions. */
	
	/* Inputs: */
	/* 	- x_n: vector of generalzied positions at time t_n */
	/* 	- v_nHalf: vector of generalized velocities at time t_n+1/2 */
	/* 	- h: integration time step */ 

	int deg_of_freedom = x_n.size();
	vector<double> x_np1;

	for (int i = 0; i < deg_of_freedom; i++){
		x_np1.push_back(x_n[i] + h*v_nHalf[i]);
	}

	return x_np1;
}

vector<double> leap :: upd_v_nHalf(vector<double> (*F)(vector<double>), vector<double> x_n, vector<double> v_n, double h) {
	/* Function that updates the n generalized half velocities. */
	
	/* Inputs: */
	/* 	- F(x): ODE system */
	/* 	- x_n: vector of generalzied positions at time t_n */
	/* 	- v_nHalf: vector of generalized velocities at time t_n+1/2 
	 * 	- v_n: vetor of generalized velocites at time t_n */
	/* 	- h: integration time step */ 

	int deg_of_freedom = x_n.size();
	vector<double> v_nHalf;

	for (int i = 0; i < deg_of_freedom; i++){
		v_nHalf.push_back(v_n[i] + h*F(x_n)[i]*0.5);
	}

	return v_nHalf;
}

vector<double> leap :: upd_v_n(vector<double> (*F)(vector<double>), vector<double> x_np1, vector<double> v_nHalf, double h) {
	/* Function that updates the n generalized half velocities. */
	
	/* Inputs: */
	/* 	- F(x): ODE system */
	/* 	- x_np1: vector of generalzied positions at time t_n+1 */
	/* 	- v_nHalf: vector of generalized velocities at time t_n+1/2 
	 * 	- v_n: vetor of generalized velocites at time t_n */
	/* 	- h: integration time step */ 

	int deg_of_freedom = x_np1.size();
	vector<double> v_np1;

	for (int i = 0; i < deg_of_freedom; i++){
		v_np1.push_back(v_nHalf[i] + h*F(x_np1)[i]*0.5);
	}

	return v_np1;
}

vector<vector<double>> leap :: leap_step_n(vector<double> (*F) (vector<double>), vector<double> x_n, vector<double> v_n, double h) {
	/* Function to integrate a step on all the degrees of freedom of the system
	 * Inputs:
	 * 		- F(x): ODE system
	 * 		- x_n: vector of generalized positions at t_n
	 * 		- v_nHalf: vector of generalized velocities at t_n+1/2
	 * 		- v_n; vector of generalized velocities at t_n
	 * 		- h: integratoin time step
	 *
	 * */

	vector<double> v_nHalf = upd_v_nHalf(F, x_n, v_n, h);
	vector<double> x_np1 = upd_x_n(x_n, v_nHalf, h);
	vector<double> v_np1 = upd_v_n(F, x_np1, v_nHalf, h);

	vector<vector<double>> sol_tnp1 = {x_np1, v_np1};

	return sol_tnp1;
}	
