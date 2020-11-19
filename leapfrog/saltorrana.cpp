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

vector<double> leap :: leap_step(double (*F)(double), double x_n, double v_nHalf, double h, bool is_start = false){
	/* Leapfrog integration step.
	 * 
	 * Inputs:
	 * 		- F(x): differential equation function for the velocity
	 * 		- x_n: position initial condition of present step 
	 * 		- v_nHalf: velocity half-step. If is_start = true, this is the initial condition 
	 * 		           velocity v_0
	 * 		- is_start: true if this is the first ever step computed, if this is the case
	 * 		            then v_nHalf is calculated from v_0, which replaces the input v_nHalf.
	 * 		            Otherwise, this entry is false
	 * 		            
	 *
	 * Outputs:
	 * 		- A double vector containing the next full step for x and v, and the next half
	 * 		  step for v, like this
	 *
	 * 		  {x_n+1, v_n+1, v_n+1/2}
	 *
	 * */

	if (is_start){
		v_nHalf = v_nHalf + h*0.5*F(x_n);
	}	

	double x_np1 = x_n + h*v_nHalf;
	double v_npHalf = v_nHalf + h*F(x_np1);
	double v_np1 = v_nHalf + h*0.5*F(x_np1);

	vector<double> n_state = {x_np1, v_np1, v_npHalf};

	return n_state;
}

void leap :: leap_rangeInt(double (*F)(double), double x_0, double v_0, double h, double t_f){
	/* Leapfrog integration over a range of time.
	 *
	 * Inputs:
	 * 		- F(x): differential equation function for the velocity
	 * 		- x_0: position initial condition
	 * 		- v_0: velocity initial condition
	 * 		- h: integration time-step
	 * 		- t_f: final integration range
	 *
	 * note: t_f should be an integral multiple of h, otherwise, final integration time wont
	 * be reached.
	 *
	 * */

	/* Initial condition vector. v_0 is repeated because v_1/2 has not been computed yet*/
	vector<double> sols = {x_0, v_0, v_0};
	print_vector(sols, 0.);

	/* Compute first leapfrog step. Last argument is true because v_1/2 must be computed*/ 
	sols = leap :: leap_step(*F, sols[0], sols[2], h, true);
	print_vector(sols, h);

	/* Compute the steps that follow until the final computation time is reached */
	for (double t = 2.*h; t <= t_f; t += h){
		sols = leap :: leap_step(F, sols[0], sols[2], h, false);
		print_vector(sols, t);
	}	
}

void print_vector(const vector<double> pvec, double t){
	/* Function to print each integration-step state in terminal. This funciton is
	 * only used by functions from the namespace and cannot be used outside this
	 * script
	 *
	 * */
	
	std::cout << t << ", ";

	int v_size = pvec.size() - 1;

	for (int i = 0; i < v_size - 1; i++) {
		std::cout << pvec[i] << ",";
	}

	std::cout << pvec[v_size-1] << std::endl;
}


/* Verlet algorithm for n degrees of freedom */

vector<double> leap :: leap_step(double (*F)(double<vector>), vector<double> x_n, double<vector> v_nHalf, double h, bool is_start = false){
	/* Leapfrog integration step.
	 * 
	 * Inputs:
	 * 		- F(x): differential equation function for the velocity
	 * 		- x_n: position initial condition of present step 
	 * 		- v_nHalf: velocity half-step. If is_start = true, this is the initial condition 
	 * 		           velocity v_0
	 * 		- is_start: true if this is the first ever step computed, if this is the case
	 * 		            then v_nHalf is calculated from v_0, which replaces the input v_nHalf.
	 * 		            Otherwise, this entry is false
	 * 		            
	 *
	 * Outputs:
	 * 		- A double vector containing the next full step for x and v, and the next half
	 * 		  step for v, like this
	 *
	 * 		  {x_n+1, v_n+1, v_n+1/2}
	 *
	 * */

	if (is_start){
		v_nHalf = v_nHalf + h*0.5*F(x_n);
	}	

	double x_np1 = x_in + h*v_nHalf;
	double v_npHalf = v_inHalf + h*F(x_np1);
	double v_np1 = v_inHalf + h*0.5*F(x_np1);

	vector<double> n_state = {x_np1, v_np1, v_npHalf};

	return n_state;
}

vector<double> update_vector(double (*F)(double<vector>), vector<double> x_n, double<vector> v_nHalf, double h)
