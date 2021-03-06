#include "Three_body.h"

//~~~~~~~~~~ Construction of the class ~~~~~~~~~~//

Three_body::Three_body(double start_time , double UGC)
{
    time = start_time;
    G = UGC;
}

void Three_body::set_Masses(vector<double> Masses)
{
    // 4-vector object M

    M = {Masses[0] + Masses[1] + Masses[2]};
    M.insert(M.end(),Masses.begin(),Masses.end());
    
    //3+1-vector objet m
    
    m = {1,0,0,0};
    
    m[1] = M[2]*M[3]/M[0];  
    m[2] = M[1]*M[3]/M[0];
    m[3] = M[1]*M[2]/M[0];
    
}

void Three_body::set_Radii(vector<double> Radii)
{
    // 3+1-vector object
    
    R = {0};  // Dummy radius for accurate correspondence  
    R.insert(R.end(),Radii.begin(),Radii.end());
}

void Three_body::set_Initial_R(vector<double> IR)
{
    /*
    The input is a 6-vector object with the initial positions and velocities of the system.
    Instance must be as follows: {x1,y1,x2,y2,x3,y3}
    */
    
    vector<double> V1(&IR[0],&IR[2]);   // {x1,y1}
    vector<double> V2(&IR[2],&IR[4]);   // {x2,y2}
    vector<double> V3(&IR[4],&IR[6]);   // {x3,y3}
    
    // Determing center of mass coordinates
    
    double x_CM = (M[1]*V1[0] + M[2]*V2[0] + M[3]*V3[0])/M[0];
    double y_CM = (M[1]*V1[1] + M[2]*V2[1] + M[3]*V3[1])/M[0];
    
    R_CM = {x_CM , y_CM};
    
    // Filling the relative coordinate system vectors
    
    vector <double> q1 {}; // q1 = r2 - r3
    vector <double> q2 {}; // q2 = r3 - r1
    vector <double> q3 {}; // q3 = r1 - r2
    
    for (int i = 0; i < 2 ; i++)
    {
        q1.push_back(V2[i] - V3[i]);
        q2.push_back(V3[i] - V1[i]);
        q3.push_back(V1[i] - V2[i]);
    }
        
    q = {q1[0],q1[1],q2[0],q2[1],q3[0],q3[1]}; 
}

void Three_body::set_Initial_K(vector<double> IK)
{
    /*
    The input is a 6-vector object with the initial positions and velocities of the system.
    Instance must be as follows: {kx1,ky1,kx2,ky2,kx3,ky3}
    */
    
    vector<double> P1(&IK[0],&IK[2]);   // {kx1,ky1}
    vector<double> P2(&IK[2],&IK[4]);   // {kx2,ky2}
    vector<double> P3(&IK[4],&IK[6]);   // {kx3,ky3}
    
    // Determing center of mass momentum
    
    double kx_CM = (P1[0] + P2[0] + P3[0])/3;
    double ky_CM = (P1[1] + P2[1] + P3[1])/3;
    
    K_CM = {kx_CM , ky_CM};
    
    // Filling the relative coordinate system vectors
    
    vector <double> p1 {}; // p1 = m1/M2 k2 - m1/M3 k3
    vector <double> p2 {}; // p2 = m2/M3 k3 - m2/M1 k1
    vector <double> p3 {}; // p3 = m3/M1 k1 - m3/M2 k2
    
    for (int i = 0; i < 2 ; i++)
    {
        p1.push_back(m[1] * (P2[i]/M[2] - P3[i]/M[3]));
        p2.push_back(m[2] * (P3[i]/M[3] - P1[i]/M[1]));
        p3.push_back(m[3] * (P1[i]/M[1] - P2[i]/M[2]));
    }
    
    p = {p1[0],p1[1],p2[0],p2[1],p3[0],p3[1]};
    
}


//~~~~~~~~~~ State functions ~~~~~~~~~~//


vector<double> Three_body::r(vector<double> Q)
{
    
    double r23 = sqrt( pow(Q[0],2) + pow(Q[1],2)); // Distance between bodies 2 and 3
    double r13 = sqrt( pow(Q[2],2) + pow(Q[3],2)); // Distance between bodies 1 and 3
    double r12 = sqrt( pow(Q[4],2) + pow(Q[5],2)); // Distance between bodies 1 and 2
        
    return {r23,r13,r12};
}


vector<double> Three_body::g(vector<double> Q)
{
    
    double gx = Q[0]/pow(r(Q)[0],3) + Q[2]/pow(r(Q)[1],3) + Q[4]/pow(r(Q)[2],3);
    double gy = Q[1]/pow(r(Q)[0],3) + Q[3]/pow(r(Q)[1],3) + Q[5]/pow(r(Q)[2],3);
    
    return {gx,gy};
}


vector<double> Three_body::f(vector<double> Q) 
{
 
    double f0 = m[1]*G*( -M[0]*Q[0]/pow(r(Q)[0],3) + M[1]*g(Q)[0] ); // d/dt ( p0 )
    double f1 = m[1]*G*( -M[0]*Q[1]/pow(r(Q)[0],3) + M[1]*g(Q)[1] ); // d/dt ( p1 )
    
    double f2 = m[2]*G*( -M[0]*Q[2]/pow(r(Q)[1],3) + M[2]*g(Q)[0] ); // d/dt ( p2 )
    double f3 = m[2]*G*( -M[0]*Q[3]/pow(r(Q)[1],3) + M[2]*g(Q)[1] ); // d/dt ( p3 )
    
    double f4 = m[3]*G*( -M[0]*Q[4]/pow(r(Q)[2],3) + M[3]*g(Q)[0] ); // d/dt ( p4 )
    double f5 = m[3]*G*( -M[0]*Q[5]/pow(r(Q)[2],3) + M[3]*g(Q)[1] ); // d/dt ( p5 )
    
    return {f0,f1,f2,f3,f4,f5};
    
}

// Integrators

                                                
void Three_body::evol_system(double number_of_steps, double time_step, ofstream *file, char opt)
{
    // Writing file header
    if (file -> is_open())
    {
	  *file << M[1] << "," << M[2] << "," << M[3] << "," << R[1] << "," << R[2] << "," << R[3] << "," << R_CM[0] << "," << R_CM[1] << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "," << 0 << "\n";
    }
    

    for(double i = 0 ; i <= number_of_steps ; i ++)
    {
        // Writes on file the current state of system
        write_on_file(file); 
        
	    // Set state
        vector<vector<double>> State {q,p};
        vector <vector<double>> Updated_State; 
   	
	    // Integrator election 
	    switch (opt)
	    {	
        case 'v': 
		  Updated_State = vel_verlet(time_step, State);
		  break;
               
        case 'e':
		  Updated_State = euler_integrator(time_step, State);
		  break;
               
        case 'r':
		   Updated_State = RK4(time_step, State);
           break;
               
        default:
		  Updated_State = euler_integrator(time_step, State);
        } 
	    
        // Updating generalized coordinate member variables
        
        q = Updated_State[0];
        p = Updated_State[1];
        time += time_step;
        
        // Checking if no collisions have happened
        if ( check_radii(q) ){ break; }
    }
}


vector<vector<double>> Three_body::euler_integrator(double time_step, vector<vector<double>> State)
{
    double h = time_step;
    
    vector<double> Q = State[0];
    vector<double> P = State[1];
    
    vector<double> updated_Q {};
    vector<double> updated_P {};
    
    for (int i = 0; i < Q.size(); i++)
    {
        updated_Q.push_back( Q[i] + h*(P[i]/m[int(i/2. + 1)]));
        updated_P.push_back( P[i] + h*f(q)[i] );  
    }
    
    return {updated_Q , updated_P};
}


vector<vector<double>> Three_body::vel_verlet(double time_step, vector<vector<double>> State) {
    double h = time_step;
    
    vector<double> Q = State[0];
    vector<double> P = State[1];
    
    vector<double> updated_Q {};
    vector<double> updated_P {};
    
    vector<double> P_nHalf;

    for (int i = 0; i < Q.size(); i++) {
      P_nHalf.push_back(P[i] + h*f(q)[i]*0.5);
      updated_Q.push_back(Q[i] + h*P_nHalf[i]/m[int(i/2. + 1)]);
    }

    for (int i=0; i < Q.size(); i++) {
      updated_P.push_back(P_nHalf[i] + h*0.5*f(Q)[i]);
    }
    
    return {updated_Q , updated_P};
}
vector<vector<double>> Three_body:: RK4 (double h, vector<vector<double>> State)
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
  vector <double> vector_state_q = State[0];
  vector <double> vector_state_p = State[1];
  for (int i = 0; i<4;i++)
    {
      vector_state_q =  v_q(h,c[i], vector_state_q,vector_state_p);
      vector_state_p = v_p(h,d[i], vector_state_p);	      
    }


  return { vector_state_q,  vector_state_p};
}


//update states of q for one order 
vector <double> Three_body :: v_q(double h, double c,vector <double> vector_state_q,vector <double> vector_state_p)
{
  int l = vector_state_q.size();
  vector <double> q1_vec;
  for (int i=0; i<l;i++)
    {
    double q1 = vector_state_q[i] + c*h*(vector_state_p[i]/m[int(i/2. + 1)]);
    q1_vec.push_back(q1);
    }
  return q1_vec;
}

//update states of p for one order
vector <double> Three_body :: v_p(double h, double d,vector <double> vector_state_p)
{
  int l = vector_state_p.size();
  //cout << l <<endl;
  vector <double> p_vec;
  for (int i=0; i<l;i++)
    {
      double p1 = vector_state_p[i] + d*h*f(q)[i];
      p_vec.push_back(p1);
    }
  return  p_vec; 
}

  


                                                    
//~~~~~~~~~~~ Tools ~~~~~~~~~~//
    
void Three_body::write_on_file(ofstream *file)
{
    /*
    Recovering the absolute coordinate system:
    The matrix which determines the change of the system is not inversible,
    therefore one must stablish a second condition.
    In this case, is proposed to fix the center of mass of the system:
    
                        M1*x1 + M2*x2 + M3*x3 = M*x_CM
                        M1*y1 + M2*y2 + M3*y3 = M*y_CM
                        
                        M1*kx1 + M2*kx2 + M3*kx3 = 0
                        M1*ky1 + M2*ky2 + M3*ky3 = 0
                        
    And also:
    
                          q1 =  0*r1 + 1*r2 - 1*r3
                          q2 = -1*r1 + 0*r2 + 1*r3
                          q3 =  1*r1 - 1*r2 + 0*r3
    */
        
    double x1 = (-M[3]*q[2] + M[2]*q[4])/M[0] + R_CM[0] ;//+ K_CM[0]*time/M[1];
    double y1 = (-M[3]*q[3] + M[2]*q[5])/M[0] + R_CM[1] ;//+ K_CM[1]*time/M[1];
    
    
    double kx1 = - p[2] + p[4] ;// + K_CM[0];
    double ky1 = - p[3] + p[5] ;// + K_CM[1];
            
    double x2 = (-M[3]*q[2] - (M[1] + M[3])*q[4])/M[0] + R_CM[0] ;//+ K_CM[0]*time/M[2];
    double y2 = (-M[3]*q[3] - (M[1] + M[3])*q[5])/M[0] + R_CM[1] ;//+ K_CM[1]*time/M[2];
    
    
    double kx2 = -m[1]/m[2] * p[2] - (1 + m[1]/m[3]) * p[4] ;//+ K_CM[0];
    double ky2 = -m[1]/m[2] * p[3] - (1 + m[1]/m[3]) * p[5] ;//+ K_CM[1];
    
    double x3 = ((M[1] + M[2])*q[2] + M[2]*q[4])/M[0] + R_CM[0] ;// + K_CM[0]*time/M[3];
    double y3 = ((M[1] + M[2])*q[3] + M[2]*q[5])/M[0] + R_CM[1] ;// + K_CM[1]*time/M[3];
    
    double kx3 = (1 + m[1]/m[2])*p[2] + m[1]/m[3] * p[4] ;//+ K_CM[0];
    double ky3 = (1 + m[1]/m[2])*p[3] + m[1]/m[3] * p[5] ;//+ K_CM[1];
    
    
    if (file -> is_open())
    {
	  *file << x1 << "," << y1 << "," << x2 << "," << y2 << "," << x3 << "," << y3 << "," << kx1 << "," << ky1 << "," << kx2 << "," << ky2 << "," << kx3 << "," << ky3 << "," << time << "\n";
    }

}
    

bool Three_body::check_radii(vector <double> Q)
{
    bool sp12 = R[1] + R[2] > r(Q)[2];  // R1 + R2 > R12
    bool sp13 = R[1] + R[3] > r(Q)[1];  // R1 + R3 > R23
    bool sp23 = R[2] + R[3] > r(Q)[0];  // R2 + R3 > R23
    
    if (sp12 || sp13 || sp23)
    {
        cout << "There was a collision" << endl;
    }
    return sp12 || sp13 || sp23;
}







