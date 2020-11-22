#include "Three_body.h"

//~~~~~~~~~~ Construction of the class ~~~~~~~~~~//

Three_body::Three_body(double start_time , double UGC)
{
    time = start_time;
    G = UGC;
}

void Three_body::set_Masses(vector<double> M)
{
    // 4-vector object

    Masses = {M[0] + M[1] + M[2]};
    Masses.insert(Masses.end(),M.begin(),M.end());
}

void Three_body::set_Radii(vector<double> R)
{
    // 3-vector object
    Radii = R;
}

void Three_body::set_Initial_Conditions(vector<double> IC)
{
    /*
    The input is a 12-vecttor object with the initial positions and velocities of the system.
    Instance must be as follows: {x1,y1,vx1,vy1,x2,y2,vx2,vy2,x3,y3,vx3,vy3}
    */
    
    vector<double> V1(&IC[0],&IC[4]);   // {x1,y1,vx1,vy1}
    vector<double> V2(&IC[4],&IC[8]);   // {x2,y2,vx2,vy2}
    vector<double> V3(&IC[8],&IC[12]);  // {x3,y3,vx3,vy3}
    
    // Determing center of mass coordinates
    
    double x_CM = (Masses[1]*V1[0] + Masses[2]*V2[0] + Masses[3]*V3[0])/Masses[0];
    double y_CM = (Masses[1]*V1[1] + Masses[2]*V2[1] + Masses[3]*V3[1])/Masses[0];
    
    CM = {x_CM , y_CM};
    
    // Filling the relative coordinate system vectors
     
    vector <double> s1 {}; // s1 = r2 - r3
    vector <double> s2 {}; // s2 = r3 - r1
    vector <double> s3 {}; // s3 = r1 - r2
    
    for (int i = 0; i < 4 ; i++)
    {
        s1.push_back(V2[i] - V3[i]);
        s2.push_back(V3[i] - V1[i]);
        s3.push_back(V1[i] - V2[i]);
    }
    
    q = {s1[0],s1[1],s2[0],s2[1],s3[0],s3[1]};
    p = {s1[2],s1[3],s2[2],s2[3],s3[2],s3[3]};

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
 
    double f0 = G*( -Masses[0]*Q[0]/pow(r(Q)[0],3) + Masses[1]*g(Q)[0] ); // Force over sx1
    double f1 = G*( -Masses[0]*Q[1]/pow(r(Q)[0],3) + Masses[1]*g(Q)[1] ); // Force over sy1
    
    double f2 = G*( -Masses[0]*Q[2]/pow(r(Q)[1],3) + Masses[2]*g(Q)[0] ); // Force over sx2
    double f3 = G*( -Masses[0]*Q[3]/pow(r(Q)[1],3) + Masses[2]*g(Q)[1] ); // Force over sy2
    
    double f4 = G*( -Masses[0]*Q[4]/pow(r(Q)[2],3) + Masses[3]*g(Q)[0] ); // Force over sx3
    double f5 = G*( -Masses[0]*Q[5]/pow(r(Q)[2],3) + Masses[3]*g(Q)[1] ); // Force over sy3
    
    return {f0,f1,f2,f3,f4,f5};
    
}

// Integrators

                                                

void Three_body::evol_system(double number_of_steps, double time_step, ofstream *file)
{
    // Writing file header
    if (file -> is_open())
    {
	  *file << Masses[0] << "," << Masses[1] << "," << Masses[2] << "," << Radii[0] << "," << Radii[1] << "," << Radii[2] << "," << 0 << "\n";
    }
    

    for(double i = 0 ; i <= number_of_steps ; i ++)
    {
        // Writes on file the current state of system
        write_on_file(file); 
        
        vector<vector<double>> State {q,p};
                
        vector <vector<double>> Updated_State; 
        
        /*
        El problema está aquí, el euler_integrator debería tomar también como entrada la referencia de la funcion f donde están evaluadas las fuerzas del sistema. Sin embargo, como no fui capaz, entonces euler mira directamene f, sin tomarla como entrada.
        */ 
        
        Updated_State = euler_integrator(time_step, State);
            
        // Updating generalized coordinate member variables
        
        q = Updated_State[0];
        p = Updated_State[1];
        
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
        updated_Q.push_back( Q[i] + h*P[i]);
        updated_P.push_back( P[i] + h*f(q)[i] );  // Aquí f debería ser una funcion por referencia y no esa f es una función miembro de la clase.
    }
    
    return {updated_Q , updated_P};
}
    
                                                    
//~~~~~~~~~~~ Tools ~~~~~~~~~~//
    
void Three_body::write_on_file(ofstream *file)
{
    /*
    Recovering the absolute coordinate system:
    The matrix which determines the change of the system is not inversible,
    therefore one must stablish a second condition.
    In this case, is proposed to fix the center of mass of the system:
    
                        m1x1 + m2x2 + m3x3 = M*x_CM
                        m1y1 + m2y2 + m3y3 = M*y_CM
                        
    And also:
    
                          s1 =  0*r1 + 1*r2 - 1*r3
                          s2 = -1*r1 + 0*r2 + 1*r3
                          s3 =  1*r1 - 1*r2 + 0*r3
    */
        
    double x1 = (-Masses[3]*q[2] + Masses[2]*q[4])/Masses[0] + CM[0];
    double y1 = (-Masses[3]*q[3] + Masses[2]*q[5])/Masses[0] + CM[1];
        
    double x2 = (-Masses[3]*q[2] - (Masses[1] + Masses[3])*q[4])/Masses[0] + CM[0];
    double y2 = (-Masses[3]*q[3] - (Masses[1] + Masses[3])*q[5])/Masses[0] + CM[1];
    
    double x3 = ((Masses[1] + Masses[2])*q[2] + Masses[2]*q[4])/Masses[0] + CM[0];
    double y3 = ((Masses[1] + Masses[2])*q[3] + Masses[2]*q[5])/Masses[0] + CM[1];
    
    if (file -> is_open())
    {
	  *file << x1 << "," << y1 << "," << x2 << "," << y2 << "," << x3 << "," << y3 << "," << time << "\n";
    }

}
    

bool Three_body::check_radii(vector <double> Q)
{
    bool sp12 = Radii[0] + Radii[1] > r(Q)[0];  // R1 + R2 > R12
    bool sp13 = Radii[0] + Radii[2] > r(Q)[1];  // R1 + R3 > R23
    bool sp23 = Radii[1] + Radii[2] > r(Q)[2];  // R2 + R3 > R23
    
    if (sp12 || sp13 || sp23)
    {
        cout << "There was a collition" << endl;
    }
    return sp12 || sp13 || sp23;
}







