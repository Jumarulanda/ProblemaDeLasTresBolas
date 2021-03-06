#include "Three_body.h"

int main()
{   
    /*
    The Lagrange solutions correspond to elliptic periodic orbits that maintains a equilateral triangle disposition
    */
    
    double r = 30;
    double theta = 0;
    double const pi = 3.1415;
    
    double x1 = r*cos(theta);
    double y1 = r*sin(theta);
    
    double x2 = r*cos(theta + 2*pi/3.);
    double y2 = r*sin(theta + 2*pi/3.);
    
    double x3 = r*cos(theta + 4*pi/3.);
    double y3 = r*sin(theta + 4*pi/3.);

    // System parameters
    
    vector <double> M = {10,1,1};  // Masses
    vector <double> R = {10,1,1};  // Radii
    
    // Initial conditions as {x,y,vx,vy}
    
    vector <double> IC1 = {x1,y1,0,0};
    vector <double> IC2 = {x2,y2,-0.2,-0.2};
    vector <double> IC3 = {x3,y3,0.2,-0.2};
    

    vector <double> IC(IC1);
    IC.insert(IC.end(),IC2.begin(),IC2.end());
    IC.insert(IC.end(),IC3.begin(),IC3.end());
    
    // Instance the Three_body object
    
    Three_body TB(0,1); // Initial_time , Value of G
    
    TB.set_Masses(M);
    TB.set_Radii(R);
    TB.set_Initial_Conditions(IC);
    
    ofstream file ("prueba1.txt"); // File for storing solutions
    
    TB.evol_system(10000,0.1,&file); // Numer of steps, step size
    
    return 0;
}