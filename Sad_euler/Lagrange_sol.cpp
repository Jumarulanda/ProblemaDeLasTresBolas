#include "Three_body.h"

int main()
{   
    /*
    The Lagrange solutions correspond to elliptic periodic orbits that maintains a equilateral triangle disposition
    */
        
    double PI = 3.14159265359;
    double G = 39.43279791722677;
    
    double r = 30;
    double theta = 0;
    
    double x1 = r*cos(theta);
    double y1 = r*sin(theta);
    
    double x2 = r*cos(theta + 2*PI/3.);
    double y2 = r*sin(theta + 2*PI/3.);
    
    double x3 = r*cos(theta + 4*PI/3.);
    double y3 = r*sin(theta + 4*PI/3.);

    // System parameters
    
    vector <double> M = {1,1,1};  // Masses
    vector <double> R = {0.1,0.1,0.1};  // Radii
    
    // Initial conditions 
    
    vector <double> IR = {-10,0,0,0,10,0};
    vector <double> IK = {0,1,0,-1,0,1}; 
    
    // Instance the Three_body object
    
    Three_body TB(0,G); // Initial_time , Value of G
    
    TB.set_Masses(M);
    TB.set_Radii(R);
    
    TB.set_Initial_R(IR);
    TB.set_Initial_K(IK);
            
    ofstream file ("prueba1.txt");      // File for storing solutions
    TB.evol_system(500,0.1,&file,'r');  // Numer of steps, step size
    
    
    return 0;
}
