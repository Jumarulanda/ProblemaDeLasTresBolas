#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>

#include "eu_int.h"

using namespace std::placeholders;
using namespace std;

class Three_body
{   
public:
        
    Three_body(double, double);
    
    void set_Masses(vector<double>);
    void set_Radii(vector<double>);
    void set_Initial_Conditions(vector<double>);
    
    void evol_system(double, double, ofstream (*));

//private:

    // Instance parameters
    
    double G;      // Universal Gravitational constant
    double time;   // Current time of the system
    
    // System parameters

    vector<double> Masses;  // {M,m1,m2,m3}
    vector<double> Radii;   // {R1,R2,R3}
    vector<double> CM;      // {x_CM , y_CM}
    
    // Generalized coordinate system
    
    vector <double> q; // {sx1,sy1,sx2,sy2,sx3,sy3} : Relative positions
    vector <double> p; // {ux1,uy1,ux2,uy2,ux3,uy3} : Relative velocities
    
    // State functions
    
    vector <double> r( vector <double>); // Distances between the bodies
    vector <double> g( vector <double>); // g forace vector
        
    vector<double> f( vector <double>); // Forces upon q   
    
    // Tools
    
    vector <vector<double>> euler_integrator( double, vector<vector<double>>);                    // Integrator: NEEDS REVISION
    
    bool check_radii(vector <double>);  // Checking if no collisions have happened
    void write_on_file(ofstream (*));   // Writing in file: {x1,y1,x2,y2,x3,y3,t} 
    
	/* Friend functions */

	friend vector<vector<double>> eu_int :: eu_int_step(vector<double> (*) (vector<double>), double, vector<vector<double>>);
};
    
