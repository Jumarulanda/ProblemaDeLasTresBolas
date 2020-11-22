#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>

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

private:

    // Instance parameters
    
    double G;      // Universal Gravitational constant
    double time;   // Current time of the system
    
    // System parameters

    vector<double> M;   // {M_t,M1,M2,M3}
    vector<double> m;   // Relative masses
    vector<double> R;   // {R1,R2,R3}
    vector<double> CM;  // {x_CM , y_CM}
    
    // Generalized coordinate system
    
    vector <double> q; // {qx1,qy1,qx2,qy2,qx3,qy3} : Relative positions
    vector <double> p; // {px1,py1,px2,py2,px3,py3} : Relative momenta
    
    // State functions
    
    vector<double> r( vector <double>); // Distances between the bodies
    vector<double> g( vector <double>); // g forace vector
    vector<double> f( vector <double>); // Forces upon q   
    
    // Tools
    
    vector <vector<double>> euler_integrator( double, vector<vector<double>>);                    // Integrator: NEEDS REVISION
    
    bool check_radii(vector <double>);  // Checking if no collisions have happened
    void write_on_file(ofstream (*));   // Writing in file: {x1,y1,x2,y2,x3,y3,t} 
    
};
    