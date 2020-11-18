#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

class Three_body
{   
public:
        
    Three_body();
    void set_Masses(vector<double>);
    void set_Radii(vector<double>);
    void set_Initial_Conditions(vector<double>);

// private:

    // Universal gravitational constant
    double G = 1;

    vector<double> Masses;
    double TM; // Total mass of the system
    vector<double> Radii;  
    
    // Relative coordinate system vectors
    
    vector <double> s1; 
    vector <double> s2;
    vector <double> s3;
    
    // And it's magnitudes
    
    double S1; double S2 ; double S3 ;
    
    // System functions
    
    vector<double> F1 ();
    vector<double> F2 ();
    vector<double> F3 ();

    vector<double> g (vector<double>,vector<double>,vector<double>);
    
    // Time tools
    
    double time = 0;
    void update_distances();

    
};
    