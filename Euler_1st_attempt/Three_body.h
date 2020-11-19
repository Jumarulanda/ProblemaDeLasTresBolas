#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
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

    // Universal gravitational constant
    double G;
    
    // System parameters

    vector<double> Masses;  // {M,m1,m2,m3}
    vector<double> Radii;  
    vector<double> CM;
    
    // Relative coordinate system vectors
    
    vector <double> s1; // s1 = r2 - r3
    vector <double> s2; // s2 = r3 - r1
    vector <double> s3; // s3 = r1 - r2
    
    // And it's magnitudes:
    
    double S1; double S2 ; double S3 ;
    
    // Forces upon relative coordinate system
    
    vector<double> F1;
    vector<double> F2;
    vector<double> F3;
    
    vector<double> g (vector<double>,vector<double>,vector<double>);
    
    // Time tools
    
    double time;
    void update(double);
    void euler_integrator(double);
    bool check_radii();
    
    //Writing in file: {x1,y1,x2,y2,x3,y3,t} 

    void write_on_file(ofstream (*));
    
    
};
    