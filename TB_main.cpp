#include "Three_body.h"

int main()
{   
    // System parameters
    
    vector <double> M = {1,1,1};
    vector <double> R = {1,1,1};
    
    // Initial conditions as {x,y,vx,vy}
    
    vector <double> IC1 = {1,1,0,0};
    vector <double> IC2 = {0,0,0,0};
    vector <double> IC3 = {0,1,0,0};
    
    vector <double> IC(IC1);
    IC.insert(IC.end(),IC2.begin(),IC2.end());
    IC.insert(IC.end(),IC3.begin(),IC3.end());
    
    // Instance the Three_body object
    
    Three_body TB; 
    
    TB.set_Masses(M);
    TB.set_Radii(R);
    TB.set_Initial_Conditions(IC);
    
    
    cout << TB.dy1() << endl;

    
    return 0;
}