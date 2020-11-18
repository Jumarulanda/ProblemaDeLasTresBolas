#include "Three_body.h"

Three_body::Three_body(){}

void Three_body::set_Masses(vector<double> M)
{
    // 3-vector object
    Masses = M;
    TM = M[0] + M[1] + M[2];
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
    
    vector<double> V1(&IC[0],&IC[4]);   //{x1,y1,vx1,vy1}
    vector<double> V2(&IC[4],&IC[8]);   //{x2,y2,vx2,vy2}
    vector<double> V3(&IC[8],&IC[12]);  //{x3,y3,vx3,vy3}
    
    // Filling the relative coordinate system vectors
    
    for (int i = 0; i < 4 ; i++)
    {
        s1.push_back(V2[i] - V3[i]);
        s2.push_back(V3[i] - V1[i]);
        s3.push_back(V1[i] - V2[i]);
    }
    
    update_distances();
}

vector <double> Three_body::g (vector<double> v1,vector<double> v2 , vector<double> v3)
{
    
    double gx = G* ( v1[0]/pow(S1,3) + v2[0]/pow(S2,3) + v3[0]/pow(S3,3));
    double gy = G* ( v1[1]/pow(S1,3) + v2[1]/pow(S2,3) + v3[1]/pow(S3,3));
    
    return {gx,gy};
}

// RHS of the differential equations

vector<double> F1()
{
    double fx = -TM*G*s1[0]/pow(S1,3) + Masses[0]*g(s1,s2,s3)[0];
    double fy = -TM*G*s1[1]/pow(S1,3) + Masses[0]*g(s1,s2,s3)[1];
    
    return {fx,fy};
}

vector<double> F2()
{
    double fx = -TM*G*s2[0]/pow(S2,3) + Masses[1]*g(s1,s2,s3)[0];
    double fy = -TM*G*s2[1]/pow(S2,3) + Masses[1]*g(s1,s2,s3)[1];
    
    return {fx,fy};
}

vector<double> F2()
{
    double fx = -TM*G*s3[0]/pow(S3,3) + Masses[2]*g(s1,s2,s3)[0];
    double fy = -TM*G*s3[1]/pow(S3,3) + Masses[2]*g(s1,s2,s3)[1];
    
    return {fx,fy};
}


void Three_body::update_distances()
{   
    // Distance magnitudes between the bodies
    S1 = sqrt(pow(s1[0],2) + pow(s1[1],2));
    S2 = sqrt(pow(s2[0],2) + pow(s2[1],2));
    S3 = sqrt(pow(s3[0],2) + pow(s3[1],2));
        
}





