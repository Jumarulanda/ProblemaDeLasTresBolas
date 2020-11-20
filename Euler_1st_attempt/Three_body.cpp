#include "Three_body.h"

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
    
    vector<double> V1(&IC[0],&IC[4]);   //{x1,y1,vx1,vy1}
    vector<double> V2(&IC[4],&IC[8]);   //{x2,y2,vx2,vy2}
    vector<double> V3(&IC[8],&IC[12]);  //{x3,y3,vx3,vy3}
    
    // Determing center of mass coordinates
    
    double x_CM = (Masses[1]*V1[0] + Masses[2]*V2[0] + Masses[3]*V3[0])/Masses[0];
    double y_CM = (Masses[1]*V1[1] + Masses[2]*V2[1] + Masses[3]*V3[1])/Masses[0];
    
    CM = {x_CM , y_CM};
    
    // Filling the relative coordinate system vectors
    
    for (int i = 0; i < 4 ; i++)
    {
        s1.push_back(V2[i] - V3[i]);
        s2.push_back(V3[i] - V1[i]);
        s3.push_back(V1[i] - V2[i]);
    }
    
    update(0);
}

vector <double> Three_body::g (vector<double> v1,vector<double> v2 , vector<double> v3)
{
    double gx = v1[0]/pow(S1,3) + v2[0]/pow(S2,3) + v3[0]/pow(S3,3);
    double gy = v1[1]/pow(S1,3) + v2[1]/pow(S2,3) + v3[1]/pow(S3,3);
    
    return {gx,gy};
}


void Three_body::update(double time_step)
{   
    // Distance magnitudes between the bodies
    
    S1 = sqrt(pow(s1[0],2) + pow(s1[1],2));
    S2 = sqrt(pow(s2[0],2) + pow(s2[1],2));
    S3 = sqrt(pow(s3[0],2) + pow(s3[1],2));
    
    //Forces over the relative coordinate system
    
    double fx1 = G*( -Masses[0]*s1[0]/pow(S1,3) + Masses[1]*g(s1,s2,s3)[0]);
    double fy1 = G*( -Masses[0]*s1[1]/pow(S1,3) + Masses[1]*g(s1,s2,s3)[1]);
    
    F1 = {fx1,fy1};
    
    double fx2 = G*( -Masses[0]*s2[0]/pow(S2,3) + Masses[2]*g(s1,s2,s3)[0]);
    double fy2 = G*( -Masses[0]*s2[1]/pow(S2,3) + Masses[2]*g(s1,s2,s3)[1]);
    
    F2 = {fx2,fy2};
    
    double fx3 = G*( -Masses[0]*s3[0]/pow(S3,3) + Masses[3]*g(s1,s2,s3)[0]);
    double fy3 = G*( -Masses[0]*s3[1]/pow(S3,3) + Masses[3]*g(s1,s2,s3)[1]);
    
    F3 =  {fx3,fy3};
    
    time += time_step ;
        
}

void Three_body::euler_integrator(double time_step)
{
    double h = time_step;
    
    // Para el cuerpo 1
    
    s1[2] += h*F1[0];   // vx1 -> vx1 + h*F1[x]
    s1[3] += h*F1[1];   // vy1 -> vy1 + h*F1[y]
    
    s1[0] += h*s1[2];   // x1 -> x1 + h*vx1
    s1[1] += h*s1[3];   // y1 -> y1 + h*vy1
        
    // Para el cuerpo 2
        
    s2[2] += h*F2[0];   // vx1 -> vx1 + h*F1[x]
    s2[3] += h*F2[1];   // vy1 -> vy1 + h*F1[y]
    
    s2[0] += h*s2[2];   // x1 -> x1 + h*vx1
    s2[1] += h*s2[3];   // y1 -> y1 + h*vy1
        
    // Para el cuerpo 3
        
    s3[2] += h*F3[0];   // vx1 -> vx1 + h*F1[x]
    s3[3] += h*F3[1];   // vy1 -> vy1 + h*F1[y]
    
    s3[0] += h*s3[2];   // x1 -> x1 + h*vx1
    s3[1] += h*s3[3];   // y1 -> y1 + h*vy1
    
    update(time_step);
    
}

void Three_body::write_on_file(ofstream *file)
{
    /*
    Recoverging the absolute coordinate system:
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
    
    double x1 = (-Masses[3]*s2[0] + Masses[2]*s3[0])/Masses[0] + CM[0];
    double y1 = (-Masses[3]*s2[1] + Masses[2]*s3[1])/Masses[0] + CM[1];
        
    double x2 = (-Masses[3]*s2[0] - (Masses[1] + Masses[3])*s3[0])/Masses[0] + CM[0];
    double y2 = (-Masses[3]*s2[1] - (Masses[1] + Masses[3])*s3[1])/Masses[0] + CM[1];
    
    double x3 = ((Masses[1] + Masses[2])*s2[0] + Masses[2]*s3[0])/Masses[0] + CM[0];
    double y3 = ((Masses[1] + Masses[2])*s2[1] + Masses[2]*s3[1])/Masses[0] + CM[1];
    
    if (file -> is_open())
    {
	  *file << x1 << "," << y1 << "," << x2 << "," << y2 << "," << x3 << "," << y3 << "," << time << "\n";
    }

}

void Three_body::evol_system(double number_of_steps, double time_step, ofstream *file)
{
    // Writing file header
    if (file -> is_open())
    {
	  *file << Masses[0] << "," << Masses[1] << "," << Masses[2] << "," << Radii[0] << "," << Radii[1] << "," << Radii[2] << "," << 0 << "\n";
    }
    
    
    for(double i = 0 ; i <= number_of_steps ; i ++)
    {
        write_on_file(file);
        euler_integrator(time_step);
        if ( check_radii()){ break;}
    }
    
}

bool Three_body::check_radii()
{
    bool sp12 = Radii[0] + Radii[1] > S3;  // R1 + R2 > S3
    bool sp13 = Radii[0] + Radii[2] > S2;  // R1 + R3 > S2
    bool sp23 = Radii[1] + Radii[2] > S1;  // R2 + R3 > S1
    
    if (sp12 || sp13 || sp23)
    {
        cout << "There was a colliton" << endl;
    }
    return sp12 || sp13 || sp23;
}







