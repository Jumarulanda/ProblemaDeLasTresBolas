#include "eu_int.h"

vector<vector<double>> eu_int :: eu_int_step(vector<double> (*F)(vector<double>), double time_step, vector<vector<double>> State) {
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
