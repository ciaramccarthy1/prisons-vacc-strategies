// Generated by covidm Thu Jul 22 10:04:57 2021

#include "user_defined.h"
#include "convenience.h"

void CppChanges(const vector<double>& x, Parameters& P)
{
    (void)x; (void)P;
    
}

double CppLogLikelihood(const vector<double>& x, Reporter& dyn)
{
    (void)dyn;
    double ll = 0.0;
    
    return ll;
}

bool CppObserver(Parameters& P, Randomizer& R, Reporter& dyn, double t, vector<double>& x)
{
    (void) P; (void) R; (void) dyn; (void) t;
    
    return true;
}
