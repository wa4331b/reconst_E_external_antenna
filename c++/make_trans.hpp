#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <boost/math/special_functions/hankel.hpp>

using namespace std;

void make_trans(int Nobs, double r_obs[], complex<double> k,
                int L, int Ndirection, double thetas[], double phis[], complex<double> alpha[]);


void make_ff(int Npol, double polAngle[],
             int Ntheta0, int Nphi0, double theta0[], double phi0[],
             complex<double> ff0_th[], complex<double> ff0_ph[],
             int Ndirection, double thetas[], double phis[],
             complex<double> ff_th[], complex<double> ff_ph[]);
