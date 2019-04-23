#include <iostream>
#include <complex>
#include <vector>
using namespace std;

complex<double> ff_interpolation2D(int Ntheta0, int Nphi0, double theta0[], double phi0[],
                                   complex<double> ff0[],
                                   double theta1, double phi1);