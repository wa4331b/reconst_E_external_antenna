#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void compute_FFmat(complex<double> FF_J_th[], complex<double> FF_J_ph[],
                   complex<double> FF_M_th[], complex<double> FF_M_ph[],
                   int N_direction, double thetas[], double phis[],
                   int N_RWG,
                   double RWGNumber_trianglesCoord_tmp[],
                   int T,
                   int triangle_NRWG_tmp[],
                   int triangle_RWGNumber_tmp[],
                   int triangle_signInRWG_tmp[],
                   int triangle_surfaces_tmp[],
                   double w, complex<double> eps_r,
                   complex<double> mu_r);
