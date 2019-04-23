#include <iostream>
#include <complex>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

#include "EMConstants.hpp"
#include "GK_triangle.hpp"
#include "triangle_int.hpp"
#include "dictionary.hpp"

void G_EJ_G_HJ(vector< vector< complex<double> > >& G_EJ,
               vector< vector< complex<double> > >& G_HJ,
               const double r_obs[],
               const double r_dip[],
               const complex<double>& eps,
               const complex<double>& mu,
               const complex<double>& k);

void compute_Z_obs(ArrayXXcd& Z_EJ, ArrayXXcd& Z_EM,
                   ArrayXXcd& Z_HJ, ArrayXXcd& Z_HM,
                   const Array3d& r_obs,
                   int N_RWG_test,
                   const ArrayXXd& RWGNumber_trianglesCoord,
                   const ArrayXd& RWGLength,
                   const ArrayXi& triangle_NRWG,
                   const ArrayXXi& triangle_RWGNumber,
                   const ArrayXXi& triangle_signInRWG,
                   const ArrayXi& triangle_surfaces,
                   IT_points_weights& IT_near,
                   IT_points_weights& IT_far,
                   const double w,
                   const complex<double>& eps_r,
                   const complex<double>& mu_r,
                   const ArrayXi& target_surface);

void wrap_compute_Z_obs(complex<double> Z_EJ_obs_tmp[],
                        complex<double> Z_EM_obs_tmp[],
                        complex<double> Z_HJ_obs_tmp[],
                        complex<double> Z_HM_obs_tmp[],
                        int N_obs,
                        double r_obs_tmp[],
                        int N_RWG,
                        double RWGNumber_trianglesCoord_tmp[],
                        int T,
                        int triangle_NRWG_tmp[],
                        int triangle_RWGNumber_tmp[],
                        int triangle_signInRWG_tmp[],
                        int triangle_surfaces_tmp[],
                        double w, complex<double> eps_r,
                        complex<double> mu_r,
                        int N_target_surface,
                        int target_surface_tmp[]);
