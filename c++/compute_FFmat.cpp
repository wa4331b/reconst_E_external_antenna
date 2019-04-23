#include "compute_FFmat.hpp"
#include "EMConstants.hpp"
#include "GK_triangle.hpp"
#include "triangle_int.hpp"

void compute_FFrow(double theta, double phi,
                   ArrayXXd& RWGNumber_trianglesCoord,
                   ArrayXd& RWGLength,
                   ArrayXi& triangle_NRWG,
                   ArrayXXi& triangle_RWGNumber,
                   ArrayXXi& triangle_signInRWG,
                   ArrayXi& triangle_surfaces,
                   complex<double> FF_J_th[], complex<double> FF_J_ph[],
                   complex<double> FF_M_th[], complex<double> FF_M_ph[],
                   double w,
                   complex<double>& eps_r,
                   complex<double>& mu_r,
                   int FULL_PRECISION)
{
	int T = triangle_NRWG.size();
	complex<double> I(0., 1.);

	complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps * mu);
	complex<double> I_k = I * k;
	complex<double> EJ_factor(-I * mu_0 * w * mu_r);
	complex<double> EM_factor(-I * k);

	Array3d theta_hat, phi_hat;
	theta_hat << cos(theta) * cos(phi) , cos(theta) * sin(phi) , -sin(theta);
	phi_hat << -sin(phi) , cos(phi) , 0.0;

	int N_points;
	if (FULL_PRECISION != 0)
		N_points = 6;
	else
		N_points = 1;
	IT_points_weights it(N_points);

	double k_hat[3] = {sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)};

	for (int tr = 0; tr < T; tr++) {
		double r0[3], r1[3], r2[3];
		if (triangle_signInRWG(tr, 0) == 1) {
			for (int i = 0; i < 3; i++) {
				int rwg = triangle_RWGNumber(tr, 0);
				r0[i] = RWGNumber_trianglesCoord(rwg, i);
				r1[i] = RWGNumber_trianglesCoord(rwg, i + 3);
				r2[i] = RWGNumber_trianglesCoord(rwg, i + 6);
			}
		} else if (triangle_signInRWG(tr, 0) == -1) {
			for (int i = 0; i < 3; i++) {
				int rwg = triangle_RWGNumber(tr, 0);
				r0[i] = RWGNumber_trianglesCoord(rwg, i + 6);
				r1[i] = RWGNumber_trianglesCoord(rwg, i + 3);
				r2[i] = RWGNumber_trianglesCoord(rwg, i + 9);
			}
		}

		Triangle triangle(r0, r1, r2, 0);
		const double *xi, *eta, *weights;
		xi = it.xi; eta = it.eta; weights = it.weights;

		double norm_factor = triangle.A;

		for (int j = 0; j < N_points; j++) {
			double r_src[3];
			r_src[0] = r0[0] * xi[j] + r1[0] * eta[j] + r2[0] * (1 - xi[j] - eta[j]);
			r_src[1] = r0[1] * xi[j] + r1[1] * eta[j] + r2[1] * (1 - xi[j] - eta[j]);
			r_src[2] = r0[2] * xi[j] + r1[2] * eta[j] + r2[2] * (1 - xi[j] - eta[j]);
			complex<double> G_FF = exp(I_k * (r_src[0] * k_hat[0] + r_src[1] * k_hat[1] + r_src[2] * k_hat[2]));

			for (int rwg = 0; rwg < triangle_NRWG(tr); rwg++) {
				int RWGNumber = triangle_RWGNumber(tr, rwg);
				double r_opp[3];
				if (triangle_signInRWG(tr, rwg) == 1) {
					for (int i = 0; i < 3; i++)
						r_opp[i] = RWGNumber_trianglesCoord(RWGNumber, i);
				} else if (triangle_signInRWG(tr, rwg) == -1) {
					for (int i = 0; i < 3; i++)
						r_opp[i] = RWGNumber_trianglesCoord(RWGNumber, i + 9);
				} else {
					cout << "error in compute_FF.hpp: triangle_signInRWG" << endl
					     << "tr = " << tr << ", rwg = " << rwg << endl
					     << triangle_NRWG(tr) << endl;
					exit(1);
				}
				const double fn[3] = {r_src[0] - r_opp[0], r_src[1] - r_opp[1], r_src[2] - r_opp[2]};
				complex<double> ITo_E_obs[3] = {fn[0] * G_FF, fn[1] * G_FF, fn[2] * G_FF};
				double C_rp = triangle_signInRWG(tr, rwg) * RWGLength(RWGNumber) * 0.5 / triangle.A;

				for (int m = 0; m < 3; m++)
					ITo_E_obs[m] *= weights[j] * norm_factor * C_rp;

				Array3cd FF_J, FF_M;
				for (int i = 0; i < 3; i++) {
					FF_J_th[RWGNumber] += EJ_factor * (ITo_E_obs[i] * theta_hat(i));
					FF_J_ph[RWGNumber] += EJ_factor * (ITo_E_obs[i] * phi_hat(i));
					FF_M_th[RWGNumber] += EM_factor * (ITo_E_obs[i] * phi_hat(i));
					FF_M_ph[RWGNumber] -= EM_factor * (ITo_E_obs[i] * theta_hat(i));
				}
			}
		}
	}
}

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
                   complex<double> mu_r)
{
	typedef Array<double, Dynamic, Dynamic, RowMajor> ArrayXXd_RowMajor;
	typedef Array<int, Dynamic, Dynamic, RowMajor> ArrayXXi_RowMajor;

	ArrayXXd RWGNumber_trianglesCoord = Map<ArrayXXd_RowMajor>(RWGNumber_trianglesCoord_tmp, N_RWG, 12);

	ArrayXi triangle_NRWG = Map<ArrayXi>(triangle_NRWG_tmp, T);
	ArrayXXi triangle_RWGNumber = Map<ArrayXXi_RowMajor>(triangle_RWGNumber_tmp, T, 3);
	ArrayXXi triangle_signInRWG = Map<ArrayXXi_RowMajor>(triangle_signInRWG_tmp, T, 3);
	ArrayXi triangle_surfaces = Map<ArrayXi>(triangle_surfaces_tmp, T);

	ArrayXd RWGLength(N_RWG);
	double r1[3], r2[3];
	for (int RWGNumber = 0; RWGNumber < N_RWG; RWGNumber++) {
		for (int i = 0; i < 3; i++) {
			r1[i] = RWGNumber_trianglesCoord(RWGNumber, i + 3);
			r2[i] = RWGNumber_trianglesCoord(RWGNumber, i + 6);
		}
		double r1_r2[3] = {r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2]};
		RWGLength(RWGNumber) = sqrt(dot3D(r1_r2, r1_r2));
	}

	int index = 0;
	for (int i = 0; i < N_direction; i++) {
		compute_FFrow(thetas[i], phis[i], RWGNumber_trianglesCoord, RWGLength,
		              triangle_NRWG, triangle_RWGNumber, triangle_signInRWG, triangle_surfaces,
		              &(FF_J_th[index]), &(FF_J_ph[index]), &(FF_M_th[index]), &(FF_M_ph[index]),
		              w, eps_r, mu_r, 1);
		index += N_RWG;
	}
}
