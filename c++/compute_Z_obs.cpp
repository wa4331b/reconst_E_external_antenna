#include "compute_Z_obs.hpp"

void G_EJ_G_HJ(vector< vector< complex<double> > >& G_EJ,
               vector< vector< complex<double> > >& G_HJ,
               const double r_obs[],
               const double r_dip[],
               const complex<double>& eps,
               const complex<double>& mu,
               const complex<double>& k)
/**
* This function computes the homogeneous space Green's functions
* due to an elementary electric current element.
*
* By reciprocity, we have that G_EM = -G_HJ and G_HM = eps/mu * G_EJ.
*/
{
	const double r_obs_r_dip[3] = {r_obs[0] - r_dip[0], r_obs[1] - r_dip[1], r_obs[2] - r_dip[2]};
	const double R = sqrt(dot3D(r_obs_r_dip, r_obs_r_dip));
	const complex<double> kRsquare = k * k * R * R;
	const complex<double> term_1 = 1.0 + 1.0 / (I * k * R);
	const complex<double> term_2 = 1.0 / R * term_1 + I * k / 2.0 * (term_1 - 1.0 / kRsquare);
	// JPA : tentative pour reduire le temps CPU lorsqu'on a beaucoup de dipoles en excitation
	// On decompose la ligne suivante : const complex<double> exp_ikR = exp(-I*k*R);
	const complex<double> minus_I_k_R = -I * k * R;
	double c, s, e;
	e = (minus_I_k_R.real() == 0.0) ? 1.0 : exp(minus_I_k_R.real());
	s = sin(minus_I_k_R.imag());
	c = cos(minus_I_k_R.imag());
	//  sincos(minus_I_k_R.imag(), &s, &c);
	const complex<double> exp_ikR = e * complex<double>(c, s);
	// fin JPA
	const complex<double> exp_ikR_R = sqrt(mu / eps) / (2.0 * M_PI) * exp_ikR / R;
	const double x_xp = r_obs_r_dip[0], y_yp = r_obs_r_dip[1], z_zp = r_obs_r_dip[2];
	const double ONE_R_R = 1.0 / (R * R);
	const double x_xp_R_square = (x_xp * x_xp) * ONE_R_R;
	const double y_yp_R_square = (y_yp * y_yp) * ONE_R_R;
	const double z_zp_R_square = (z_zp * z_zp) * ONE_R_R;
	G_EJ[0][1] = term_2 * y_yp / R * x_xp / R * exp_ikR_R;
	G_EJ[0][2] = term_2 * z_zp / R * x_xp / R * exp_ikR_R;
	G_EJ[1][2] = term_2 * z_zp / R * y_yp / R * exp_ikR_R;
	G_EJ[1][0] = G_EJ[0][1];
	G_EJ[2][0] = G_EJ[0][2];
	G_EJ[2][1] = G_EJ[1][2];
	G_EJ[0][0] = (x_xp_R_square * term_1 / R - (1.0 - x_xp_R_square) * I * k / 2.0 * (term_1 - 1.0 / kRsquare)) * exp_ikR_R;
	G_EJ[1][1] = (y_yp_R_square * term_1 / R - (1.0 - y_yp_R_square) * I * k / 2.0 * (term_1 - 1.0 / kRsquare)) * exp_ikR_R;
	G_EJ[2][2] = (z_zp_R_square * term_1 / R - (1.0 - z_zp_R_square) * I * k / 2.0 * (term_1 - 1.0 / kRsquare)) * exp_ikR_R;

	complex<double> G_i = exp_ikR / (4.0 * M_PI) * (1.0 + I * k * R) / (R * R * R);
	G_HJ[0][1] = (z_zp) * G_i;
	G_HJ[1][0] = -G_HJ[0][1];
	G_HJ[2][0] = (y_yp) * G_i;
	G_HJ[2][1] = -1.0 * (x_xp) * G_i;
	G_HJ[0][2] = -G_HJ[2][0];
	G_HJ[1][2] = -G_HJ[2][1];
	for (int i = 0; i < 3; i++) G_HJ[i][i] = 0.0;
}

/*****************************************
 * computation of the observation fields *
 *****************************************/
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
                   const ArrayXi& target_surface)
{
	// def of k, mu_i, eps_i
	int T = triangle_NRWG.size();
	complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r, k = w * sqrt(eps * mu);
	complex<double> imp2 = mu / eps;

	// geometrical entities
	vector< vector< complex<double> > > G_EJ, G_HJ;
	G_EJ.resize(3);
	G_HJ.resize(3);
	for (int i = 0; i < 3; i++)
	{
		G_EJ[i].resize(3);
		G_HJ[i].resize(3);
	}
	double rObs[3];
	for (int i = 0; i < 3; ++i) rObs[i] = r_obs(i);

	double r0[3], r1[3], r2[3];

	Z_EJ = ArrayXXcd::Zero(3, N_RWG_test);
	Z_EM = ArrayXXcd::Zero(3, N_RWG_test);
	Z_HJ = ArrayXXcd::Zero(3, N_RWG_test);
	Z_HM = ArrayXXcd::Zero(3, N_RWG_test);

	for (int tr = 0; tr < T; tr++)
	{
		if ((target_surface != triangle_surfaces(tr)).all()) continue;
		if (triangle_signInRWG(tr, 0) == 1)
		{
			for (int i = 0; i < 3; i++)
			{
				int rwg = triangle_RWGNumber(tr, 0);
				r0[i] = RWGNumber_trianglesCoord(rwg, i);
				r1[i] = RWGNumber_trianglesCoord(rwg, i + 3);
				r2[i] = RWGNumber_trianglesCoord(rwg, i + 6);
			}
		}
		else if (triangle_signInRWG(tr, 0) == -1)
		{
			for (int i = 0; i < 3; i++)
			{
				int rwg = triangle_RWGNumber(tr, 0);
				r0[i] = RWGNumber_trianglesCoord(rwg, i + 6);
				r1[i] = RWGNumber_trianglesCoord(rwg, i + 3);
				r2[i] = RWGNumber_trianglesCoord(rwg, i + 9);
			}
		}

		Triangle triangle(r0, r1, r2, 0);

		const double rGrav_rObs[3] = {triangle.r_grav[0] - rObs[0],
			triangle.r_grav[1] - rObs[1],
			triangle.r_grav[2] - rObs[2]
		};

		double R_os = sqrt(dot3D(rGrav_rObs, rGrav_rObs));
		bool IS_NEAR = (R_os - 1.5 * triangle.R_max <= 0.0);
		double sum_weigths = 1.;
		const double *xi, *eta, *weights;
		int N_points;
		if (IS_NEAR)
		{
			N_points = IT_near.Nquads;
			xi = IT_near.xi;
			eta = IT_near.eta;
			weights = IT_near.weights;
		}
		else
		{
			N_points = IT_far.Nquads;
			xi = IT_far.xi;
			eta = IT_far.eta;
			weights = IT_far.weights;
		}
		//		IT_points(xi, eta, weights, sum_weigths, N_points);

		const double norm_factor = triangle.A / sum_weigths;

		for (int j = 0; j < N_points; ++j)
		{
			double r_src[3];
			r_src[0] = r0[0] * xi[j] + r1[0] * eta[j] + r2[0] * (1 - xi[j] - eta[j]);
			r_src[1] = r0[1] * xi[j] + r1[1] * eta[j] + r2[1] * (1 - xi[j] - eta[j]);
			r_src[2] = r0[2] * xi[j] + r1[2] * eta[j] + r2[2] * (1 - xi[j] - eta[j]);
			G_EJ_G_HJ(G_EJ, G_HJ, rObs, r_src, eps, mu, k);

			for (int rwg = 0; rwg < triangle_NRWG(tr); rwg++)
			{
				int RWGNumber = triangle_RWGNumber(tr, rwg);
				double r_opp[3];
				if (triangle_signInRWG(tr, rwg) == 1)
				{
					for (int i = 0; i < 3; i++)
						r_opp[i] = RWGNumber_trianglesCoord(RWGNumber, i);
				}
				else if (triangle_signInRWG(tr, rwg) == -1)
				{
					for (int i = 0; i < 3; i++)
						r_opp[i] = RWGNumber_trianglesCoord(RWGNumber, i + 9);
				}
				else
				{
					cout << "error in compute_Z_obs.hpp: triangle_signInRWG" << endl;
					exit(1);
				}
				// computation of E_obs
				const double fn[3]
					= {r_src[0] - r_opp[0], r_src[1] - r_opp[1], r_src[2] - r_opp[2]};

				complex<double> ITo_E_obs[3] = {0.0, 0.0, 0.0};
				complex<double> ITo_H_obs[3] = {0.0, 0.0, 0.0};
				for (int m = 0; m < 3; m++)
					ITo_E_obs[m]
						+= (G_EJ[m][0] * fn[0] + G_EJ[m][1] * fn[1] + G_EJ[m][2] * fn[2])
						* weights[j];

				// computation of H_obs
				for (int m = 0; m < 3; m++)
				{
					ITo_H_obs[m]
						+= (G_HJ[m][0] * fn[0] + G_HJ[m][1] * fn[1] + G_HJ[m][2] * fn[2])
						* weights[j];
				}

				for (int i = 0; i < 3; ++i) ITo_E_obs[i] *= norm_factor;
				for (int i = 0; i < 3; ++i) ITo_H_obs[i] *= norm_factor;

				double C_rp = triangle_signInRWG(tr, rwg) * RWGLength(RWGNumber)
					* 0.5 / triangle.A;

				for (int i = 0; i < 3; ++i)
				{
					Z_EJ(i, RWGNumber) += ITo_E_obs[i] * C_rp;
					Z_HJ(i, RWGNumber) += ITo_H_obs[i] * C_rp;
					Z_EM(i, RWGNumber) -= ITo_H_obs[i] * C_rp;
					Z_HM(i, RWGNumber) += ITo_E_obs[i] * C_rp / imp2;
				}
			}
		}
	}
}

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
                        int target_surface_tmp[])
{
	typedef Array<double, Dynamic, Dynamic, RowMajor> ArrayXXd_RowMajor;
	typedef Array<int, Dynamic, Dynamic, RowMajor> ArrayXXi_RowMajor;

	ArrayXXd r_obs = Map<ArrayXXd_RowMajor>(r_obs_tmp, N_obs, 3);
	ArrayXi target_surface = Map<ArrayXi>(target_surface_tmp, N_target_surface);

	ArrayXXd RWGNumber_trianglesCoord
		= Map<ArrayXXd_RowMajor>(RWGNumber_trianglesCoord_tmp, N_RWG, 12);

	ArrayXi triangle_NRWG = Map<ArrayXi>(triangle_NRWG_tmp, T);
	ArrayXXi triangle_RWGNumber
		= Map<ArrayXXi_RowMajor>(triangle_RWGNumber_tmp, T, 3);
	ArrayXXi triangle_signInRWG
		= Map<ArrayXXi_RowMajor>(triangle_signInRWG_tmp, T, 3);
	ArrayXi triangle_surfaces = Map<ArrayXi>(triangle_surfaces_tmp, T);

	ArrayXd RWGLength(N_RWG);
	double r1[3], r2[3];
	for (int RWGNumber = 0; RWGNumber < N_RWG; RWGNumber++)
	{
		for (int i = 0; i < 3; i++)
		{
			r1[i] = RWGNumber_trianglesCoord(RWGNumber, i + 3);
			r2[i] = RWGNumber_trianglesCoord(RWGNumber, i + 6);
		}
		double r1_r2[3] = {r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2]};
		RWGLength(RWGNumber) = sqrt(dot3D(r1_r2, r1_r2));
	}

	IT_points_weights IT_near(25), IT_far(6);
	int index = 0;
	for (int i = 0; i < N_obs; i++)
	{
		ArrayXXcd Z_EJ(3, N_RWG), Z_EM(3, N_RWG), Z_HJ(3, N_RWG), Z_HM(3, N_RWG);
		Array3d r_obs_row = r_obs.row(i);
		compute_Z_obs(Z_EJ, Z_EM, Z_HJ, Z_HM,
		              r_obs_row, N_RWG,
		              RWGNumber_trianglesCoord, RWGLength,
		              triangle_NRWG, triangle_RWGNumber,
		              triangle_signInRWG, triangle_surfaces,
		              IT_near, IT_far, w, eps_r, mu_r, target_surface);

		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < N_RWG; k++)
			{
				Z_EJ_obs_tmp[index] = Z_EJ(j, k);
				Z_EM_obs_tmp[index] = Z_EM(j, k);
				Z_HJ_obs_tmp[index] = Z_HJ(j, k);
				Z_HM_obs_tmp[index] = Z_HM(j, k);
				index++;
			}
		}
	}
}
