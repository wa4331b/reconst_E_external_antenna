#include <iostream>
#include <iomanip>
#include <complex>
#include "poynting.hpp"
#include "triangle_int.hpp"

void poyntingMatrix(int N_RWG,
                    double RWGNumber_trianglesCoord_tmp[],
                    int T,
                    int triangle_NRWG[],
                    int triangle_RWGNumber_tmp[],
                    int triangle_signInRWG_tmp[],
                    int triangle_surfaces[],
                    int N_target_surface,
                    int target_surface_tmp[],
                    complex<double> poyntingMatrix_tmp[])
{
	typedef Array<double, Dynamic, Dynamic, RowMajor> ArrayXXd_RowMajor;
	typedef Array<int, Dynamic, Dynamic, RowMajor> ArrayXXi_RowMajor;

	ArrayXi target_surface = Map<ArrayXi>(target_surface_tmp, N_target_surface);
	ArrayXXd RWGNumber_trianglesCoord
		= Map<ArrayXXd_RowMajor>(RWGNumber_trianglesCoord_tmp, N_RWG, 12);
	ArrayXXi triangle_RWGNumber = Map<ArrayXXi_RowMajor>(triangle_RWGNumber_tmp, T, 3);
	ArrayXXi triangle_signInRWG = Map<ArrayXXi_RowMajor>(triangle_signInRWG_tmp, T, 3);
	ArrayXXcd poyntingMatrix = ArrayXXcd::Zero(N_RWG, N_RWG);

	ArrayXd RWGLength(N_RWG);
	for (int RWGNumber = 0; RWGNumber < N_RWG; RWGNumber++)
	{
		double r1[3], r2[3];
		for (int i = 0; i < 3; i++)
		{
			r1[i] = RWGNumber_trianglesCoord(RWGNumber, i + 3);
			r2[i] = RWGNumber_trianglesCoord(RWGNumber, i + 6);
		}
		double r1_r2[3] = {r1[0] - r2[0], r1[1] - r2[1], r1[2] - r2[2]};
		RWGLength(RWGNumber) = sqrt(dot3D(r1_r2, r1_r2));
	}

	for (int tr = 0; tr < T; tr++)
	{
		if ((target_surface != triangle_surfaces[tr]).all()) continue;

		double r0[3], r1[3], r2[3];
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
		// for (int i=0; i<3; i++) triangle.n_hat[i] *= -1;

		double IT_r_square; // n_hat_X_r_p_dot_IT_r;
		double IT_r[3], IT_n_hat_X_r[3];
		IT_fm_fn(IT_r_square, IT_r, triangle);
		cross3D(IT_n_hat_X_r, triangle.n_hat, IT_r);

		for (int p = 0; p < triangle_NRWG[tr]; p++)
		{
			int RWGNumber_p = triangle_RWGNumber(tr, p);
			double sign_edge_p = triangle_signInRWG(tr, p);
			complex<double> C_p = sign_edge_p * RWGLength(RWGNumber_p) * 0.5 / triangle.A;

			double r_p[3], n_hat_X_r_p[3];
			if (triangle_signInRWG(tr, p) == 1)
			{
				for (int i = 0; i < 3; i++)
					r_p[i] = RWGNumber_trianglesCoord(RWGNumber_p, i);
			}
			else if (triangle_signInRWG(tr, p) == -1)
			{
				for (int i = 0; i < 3; i++)
					r_p[i] = RWGNumber_trianglesCoord(RWGNumber_p, i + 9);
			}
			else
			{
				cout << "Error at poynting.cpp. Terminated...";
				exit(1);
			}

			cross3D(n_hat_X_r_p, triangle.n_hat, r_p);

			const double n_hat_X_rp_dot_IT_r(
				n_hat_X_r_p[0] * IT_r[0] + n_hat_X_r_p[1] * IT_r[1] + n_hat_X_r_p[2] * IT_r[2]
			);

			for (int q = 0; q < triangle_NRWG[tr]; q++)
			{
				int RWGNumber_q = triangle_RWGNumber(tr, q);
				double sign_edge_q = triangle_signInRWG(tr, q);
				complex<double> C_pq = C_p * sign_edge_q * RWGLength(RWGNumber_q) * 0.5 / triangle.A;

				double r_q[3];
				if (triangle_signInRWG(tr, q) == 1)
				{
					for (int i = 0; i < 3; i++)
						r_q[i] = RWGNumber_trianglesCoord(RWGNumber_q, i);
				}
				else if (triangle_signInRWG(tr, q) == -1)
				{
					for (int i = 0; i < 3; i++)
						r_q[i] = RWGNumber_trianglesCoord(RWGNumber_q, i + 9);
				}
				else
				{
					cout << "Error at poynting.cpp. Terminated...";
					exit(1);
				}

				poyntingMatrix(RWGNumber_q, RWGNumber_p) += -C_pq * (
					n_hat_X_rp_dot_IT_r
					+ (r_q[0] * IT_n_hat_X_r[0]
						+ r_q[1] * IT_n_hat_X_r[1] + r_q[2] * IT_n_hat_X_r[2])
					- (r_q[0] * n_hat_X_r_p[0]
						+ r_q[1] * n_hat_X_r_p[1] + r_q[2] * n_hat_X_r_p[2]) * triangle.A
				);
			}
		}
	}

	int index = 0;
	for (int i = 0; i < N_RWG; i++)
	{
		for (int j = 0; j < N_RWG; j++)
		{
			poyntingMatrix_tmp[index] = poyntingMatrix(i, j);
			index++;
		}
	}
}
