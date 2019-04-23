/**********************************************************************
 *
 * Z_EJ_Z_HJ_FS_triangles_arrays.cpp
 *
 * Copyright (C) 2014 Idesbald Van den Bosch
 *
 * This file is part of Puma-EM.
 *
 * Puma-EM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Puma-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Puma-EM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions/bugs : <vandenbosch.idesbald@gmail.com>
 *
 **********************************************************************/

#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include <algorithm>

#include "Z_EJ_Z_HJ.hpp"
#include "EMConstants.hpp"
#include "triangle_int.hpp"
#include "GK_triangle.hpp"
#include "dictionary.hpp"

void Z_CFIE_computation(ArrayXXcd& Z_CFIE_EJ,
                        ArrayXXcd& Z_CFIE_EM,
                        ArrayXXcd& Z_CFIE_HJ,
                        ArrayXXcd& Z_CFIE_HM,
                        const complex<double> CFIE[],
                        const bool CFIE_bool[], const bool M_current,
                        const double signSurfObs,
                        const double signSurfSrc,
                        const ArrayXi& numbers_RWG_src,
                        const ArrayXi& numbers_RWG_test,
                        const int srcRWGNumber_surface[],
                        const ArrayXXi& RWGNumber_signedTriangles,
                        const ArrayXXi& RWGNumber_nodes,
                        const ArrayXXd& nodesCoord,
                        const double w,
                        const complex<double>& eps_r,
                        const complex<double>& mu_r,
                        const ArrayXi& target_surface)
{
	// half RWGs construction
	int N_RWG_src = numbers_RWG_src.size(), N_RWG_test = numbers_RWG_test.size();
	vector<RWG> src_RWGs, test_RWGs;
	src_RWGs.reserve(N_RWG_src);
	test_RWGs.reserve(N_RWG_test);

	for (int i = 0; i < N_RWG_src; ++i)
	{
		int RWGNumber = numbers_RWG_src(i);

		int triangle_numbers[2], triangle_signs[2];
		triangle_numbers[0] = abs(RWGNumber_signedTriangles(RWGNumber, 0));
		triangle_numbers[1] = abs(RWGNumber_signedTriangles(RWGNumber, 1));
		triangle_signs[0] = 1;
		triangle_signs[1] = -1;

		// the nodes of the RWG
		const int n0 = RWGNumber_nodes(RWGNumber, 0);
		const int n1 = RWGNumber_nodes(RWGNumber, 1);
		const int n2 = RWGNumber_nodes(RWGNumber, 2);
		const int n3 = RWGNumber_nodes(RWGNumber, 3);
		const double r0[3]
			= {nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
		const double r1[3]
			= {nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
		const double r2[3]
			= {nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
		const double r3[3]
			= {nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};

		src_RWGs.push_back(RWG(RWGNumber, srcRWGNumber_surface[RWGNumber],
		                       triangle_numbers, triangle_signs, r0, r1, r2, r3));
	}

	for (int i = 0; i < N_RWG_test; ++i)
	{
		int RWGNumber = numbers_RWG_test(i);

		int triangle_numbers[2], triangle_signs[2];

		triangle_numbers[0] = abs(RWGNumber_signedTriangles(RWGNumber, 0));
		triangle_numbers[1] = abs(RWGNumber_signedTriangles(RWGNumber, 1));
		triangle_signs[0] = 1;
		triangle_signs[1] = -1;

		// the nodes of the RWG
		const int n0 = RWGNumber_nodes(RWGNumber, 0);
		const int n1 = RWGNumber_nodes(RWGNumber, 1);
		const int n2 = RWGNumber_nodes(RWGNumber, 2);
		const int n3 = RWGNumber_nodes(RWGNumber, 3);
		const double r0[3] = {nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
		const double r1[3] = {nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
		const double r2[3] = {nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
		const double r3[3] = {nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};

		test_RWGs.push_back(RWG(RWGNumber, srcRWGNumber_surface[RWGNumber],
		                        triangle_numbers, triangle_signs, r0, r1, r2, r3));
	}

	// triangles
	vector< Dictionary<int, int> > srcTriangleToRWG, testTriangleToRWG;
	srcTriangleToRWG.reserve(N_RWG_src * 2);
	testTriangleToRWG.reserve(N_RWG_test * 2);
	for (unsigned int i = 0; i < src_RWGs.size(); ++i)
	{
		srcTriangleToRWG.push_back(Dictionary<int, int>(src_RWGs[i].triangleNumbers[0], i));
		srcTriangleToRWG.push_back(Dictionary<int, int>(src_RWGs[i].triangleNumbers[1], i));
	}

	for (unsigned int i = 0; i < test_RWGs.size(); ++i)
	{
		testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers[0], i));
		testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers[1], i));
	}

	sort(srcTriangleToRWG.begin(), srcTriangleToRWG.end());
	sort(testTriangleToRWG.begin(), testTriangleToRWG.end());
	vector<Triangle> triangles_src, triangles_test;
	constructVectorTriangles(triangles_src, src_RWGs, srcTriangleToRWG);
	constructVectorTriangles(triangles_test, test_RWGs, testTriangleToRWG);

	// Z_CFIE computation
	// def of k, mu_i, eps_i
	int EXTRACT_1_R, EXTRACT_R;
	const complex<double> mu = mu_0 * mu_r, eps = eps_0 * eps_r;
	const complex<double> k = w * sqrt(eps * mu), k_square = k * k;
	const complex<double> z_pq_factor = 1.0 / (I * w * eps);

	Z_CFIE_EJ = ArrayXXcd::Zero(N_RWG_test, N_RWG_src);
	Z_CFIE_EM = ArrayXXcd::Zero(N_RWG_test, N_RWG_src);
	Z_CFIE_HJ = ArrayXXcd::Zero(N_RWG_test, N_RWG_src);
	Z_CFIE_HM = ArrayXXcd::Zero(N_RWG_test, N_RWG_src);

	// const bool tHM = false, nHM = false, tEM = true, nEM = false;
	bool tEJ = CFIE_bool[0], nEJ = CFIE_bool[1], tHJ = CFIE_bool[2], nHJ = CFIE_bool[3];
	bool tEM, nEM, tHM, nHM;
	if (M_current)
	{
		tEM = CFIE_bool[0];
		nEM = CFIE_bool[1];
		tHM = CFIE_bool[2];
		nHM = CFIE_bool[3];
	}
	else
	{
		tEM = false;
		nEM = false;
		tHM = false;
		nHM = false;
	}

	// possively 1, 3, 4, 6, 7, 12, 13, 16, 19, 25, 27, 33, 37, 42, 48, 52, 61, 70, 73, 79
	IT_points_weights IT_3(3), IT_6(6), IT_13(13);

	ArrayXd omega_eq(triangles_test.size());
	for (unsigned int r = 0; r < triangles_test.size(); ++r)
	{
		double* n_hat0 = triangles_test[r].n_hat;
		double* r_grav0 = triangles_test[r].r_grav;
		vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
		double denom = 0., numer = 0.;
		for (unsigned int p = 0; p < RWGsIndexes_test.size(); ++p)
		{
			const int index_p = RWGsIndexes_test[p];
			int triangleNumber1 = test_RWGs[index_p].triangleNumbers[0];
			if (triangleNumber1 == triangles_test[r].number)
				triangleNumber1 = test_RWGs[index_p].triangleNumbers[1];
			int triangleIndex1 = -1;
			double* n_hat1 = triangles_test[0].n_hat;
			double* r_grav1 = triangles_test[0].r_grav;
			for (unsigned int rp = 0; rp < triangles_test.size(); ++rp)
			{
				if (triangles_test[rp].number == triangleNumber1)
				{
					triangleIndex1 = rp;
					n_hat1 = triangles_test[rp].n_hat;
					r_grav1 = triangles_test[rp].r_grav;
					break;
				}
			}
			if (triangleIndex1 < 0)
			{
				cout << "Triangle index search error, terminated..." << endl;
				exit(1);
			}
			double n_hat_sub[3], r_grav_sub[3];
			for (int i = 0; i < 3; i++)
			{
				n_hat_sub[i] = n_hat1[i] - n_hat0[i];
				r_grav_sub[i] = r_grav1[i] - r_grav0[i];
			}
			double omega = M_PI - acos(dot3D(n_hat1, n_hat0));
			if (dot3D(r_grav_sub, n_hat_sub) >= 0.) omega = 2 * M_PI - omega;
			omega *= 2;

			denom += test_RWGs[index_p].length;
			numer += test_RWGs[index_p].length * omega;
		}
		omega_eq(r) = numer / denom;
	}

	// loop on the observation RWGs
	for (unsigned int r = 0; r < triangles_test.size(); ++r)
	{
		if ((target_surface != triangles_test[r].surface_num).all())
			continue;

		// computation of the triangle-to-triangle terms
		double* n_hat;
		n_hat = triangles_test[r].n_hat;
		double IT_r_square; // n_hat_X_r_p_dot_IT_r;
		double IT_r[3], IT_n_hat_X_r[3];

		// double solid_angle_factor = omega_eq(r) / (4 * M_PI) / (2 * M_PI);
		double solid_angle_factor = 1.;

		// serves for <f_m ; f_n> and <f_m ; n x f_n>
		IT_fm_fn(IT_r_square, IT_r, triangles_test[r]);
		// serves for <f_m ; n x f_n>
		cross3D(IT_n_hat_X_r, n_hat, IT_r);

		// the RWGs concerned by the test triangle
		vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
		vector<int> triangleTest_indexesInRWGs(triangles_test[r].indexesInRWGs);
		vector<double> triangleTest_signsInRWGs(triangles_test[r].signInRWG);

		// we now start the loop on the src triangles
		// loop on the source RWGs
		for (unsigned int s = 0; s < triangles_src.size(); ++s)
		{
			if ((target_surface != triangles_src[s].surface_num).all())
				continue;

			// the RWGs concerned by the source triangle
			vector<int> RWGsIndexes_src(triangles_src[s].RWGIndexes);
			vector<int> triangleSrc_indexesInRWGs(triangles_src[s].indexesInRWGs);
			vector<double> triangleSrc_signsInRWGs(triangles_src[s].signInRWG);

			double r_grav_obs_r_grav_src[3]
				= {triangles_test[r].r_grav[0] - triangles_src[s].r_grav[0],
					triangles_test[r].r_grav[1] - triangles_src[s].r_grav[1],
					triangles_test[r].r_grav[2] - triangles_src[s].r_grav[2]
				};
			double R_os = sqrt(dot3D(r_grav_obs_r_grav_src, r_grav_obs_r_grav_src));

			const bool IS_TOUCH
				= ((R_os - triangles_test[r].R_max - triangles_src[s].R_max) <= 0.0);
			const bool IS_NEAR
				= ((R_os - 1.5 * triangles_test[r].R_max - 1.5 * triangles_src[s].R_max) <= 0.0);
			const bool IS_SAME_TR
				= (triangles_test[r].number == triangles_src[s].number);

			IT_points_weights *IT_obs, *IT_src;
			if (IS_SAME_TR)
			{
				EXTRACT_1_R = (EXTRACT_R = 1);
				IT_obs = &IT_13;
				IT_src = &IT_13;
			}
			else if (IS_TOUCH)
			{
				EXTRACT_1_R = (EXTRACT_R = 1);
				IT_obs = &IT_13;
				IT_src = &IT_13;
			}
			else if (IS_NEAR)
			{
				EXTRACT_1_R = (EXTRACT_R = 1);
				IT_obs = &IT_6;
				IT_src = &IT_6;
			}
			else
			{
				EXTRACT_1_R = (EXTRACT_R = 0);
				IT_obs = &IT_6;
				IT_src = &IT_6;
			}

			// declaration of the scalars and vectors needed in the integrations
			complex<double> ITo_ITs_G, ITo_r_dot_ITs_G_rprime,
			                ITo_n_hat_X_r_dot_ITs_G_rprime, IDTo_l_hat_dot_r_ITs_G;
			complex<double> ITo_n_hat_X_r_dot_r_X_ITs_grad_G;
			complex<double> ITo_r_ITs_G[3],
			                ITo_ITs_G_rprime[3], IDTo_l_hat_ITs_G[3];
			complex<double> ITo_ITs_grad_G[3],
			                ITo_r_X_ITs_grad_G[3], ITo_n_hat_X_r_X_ITs_grad_G[3];

			ITo_ITs_free(ITo_ITs_G, ITo_r_ITs_G, ITo_ITs_G_rprime,
			             ITo_r_dot_ITs_G_rprime, ITo_n_hat_X_r_dot_ITs_G_rprime,
			             ITo_ITs_grad_G, ITo_r_X_ITs_grad_G,
			             ITo_n_hat_X_r_dot_r_X_ITs_grad_G,
			             ITo_n_hat_X_r_X_ITs_grad_G, triangles_test[r],
			             triangles_src[s], k, *IT_obs, *IT_src, EXTRACT_1_R,
			             EXTRACT_R);

			const complex<double>
				ITo_n_hat_X_r_ITs_G[3] = {n_hat[1] * ITo_r_ITs_G[2] - n_hat[2] * ITo_r_ITs_G[1],
					n_hat[2] * ITo_r_ITs_G[0] - n_hat[0] * ITo_r_ITs_G[2],
					n_hat[0] * ITo_r_ITs_G[1] - n_hat[1] * ITo_r_ITs_G[0]
				};
			const complex<double>
				n_hat_dot_ITo_r_X_ITs_grad_G(n_hat[0] * ITo_r_X_ITs_grad_G[0] + n_hat[1] * ITo_r_X_ITs_grad_G[1] + n_hat[2] * ITo_r_X_ITs_grad_G[2]);

			if ((IS_TOUCH) && (nEJ || nHM))
				IDTo_ITs_free(IDTo_l_hat_dot_r_ITs_G, IDTo_l_hat_ITs_G,
				              triangles_test[r], triangles_src[s], k,
				              IT_3, *IT_src, EXTRACT_1_R, EXTRACT_R);

			for (unsigned int p = 0; p < RWGsIndexes_test.size(); ++p)
			{
				const int index_p = RWGsIndexes_test[p];
				const int local_number_edge_p = index_p;
				const double l_p = test_RWGs[index_p].length;
				const double sign_edge_p = triangleTest_signsInRWGs[p];
				const double C_p = sign_edge_p * l_p * 0.5 / triangles_test[r].A;
				double *r_p, n_hat_X_r_p[3];

				if (triangleTest_indexesInRWGs[p] == 0)
					r_p = test_RWGs[index_p].vertexesCoord_0;
				else
					r_p = test_RWGs[index_p].vertexesCoord_3;
				cross3D(n_hat_X_r_p, n_hat, r_p);

				// temporary elements for Z_tE_J
				const complex<double>
					p_dot_ITo_ITs_G_rprime(r_p[0] * ITo_ITs_G_rprime[0] + r_p[1] * ITo_ITs_G_rprime[1] + r_p[2] * ITo_ITs_G_rprime[2]);

				// temporary elements for Z_nE_J
				const complex<double>
					n_hat_X_r_p_dot_ITo_ITs_G_rprime(n_hat_X_r_p[0] * ITo_ITs_G_rprime[0] + n_hat_X_r_p[1] * ITo_ITs_G_rprime[1] + n_hat_X_r_p[2] * ITo_ITs_G_rprime[2]);

				const complex<double>
					n_hat_X_r_p_dot_ITo_ITs_grad_G(n_hat_X_r_p[0] * ITo_ITs_grad_G[0] + n_hat_X_r_p[1] * ITo_ITs_grad_G[1] + n_hat_X_r_p[2] * ITo_ITs_grad_G[2]);

				// temporary elements for Z_tH_J
				const complex<double>
					r_p_X_ITo_ITs_grad_G[3] = {r_p[1] * ITo_ITs_grad_G[2] - r_p[2] * ITo_ITs_grad_G[1],
						r_p[2] * ITo_ITs_grad_G[0] - r_p[0] * ITo_ITs_grad_G[2],
						r_p[0] * ITo_ITs_grad_G[1] - r_p[1] * ITo_ITs_grad_G[0]
					};

				const complex<double>
					ITo_r_rp_X_ITs_grad_G[3] = {ITo_r_X_ITs_grad_G[0] - r_p_X_ITo_ITs_grad_G[0],
						ITo_r_X_ITs_grad_G[1] - r_p_X_ITo_ITs_grad_G[1],
						ITo_r_X_ITs_grad_G[2] - r_p_X_ITo_ITs_grad_G[2]
					};

				const double
					n_hat_X_rp_dot_IT_r(n_hat_X_r_p[0] * IT_r[0] + n_hat_X_r_p[1] * IT_r[1] + n_hat_X_r_p[2] * IT_r[2]);

				// temporary elements for Z_nH_J
				const complex<double>
					n_hat_X_r_p_X_ITo_ITs_grad_G[3] = {n_hat_X_r_p[1] * ITo_ITs_grad_G[2]
						- n_hat_X_r_p[2] * ITo_ITs_grad_G[1],
						n_hat_X_r_p[2] * ITo_ITs_grad_G[0]
						- n_hat_X_r_p[0] * ITo_ITs_grad_G[2],
						n_hat_X_r_p[0] * ITo_ITs_grad_G[1]
						- n_hat_X_r_p[1] * ITo_ITs_grad_G[0]
					};

				const complex<double>
					ITo_n_hat_X_r_rp_X_ITs_grad_G[3] = {ITo_n_hat_X_r_X_ITs_grad_G[0] - n_hat_X_r_p_X_ITo_ITs_grad_G[0],
						ITo_n_hat_X_r_X_ITs_grad_G[1] - n_hat_X_r_p_X_ITo_ITs_grad_G[1],
						ITo_n_hat_X_r_X_ITs_grad_G[2] - n_hat_X_r_p_X_ITo_ITs_grad_G[2]
					};

				const complex<double>
					n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G(n_hat_X_r_p[0] * ITo_r_X_ITs_grad_G[0] + n_hat_X_r_p[1] * ITo_r_X_ITs_grad_G[1] + n_hat_X_r_p[2] * ITo_r_X_ITs_grad_G[2]);

				for (unsigned int q = 0; q < RWGsIndexes_src.size(); ++q)
				{
					const int index_q = RWGsIndexes_src[q];
					const int local_number_edge_q = index_q;
					const double l_q = src_RWGs[index_q].length;
					const double sign_edge_q = triangleSrc_signsInRWGs[q];
					const double C_pq = C_p * sign_edge_q * l_q * 0.5 / triangles_src[s].A;
					double* r_q;
					if (triangleSrc_indexesInRWGs[q] == 0)
						r_q = src_RWGs[index_q].vertexesCoord_0;
					else
						r_q = src_RWGs[index_q].vertexesCoord_3;

					const double
						rp_dot_rq(r_p[0] * r_q[0] + r_p[1] * r_q[1] + r_p[2] * r_q[2]);
					const double rp_plus_rq_dot_ITr((r_p[0] + r_q[0]) * IT_r[0]
						+ (r_p[1] + r_q[1]) * IT_r[1]
						+ (r_p[2] + r_q[2]) * IT_r[2]);
					complex<double> z_pq, D_mn_1, D_mn_2;

					// <f_p ; EFIE> : Z_tE_J computation. Z_tH_M = eps/mu * Z_tE_J
					if (tEJ || tHM)
					{
						D_mn_1 = (-4.0 * C_pq) * ITo_ITs_G;
						D_mn_2 = C_pq * (ITo_r_dot_ITs_G_rprime
							- (ITo_r_ITs_G[0] * r_q[0]
								+ ITo_r_ITs_G[1] * r_q[1]
								+ ITo_r_ITs_G[2] * r_q[2])
							- p_dot_ITo_ITs_G_rprime + rp_dot_rq * ITo_ITs_G);
						z_pq = z_pq_factor * (D_mn_1 + k_square * D_mn_2);

						if (tEJ)
							Z_CFIE_EJ(local_number_edge_p, local_number_edge_q)
								+= (signSurfObs * signSurfSrc) * CFIE[0] * z_pq;
						if (tHM)
							Z_CFIE_HM(local_number_edge_p, local_number_edge_q)
								+= (signSurfObs * signSurfSrc) * CFIE[2] * eps / mu * z_pq;
					}

					// <n x f_p ; EFIE> : Z_nE_J computation. Z_nH_M = eps/mu * Z_nE_J
					if (nEJ || nHM)
					{
						if (IS_TOUCH)
							D_mn_1 = 2.0 * C_pq
								* (-IDTo_l_hat_dot_r_ITs_G
									+ (r_p[0] * IDTo_l_hat_ITs_G[0]
										+ r_p[1] * IDTo_l_hat_ITs_G[1]
										+ r_p[2] * IDTo_l_hat_ITs_G[2]));
						else
							D_mn_1 = 2.0 * C_pq * (n_hat_dot_ITo_r_X_ITs_grad_G
								- n_hat_X_r_p_dot_ITo_ITs_grad_G);

						D_mn_2 = C_pq * (ITo_n_hat_X_r_dot_ITs_G_rprime
							- (ITo_n_hat_X_r_ITs_G[0] * r_q[0]
								+ ITo_n_hat_X_r_ITs_G[1] * r_q[1]
								+ ITo_n_hat_X_r_ITs_G[2] * r_q[2])
							- n_hat_X_r_p_dot_ITo_ITs_G_rprime
							+ (n_hat_X_r_p[0] * r_q[0]
								+ n_hat_X_r_p[1] * r_q[1]
								+ n_hat_X_r_p[2] * r_q[2]) * ITo_ITs_G);
						// z_pq_factor = 1.0/(I*w*eps)
						z_pq = z_pq_factor * (D_mn_1 + k_square * D_mn_2);

						if (nEJ)
							Z_CFIE_EJ(local_number_edge_p, local_number_edge_q)
								+= (signSurfObs * signSurfSrc) * CFIE[1] * z_pq;

						if (nHM)
							Z_CFIE_HM(local_number_edge_p, local_number_edge_q)
								+= (signSurfObs * signSurfSrc) * CFIE[3] * eps / mu * z_pq;
					}

					// <f_p ; MFIE> : Z_tH_J computation. Z_tE_M = -Z_tH_J
					if (tHJ || tEM)
					{
						if (!IS_SAME_TR)
						{
							z_pq = (signSurfObs * signSurfSrc * C_pq)
								* ((r_p[0] - r_q[0]) * ITo_r_rp_X_ITs_grad_G[0]
									+ (r_p[1] - r_q[1]) * ITo_r_rp_X_ITs_grad_G[1]
									+ (r_p[2] - r_q[2]) * ITo_r_rp_X_ITs_grad_G[2]);
						}
						// else we have: z_pq = 0.5 * <f_p ; n x f_q>
						else
						{
							z_pq = (signSurfObs * 0.5 * C_pq * solid_angle_factor)
								* (n_hat_X_rp_dot_IT_r
									+ (r_q[0] * IT_n_hat_X_r[0]
										+ r_q[1] * IT_n_hat_X_r[1]
										+ r_q[2] * IT_n_hat_X_r[2])
									- (r_q[0] * n_hat_X_r_p[0]
										+ r_q[1] * n_hat_X_r_p[1] + r_q[2] * n_hat_X_r_p[2])
									* triangles_test[r].A);
						}

						if (tHJ)
							Z_CFIE_HJ(local_number_edge_p, local_number_edge_q) += CFIE[2] * z_pq;

						if (tEM)
							Z_CFIE_EM(local_number_edge_p, local_number_edge_q) -= CFIE[0] * z_pq;
					}

					// <n x f_p ; MFIE> : Z_nH_J computation. Z_nE_M = -Z_nH_J
					if (nHJ || nEM)
					{
						if (!IS_SAME_TR)
							z_pq = (-signSurfObs * signSurfSrc * C_pq)
								* (ITo_n_hat_X_r_dot_r_X_ITs_grad_G
									+ (r_q[0] * ITo_n_hat_X_r_rp_X_ITs_grad_G[0]
										+ r_q[1] * ITo_n_hat_X_r_rp_X_ITs_grad_G[1]
										+ r_q[2] * ITo_n_hat_X_r_rp_X_ITs_grad_G[2])
									- n_hat_X_r_p_dot_ITo_r_X_ITs_grad_G);
						// else we have: z_pq = 0.5 * <n x f_p ; n x f_q> = 0.5 * <f_p ; f_q>
						else
						{
							z_pq = (signSurfObs * 0.5 * C_pq * solid_angle_factor)
								* (IT_r_square - rp_plus_rq_dot_ITr
									+ rp_dot_rq * triangles_test[r].A);
						}

						if (nHJ)
							Z_CFIE_HJ(local_number_edge_p, local_number_edge_q)
								+= CFIE[3] * z_pq;

						if (nEM)
							Z_CFIE_EM(local_number_edge_p, local_number_edge_q)
								+= CFIE[1] * z_pq;
					}
				} // for q
			} // for p
		} // for s
	} // for r
}

void wrap_Z_CFIE_computation(complex<double> Z_CFIE_EJ_tmp[],
                             complex<double> Z_CFIE_EM_tmp[],
                             complex<double> Z_CFIE_HJ_tmp[],
                             complex<double> Z_CFIE_HM_tmp[],
                             complex<double> CFIE[],
                             int CFIE_bool_tmp[],
                             int M_current_tmp,
                             double signSurfObs,
                             double signSurfSrc,
                             int N_RWG, int N_RWG_src, int N_RWG_test, int V,
                             int numbers_RWG_src_tmp[],
                             int numbers_RWG_test_tmp[],
                             int srcRWGNumber_surface[],
                             int RWGNumber_signedTriangles_tmp[],
                             int RWGNumber_nodes_tmp[],
                             double nodesCoord_tmp[],
                             double w,
                             complex<double> eps_r,
                             complex<double> mu_r,
                             int N_target_surface, int target_surface_tmp[])
{
	typedef Array<int, Dynamic, Dynamic, RowMajor> ArrayXXi_RowMajor;
	typedef Array<double, Dynamic, Dynamic, RowMajor> ArrayXXd_RowMajor;

	ArrayXXcd Z_CFIE_EJ(N_RWG_test, N_RWG_src), Z_CFIE_EM(N_RWG_test, N_RWG_src);
	ArrayXXcd Z_CFIE_HJ(N_RWG_test, N_RWG_src), Z_CFIE_HM(N_RWG_test, N_RWG_src);

	ArrayXi numbers_RWG_src = Map<ArrayXi>(numbers_RWG_src_tmp, N_RWG_src);
	ArrayXi numbers_RWG_test = Map<ArrayXi>(numbers_RWG_test_tmp, N_RWG_test);

	ArrayXi target_surface = Map<ArrayXi>(target_surface_tmp, N_target_surface);
	ArrayXXi RWGNumber_signedTriangles = Map<ArrayXXi_RowMajor>(RWGNumber_signedTriangles_tmp, N_RWG, 2);
	ArrayXXi RWGNumber_nodes = Map<ArrayXXi_RowMajor>(RWGNumber_nodes_tmp, N_RWG, 4);
	ArrayXXd nodesCoord = Map<ArrayXXd_RowMajor>(nodesCoord_tmp, V, 3);

	bool CFIE_bool[4];
	for (int i = 0; i < 4; i++)
	{
		if (CFIE_bool_tmp[i] == 0) CFIE_bool[i] = false;
		else CFIE_bool[i] = true;
	}
	bool M_current;
	if (M_current_tmp == 0) M_current = false;
	else M_current = true;

	Z_CFIE_computation(Z_CFIE_EJ, Z_CFIE_EM, Z_CFIE_HJ, Z_CFIE_HM,
	                   CFIE, CFIE_bool, M_current, signSurfObs, signSurfSrc,
	                   numbers_RWG_src, numbers_RWG_test, srcRWGNumber_surface,
	                   RWGNumber_signedTriangles, RWGNumber_nodes, nodesCoord,
	                   w, eps_r, mu_r, target_surface);

	int index = 0;
	for (int i = 0; i < N_RWG_test; i++)
	{
		for (int j = 0; j < N_RWG_src; j++)
		{
			Z_CFIE_EJ_tmp[index] = Z_CFIE_EJ(i, j);
			Z_CFIE_EM_tmp[index] = Z_CFIE_EM(i, j);
			Z_CFIE_HJ_tmp[index] = Z_CFIE_HJ(i, j);
			Z_CFIE_HM_tmp[index] = Z_CFIE_HM(i, j);
			index++;
		}
	}
}

void Z_CFIE_test_computation(ArrayXXcd& Z_CFIE,
                             const int N_RWG_src, const int N_RWG_test,
                             const ArrayXi& numbers_RWG_src,
                             const int srcRWGNumber_surface[],
                             const ArrayXXi& RWGNumber_signedTriangles,
                             const ArrayXXi& RWGNumber_nodes,
                             const ArrayXXd& nodesCoord,
                             const double w,
                             const complex<double>& eps_r,
                             const complex<double>& mu_r,
                             const ArrayXi& target_surface)
{
	// half RWGs construction
	vector<RWG> src_RWGs, test_RWGs;
	src_RWGs.reserve(N_RWG_src);
	test_RWGs.reserve(N_RWG_test);

	for (int i = 0; i < N_RWG_src; ++i)
	{
		int triangle_numbers[2], triangle_signs[2];
		triangle_numbers[0] = abs(RWGNumber_signedTriangles(i, 0));
		triangle_numbers[1] = abs(RWGNumber_signedTriangles(i, 1));
		triangle_signs[0] = 1;
		triangle_signs[1] = -1;

		// the nodes of the RWG
		const int n0 = RWGNumber_nodes(i, 0);
		const int n1 = RWGNumber_nodes(i, 1);
		const int n2 = RWGNumber_nodes(i, 2);
		const int n3 = RWGNumber_nodes(i, 3);
		const double r0[3]
			= {nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
		const double r1[3]
			= {nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
		const double r2[3]
			= {nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
		const double r3[3]
			= {nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};

		src_RWGs.push_back(RWG(numbers_RWG_src(i), srcRWGNumber_surface[i],
		                       triangle_numbers, triangle_signs, r0, r1, r2, r3));
	}

	for (int i = 0; i < N_RWG_test; ++i)
	{
		int triangle_numbers[2], triangle_signs[2];

		triangle_numbers[0] = abs(RWGNumber_signedTriangles(i, 0));
		triangle_numbers[1] = abs(RWGNumber_signedTriangles(i, 1));
		triangle_signs[0] = 1;
		triangle_signs[1] = -1;

		// the nodes of the RWG
		const int n0 = RWGNumber_nodes(i, 0);
		const int n1 = RWGNumber_nodes(i, 1);
		const int n2 = RWGNumber_nodes(i, 2);
		const int n3 = RWGNumber_nodes(i, 3);
		const double r0[3] = {nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
		const double r1[3] = {nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
		const double r2[3] = {nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
		const double r3[3] = {nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};

		test_RWGs.push_back(RWG(numbers_RWG_src(i), srcRWGNumber_surface[i],
		                        triangle_numbers,
		                        triangle_signs, r0, r1, r2, r3));
	}

	// triangles
	vector< Dictionary<int, int> > srcTriangleToRWG, testTriangleToRWG;
	srcTriangleToRWG.reserve(N_RWG_src * 2);
	testTriangleToRWG.reserve(N_RWG_test * 2);
	for (unsigned int i = 0; i < src_RWGs.size(); ++i)
	{
		srcTriangleToRWG.push_back(Dictionary<int, int>(src_RWGs[i].triangleNumbers[0], i));
		srcTriangleToRWG.push_back(Dictionary<int, int>(src_RWGs[i].triangleNumbers[1], i));
	}

	for (unsigned int i = 0; i < test_RWGs.size(); ++i)
	{
		testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers[0], i));
		testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers[1], i));
	}

	sort(srcTriangleToRWG.begin(), srcTriangleToRWG.end());
	sort(testTriangleToRWG.begin(), testTriangleToRWG.end());
	vector<Triangle> triangles_src, triangles_test;
	constructVectorTriangles(triangles_src, src_RWGs, srcTriangleToRWG);
	constructVectorTriangles(triangles_test, test_RWGs, testTriangleToRWG);

	Z_CFIE = ArrayXXcd::Zero(N_RWG_test, N_RWG_src);

	// loop on the observation RWGs
	for (unsigned int r = 0; r < triangles_test.size(); ++r)
	{
		if ((target_surface != triangles_test[r].surface_num).all())
			continue;

		// computation of the triangle-to-triangle terms
		double* n_hat;
		n_hat = triangles_test[r].n_hat;
		double IT_r_square; // n_hat_X_r_p_dot_IT_r;
		double IT_r[3], IT_n_hat_X_r[3];

		// serves for <f_m ; f_n> and <f_m ; n x f_n>
		IT_fm_fn(IT_r_square, IT_r, triangles_test[r]);
		// serves for <f_m ; n x f_n>
		cross3D(IT_n_hat_X_r, n_hat, IT_r);

		// the RWGs concerned by the test triangle
		vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
		vector<int> triangleTest_indexesInRWGs(triangles_test[r].indexesInRWGs);
		vector<double> triangleTest_signsInRWGs(triangles_test[r].signInRWG);

		// we now start the loop on the src triangles
		// loop on the source RWGs
		unsigned int s = r;

		// the RWGs concerned by the source triangle
		vector<int> RWGsIndexes_src(triangles_src[s].RWGIndexes);
		vector<int> triangleSrc_indexesInRWGs(triangles_src[s].indexesInRWGs);
		vector<double> triangleSrc_signsInRWGs(triangles_src[s].signInRWG);

		for (unsigned int p = 0; p < RWGsIndexes_test.size(); ++p)
		{
			const int index_p = RWGsIndexes_test[p];
			const int local_number_edge_p = index_p;
			const double l_p = test_RWGs[index_p].length;
			const double sign_edge_p = triangleTest_signsInRWGs[p];
			const double C_p = sign_edge_p * l_p * 0.5 / triangles_test[r].A;
			double* r_p;

			if (triangleTest_indexesInRWGs[p] == 0)
				r_p = test_RWGs[index_p].vertexesCoord_0;
			else
				r_p = test_RWGs[index_p].vertexesCoord_3;

			for (unsigned int q = 0; q < RWGsIndexes_src.size(); ++q)
			{
				const int index_q = RWGsIndexes_src[q];
				const int local_number_edge_q = index_q;
				const double l_q = src_RWGs[index_q].length;
				const double sign_edge_q = triangleSrc_signsInRWGs[q];
				const double C_pq = C_p * sign_edge_q * l_q * 0.5 / triangles_src[s].A;
				double* r_q;
				if (triangleSrc_indexesInRWGs[q] == 0)
					r_q = src_RWGs[index_q].vertexesCoord_0;
				else
					r_q = src_RWGs[index_q].vertexesCoord_3;

				const double
					rp_dot_rq(r_p[0] * r_q[0] + r_p[1] * r_q[1] + r_p[2] * r_q[2]);
				const double rp_plus_rq_dot_ITr((r_p[0] + r_q[0]) * IT_r[0]
					+ (r_p[1] + r_q[1]) * IT_r[1] + (r_p[2] + r_q[2]) * IT_r[2]);
				complex<double> z_pq;

				z_pq = C_pq * (IT_r_square - rp_plus_rq_dot_ITr
					+ rp_dot_rq * triangles_test[r].A);

				Z_CFIE(local_number_edge_p, local_number_edge_q) += z_pq;
			}
		}
	}
}

void wrap_Z_CFIE_test_computation(complex<double> Z_CFIE_tmp[],
                                  int N_RWG, int V,
                                  int srcRWGNumber_surface[],
                                  int RWGNumber_signedTriangles_tmp[],
                                  int RWGNumber_nodes_tmp[],
                                  double nodesCoord_tmp[],
                                  double w,
                                  complex<double> eps_r,
                                  complex<double> mu_r,
                                  int N_target_surface, int target_surface_tmp[])
{
	typedef Array<int, Dynamic, Dynamic, RowMajor> ArrayXXi_RowMajor;
	typedef Array<double, Dynamic, Dynamic, RowMajor> ArrayXXd_RowMajor;

	ArrayXXcd Z_CFIE(N_RWG, N_RWG);
	ArrayXi numbers_RWG_src = ArrayXi::LinSpaced(N_RWG, 0, N_RWG - 1);
	ArrayXi target_surface = Map<ArrayXi>(target_surface_tmp, N_target_surface);
	ArrayXXi RWGNumber_signedTriangles = Map<ArrayXXi_RowMajor>(RWGNumber_signedTriangles_tmp, N_RWG, 2);
	ArrayXXi RWGNumber_nodes = Map<ArrayXXi_RowMajor>(RWGNumber_nodes_tmp, N_RWG, 4);
	ArrayXXd nodesCoord = Map<ArrayXXd_RowMajor>(nodesCoord_tmp, V, 3);

	Z_CFIE_test_computation(Z_CFIE, N_RWG, N_RWG, numbers_RWG_src, srcRWGNumber_surface,
	                        RWGNumber_signedTriangles, RWGNumber_nodes, nodesCoord,
	                        w, eps_r, mu_r, target_surface);

	int index = 0;
	for (int i = 0; i < N_RWG; i++)
	{
		for (int j = 0; j < N_RWG; j++)
		{
			Z_CFIE_tmp[index] = Z_CFIE(i, j);
			index++;
		}
	}
}


// void Z_CFIE_test_computation(ArrayXXcd& Z_CFIE,
//                              const int N_RWG_src, const int N_RWG_test,
//                              const ArrayXi& numbers_RWG_src,
//                              const int srcRWGNumber_surface[],
//                              const ArrayXXi& RWGNumber_signedTriangles,
//                              const ArrayXXi& RWGNumber_nodes,
//                              const ArrayXXd& nodesCoord,
//                              const double w,
//                              const complex<double>& eps_r,
//                              const complex<double>& mu_r,
//                              const ArrayXi& target_surface)
// {
// 	// half RWGs construction
// 	vector<RWG> src_RWGs, test_RWGs;
// 	src_RWGs.reserve(N_RWG_src);
// 	test_RWGs.reserve(N_RWG_test);

// 	for (int i = 0; i < N_RWG_src; ++i)
// 	{
// 		int triangle_numbers[2], triangle_signs[2];
// 		triangle_numbers[0] = abs(RWGNumber_signedTriangles(i, 0));
// 		triangle_numbers[1] = abs(RWGNumber_signedTriangles(i, 1));
// 		triangle_signs[0] = 1;
// 		triangle_signs[1] = -1;

// 		// the nodes of the RWG
// 		const int n0 = RWGNumber_nodes(i, 0);
// 		const int n1 = RWGNumber_nodes(i, 1);
// 		const int n2 = RWGNumber_nodes(i, 2);
// 		const int n3 = RWGNumber_nodes(i, 3);
// 		const double r0[3]
// 			= {nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
// 		const double r1[3]
// 			= {nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
// 		const double r2[3]
// 			= {nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
// 		const double r3[3]
// 			= {nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};

// 		src_RWGs.push_back(RWG(numbers_RWG_src(i), srcRWGNumber_surface[i],
// 		                       triangle_numbers, triangle_signs, r0, r1, r2, r3));
// 	}

// 	for (int i = 0; i < N_RWG_test; ++i)
// 	{
// 		int triangle_numbers[2], triangle_signs[2];

// 		triangle_numbers[0] = abs(RWGNumber_signedTriangles(i, 0));
// 		triangle_numbers[1] = abs(RWGNumber_signedTriangles(i, 1));
// 		triangle_signs[0] = 1;
// 		triangle_signs[1] = -1;

// 		// the nodes of the RWG
// 		const int n0 = RWGNumber_nodes(i, 0);
// 		const int n1 = RWGNumber_nodes(i, 1);
// 		const int n2 = RWGNumber_nodes(i, 2);
// 		const int n3 = RWGNumber_nodes(i, 3);
// 		const double r0[3] = {nodesCoord(n0, 0), nodesCoord(n0, 1), nodesCoord(n0, 2)};
// 		const double r1[3] = {nodesCoord(n1, 0), nodesCoord(n1, 1), nodesCoord(n1, 2)};
// 		const double r2[3] = {nodesCoord(n2, 0), nodesCoord(n2, 1), nodesCoord(n2, 2)};
// 		const double r3[3] = {nodesCoord(n3, 0), nodesCoord(n3, 1), nodesCoord(n3, 2)};

// 		test_RWGs.push_back(RWG(numbers_RWG_src(i), srcRWGNumber_surface[i],
// 		                        triangle_numbers,
// 		                        triangle_signs, r0, r1, r2, r3));
// 	}

// 	// triangles
// 	vector< Dictionary<int, int> > srcTriangleToRWG, testTriangleToRWG;
// 	srcTriangleToRWG.reserve(N_RWG_src * 2);
// 	testTriangleToRWG.reserve(N_RWG_test * 2);
// 	for (unsigned int i = 0; i < src_RWGs.size(); ++i)
// 	{
// 		srcTriangleToRWG.push_back(Dictionary<int, int>(src_RWGs[i].triangleNumbers[0], i));
// 		srcTriangleToRWG.push_back(Dictionary<int, int>(src_RWGs[i].triangleNumbers[1], i));
// 	}

// 	for (unsigned int i = 0; i < test_RWGs.size(); ++i)
// 	{
// 		testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers[0], i));
// 		testTriangleToRWG.push_back(Dictionary<int, int>(test_RWGs[i].triangleNumbers[1], i));
// 	}

// 	sort(srcTriangleToRWG.begin(), srcTriangleToRWG.end());
// 	sort(testTriangleToRWG.begin(), testTriangleToRWG.end());
// 	vector<Triangle> triangles_src, triangles_test;
// 	constructVectorTriangles(triangles_src, src_RWGs, srcTriangleToRWG);
// 	constructVectorTriangles(triangles_test, test_RWGs, testTriangleToRWG);

// 	Z_CFIE = ArrayXXcd::Zero(N_RWG_test, N_RWG_src);

// 	// loop on the observation RWGs
// 	for (unsigned int r = 0; r < triangles_test.size(); ++r)
// 	{
// 		if ((target_surface != triangles_test[r].surface_num).all())
// 			continue;

// 		// computation of the triangle-to-triangle terms
// 		double* n_hat;
// 		n_hat = triangles_test[r].n_hat;
// 		double IT_r_square; // n_hat_X_r_p_dot_IT_r;
// 		double IT_r[3], IT_n_hat_X_r[3];

// 		// serves for <f_m ; f_n> and <f_m ; n x f_n>
// 		IT_fm_fn(IT_r_square, IT_r, triangles_test[r]);
// 		// serves for <f_m ; n x f_n>
// 		cross3D(IT_n_hat_X_r, n_hat, IT_r);

// 		// the RWGs concerned by the test triangle
// 		vector<int> RWGsIndexes_test(triangles_test[r].RWGIndexes);
// 		vector<int> triangleTest_indexesInRWGs(triangles_test[r].indexesInRWGs);
// 		vector<double> triangleTest_signsInRWGs(triangles_test[r].signInRWG);

// 		// we now start the loop on the src triangles
// 		// loop on the source RWGs
// 		unsigned int s = r;

// 		// the RWGs concerned by the source triangle
// 		vector<int> RWGsIndexes_src(triangles_src[s].RWGIndexes);
// 		vector<int> triangleSrc_indexesInRWGs(triangles_src[s].indexesInRWGs);
// 		vector<double> triangleSrc_signsInRWGs(triangles_src[s].signInRWG);

// 		for (unsigned int p = 0; p < RWGsIndexes_test.size(); ++p)
// 		{
// 			const int index_p = RWGsIndexes_test[p];
// 			const int local_number_edge_p = index_p;
// 			const double l_p = test_RWGs[index_p].length;
// 			const double sign_edge_p = triangleTest_signsInRWGs[p];
// 			const double C_p = sign_edge_p * l_p * 0.5 / triangles_test[r].A;
// 			double* r_p;

// 			if (triangleTest_indexesInRWGs[p] == 0)
// 				r_p = test_RWGs[index_p].vertexesCoord_0;
// 			else
// 				r_p = test_RWGs[index_p].vertexesCoord_3;

// 			for (unsigned int q = 0; q < RWGsIndexes_src.size(); ++q)
// 			{
// 				const int index_q = RWGsIndexes_src[q];
// 				const int local_number_edge_q = index_q;
// 				const double l_q = src_RWGs[index_q].length;
// 				const double sign_edge_q = triangleSrc_signsInRWGs[q];
// 				const double C_pq = C_p * sign_edge_q * l_q * 0.5 / triangles_src[s].A;
// 				double* r_q;
// 				if (triangleSrc_indexesInRWGs[q] == 0)
// 					r_q = src_RWGs[index_q].vertexesCoord_0;
// 				else
// 					r_q = src_RWGs[index_q].vertexesCoord_3;

// 				const double
// 					rp_dot_rq(r_p[0] * r_q[0] + r_p[1] * r_q[1] + r_p[2] * r_q[2]);
// 				const double rp_plus_rq_dot_ITr((r_p[0] + r_q[0]) * IT_r[0]
// 					+ (r_p[1] + r_q[1]) * IT_r[1] + (r_p[2] + r_q[2]) * IT_r[2]);
// 				complex<double> z_pq;

// 				z_pq = C_pq * (IT_r_square - rp_plus_rq_dot_ITr
// 					+ rp_dot_rq * triangles_test[r].A);

// 				Z_CFIE(local_number_edge_p, local_number_edge_q) += z_pq;
// 			}
// 		}
// 	}
// }