/**********************************************************************
 *
 * alpha_computation.h
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

#ifndef ALPHA_COMPUTATION_H
#define ALPHA_COMPUTATION_H

#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <boost/math/special_functions/hankel.hpp>
#include <Eigen/Dense>
#include <math.h>

using namespace std;
using namespace Eigen;

#include "EMConstants.hpp"
//#include "integr_1D_X_W.hpp"
//#include "GK_triangle.hpp"
//#include "./amos/zbesh/zbesh_interface.h"

template <typename T>
void P_Legendre(Array<T, Dynamic, 1>& P,
                const T z)
{
	int L = P.size() - 1;
	P(0) = 1;
	if (L >= 1)
	{
		P(1) = z;
		for (int i = 1; i < L; i++) P(i + 1) = ((2.0 * i + 1.0) * z * P(i) - i * P(i - 1)) / (i + 1.0);
	}
}

std::complex<double> alpha_computation
(const double& theta,
 const double& phi,
 const ArrayXcd& h2_sph,
 const double r_mn[],
 const int L,
 const int L_prime,
 const std::complex<double>& k)
{
	const double norm_r_mn = sqrt(r_mn[0] * r_mn[0] + r_mn[1] * r_mn[1] + r_mn[2] * r_mn[2]);
	const double sin_theta = sin(theta), cos_theta = cos(theta);
	const double r_mn_hat[3] = {r_mn[0] / norm_r_mn, r_mn[1] / norm_r_mn, r_mn[2] / norm_r_mn};
	const double s_hat[3] = {sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
	ArrayXd P_Leg(L_prime + 1);
	ArrayXcd coeff(L_prime + 1);
	P_Legendre(P_Leg, (s_hat[0] * r_mn_hat[0] + s_hat[1] * r_mn_hat[1] + s_hat[2] * r_mn_hat[2]));
//	cout << theta << " " << phi << " " << P_Leg(L_prime) << endl;
	for (int j = 0; j < L + 1; ++j) coeff(j) = pow(-I, j) * (2 * j + 1.0) * h2_sph(j);
	const std::complex<double> const_coeff_lissage = pow(-I, L) * (2 * L + 1.0) * h2_sph(L);
	for (int j = L + 1; j < L_prime + 1; ++j) coeff(j) = const_coeff_lissage * pow(cos((j - L) * M_PI / 2.0 / (L_prime - L)), 2);
	return (-I * k / (16.0 * M_PI * M_PI) * (coeff * P_Leg).sum());
}

template <typename T>
void IT_theta_IT_phi_alpha_C
(Array<std::complex<T>, Dynamic, Dynamic>& alpha, /**< OUTPUT: 2D array \f$ \alpha \f$ */
 const double r_mn[], /**< INPUT: \f$ r_{mn} = r_m - r_n \f$ */
 const std::complex<double>& k, /**< INPUT: the wavenumber */
 const int L, /**< the expansion number */
 const int L_prime, /**< the second expansion number */
 const Array<T, Dynamic, 1>& Xtheta, /**< 1D array of \f$ \theta \f$ angles */
 const Array<T, Dynamic, 1>& Xphi) /**< 1D array of \f$ \phi \f$ angles */
{
	int kode = 1, M = 2, N = 1, nz, ierr;
	const double norm_r_mn
		= sqrt(r_mn[0] * r_mn[0] + r_mn[1] * r_mn[1] + r_mn[2] * r_mn[2]);
	ArrayXcd h2_sph(L + 1);
	ArrayXd j_sph(L + 1), y_sph(L + 1);
	std::complex<double> z(k * norm_r_mn);
	//	for (int i = 0; i < L + 1; i++)
	//	{
	//		double fnu = i + 0.5;
	//		zbesh(z, fnu, kode, M, N, h2_sph(i), nz, ierr);
	//	}
	//	h2_sph *= sqrt(M_PI / (2.0 * z));
	//	gsl_sf_bessel_jl_array(L, real(z), j_sph.data());
	//	gsl_sf_bessel_yl_array(L, real(z), y_sph.data());
	//	h2_sph = j_sph - I * y_sph;
	for (int i = 0; i < L + 1; i++) h2_sph(i) = boost::math::sph_hankel_2(i, real(z));

	for (int index_theta = 0; index_theta < Xtheta.size(); index_theta++)
	{
		for (int index_phi = 0; index_phi < Xphi.size(); index_phi++)
		{
			double theta = static_cast<double>(Xtheta(index_theta));
			double phi = static_cast<double>(Xphi(index_phi));
			alpha(index_theta, index_phi)
				= static_cast<complex<T>>
				(alpha_computation(theta, phi, h2_sph, r_mn, L, L_prime, k));
		}
	}
}

template <typename T>
void IT_theta_IT_phi_alpha_C2
(Array<complex<T>, Dynamic, 1>& alpha, /**< OUTPUT: 2D array \f$ \alpha \f$ */
 const double r_mn[], /**< INPUT: \f$ r_{mn} = r_m - r_n \f$ */
 const complex<double>& k, /**< INPUT: the wavenumber */
 const int L, /**< the expansion number */
 const int L_prime, /**< the second expansion number */
 const Array<T, Dynamic, Dynamic>& thetasPhis)
{
	int kode = 1, M = 2, N = 1, nz, ierr;
	const double norm_r_mn
		= sqrt(r_mn[0] * r_mn[0] + r_mn[1] * r_mn[1] + r_mn[2] * r_mn[2]);
	ArrayXcd h2_sph(L + 1);
	ArrayXd j_sph(L + 1), y_sph(L + 1);
	complex<double> z(k * norm_r_mn);
//	for (int i = 0; i < L + 1; i++)
//	{
//		double fnu = i + 0.5;
//		zbesh(z, fnu, kode, M, N, h2_sph(i), nz, ierr);
//	}
	// h2_sph *= sqrt(M_PI/(2.0 * z));
	//	gsl_sf_bessel_jl_array(L, real(z), j_sph.data());
	//	gsl_sf_bessel_yl_array(L, real(z), y_sph.data());
	//	h2_sph = j_sph - I * y_sph;
	for (int i = 0; i < L + 1; i++) h2_sph(i) = boost::math::sph_hankel_2(i, real(z));
//	cout << h2_sph(L) << endl;
//	cout << real(z) << k << norm_r_mn << endl;

	const int N_directions = alpha.size();
	for (int i = 0; i < N_directions; ++i)
	{
		double theta = static_cast<double>(thetasPhis(i, 0));
		double phi = static_cast<double>(thetasPhis(i, 1));
		alpha(i) = static_cast<std::complex<T>>(alpha_computation(theta, phi, h2_sph, r_mn, L, L_prime, k));
	}
}

template <typename T>
void IT_theta_IT_phi_alpha(Array<std::complex<T>, Dynamic, Dynamic>& alpha,
                           const Array<std::complex<double>, Dynamic, 1>& h2_sph,
                           const double r_mn[],
                           const std::complex<double>& k,
                           const ArrayXd& Xtheta,
                           const ArrayXd& Xphi)
{
	const double norm_r_mn = sqrt(r_mn[0] * r_mn[0] + r_mn[1] * r_mn[1] + r_mn[2] * r_mn[2]);
	const int L = h2_sph.size() - 1, L_prime = h2_sph.size() - 1;
	for (int index_theta = 0; index_theta < Xtheta.size(); index_theta++)
	{
		for (int index_phi = 0; index_phi < Xphi.size(); index_phi++)
		{
			alpha(index_theta, index_phi) = static_cast<std::complex<T>>(alpha_computation(Xtheta(index_theta), Xphi(index_phi), h2_sph, r_mn, L, L_prime, k));
		}
	}
}

#endif
