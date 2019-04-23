#include <Eigen/Dense>
using namespace Eigen;

#include "make_trans.hpp"
#include "interpolation.hpp"

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

std::complex<double> alpha_computation(const double& theta,
                                       const double& phi,
                                       const ArrayXcd& h2_sph,
                                       const double r_mn[],
                                       const int L,
                                       const int L_prime,
                                       const std::complex<double>& k)
{
	const complex<double> I(0.0, 1.0);
	const double norm_r_mn = sqrt(r_mn[0] * r_mn[0] + r_mn[1] * r_mn[1] + r_mn[2] * r_mn[2]);
	const double sin_theta = sin(theta), cos_theta = cos(theta);
	const double r_mn_hat[3] = {r_mn[0] / norm_r_mn, r_mn[1] / norm_r_mn, r_mn[2] / norm_r_mn};
	const double s_hat[3] = {sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
	ArrayXd P_Leg(L_prime + 1);
	ArrayXcd coeff(L_prime + 1);
	P_Legendre(P_Leg, (s_hat[0] * r_mn_hat[0] + s_hat[1] * r_mn_hat[1] + s_hat[2] * r_mn_hat[2]));
	for (int j = 0; j < L + 1; ++j) coeff(j) = pow(-I, j) * (2 * j + 1.0) * h2_sph(j);
	const std::complex<double> const_coeff_lissage = pow(-I, L) * (2 * L + 1.0) * h2_sph(L);
	for (int j = L + 1; j < L_prime + 1; ++j)
		coeff(j) = const_coeff_lissage * pow(cos((j - L) * M_PI / 2.0 / (L_prime - L)), 2);
	return (-I * k / (16.0 * M_PI * M_PI) * (coeff * P_Leg).sum());
}

void make_trans(int Nobs, double r_obs[], complex<double> k,
                int L, int Ndirection, double thetas[], double phis[], complex<double> alpha[])
{
	typedef Array<double, Dynamic, Dynamic, RowMajor> ArrayXXd_RowMajor;
	ArrayXXd Robs = Map<ArrayXXd_RowMajor>(r_obs, Nobs, 3);
	int index = 0;
	for (int obs = 0; obs < Nobs; obs++)
	{
		double norm_robs = 0., r_mn[3];
		for (int i = 0; i < 3; i++)
		{
			norm_robs += Robs(obs, i) * Robs(obs, i);
			r_mn[i] = Robs(obs, i);
		}
		norm_robs = sqrt(norm_robs);
		ArrayXcd h2_sph(L + 1);
		complex<double> z = k * norm_robs;
		for (int i = 0; i < L + 1; i++) h2_sph(i) = boost::math::sph_hankel_2(i, real(z));
		for (int dir = 0; dir < Ndirection; dir++)
		{
			alpha[index] = alpha_computation(thetas[dir], phis[dir], h2_sph, r_mn, L, L, k);
			index++;
		}
	}
}

Matrix3d make_rotMat(double alpha, double beta, double gamma)
{
	Matrix3d rotMat_alpha, rotMat_beta, rotMat_gamma;
	rotMat_alpha << 1 , 0 , 0 ,
	             0 , cos(alpha) , sin(alpha) ,
	             0 , -sin(alpha) , cos(alpha);
	rotMat_beta << cos(beta) , 0 , -sin(beta) ,
	            0 , 1 , 0 ,
	            sin(beta) , 0 , cos(beta);
	rotMat_gamma << cos(gamma) , sin(gamma) , 0 ,
	             -sin(gamma) , cos(gamma) , 0 ,
	             0 , 0 , 1;
	return rotMat_alpha * rotMat_beta * rotMat_gamma;
}

void make_ff(int Npol, double polAngle[],
             int Ntheta0, int Nphi0, double theta0[], double phi0[],
             complex<double> ff0_th[], complex<double> ff0_ph[],
             int Ndirection, double thetas[], double phis[],
             complex<double> ff_th[], complex<double> ff_ph[])
{
	for (int j = 0; j < Ndirection; j++)
	{
		double khat[3] = {sin(thetas[j]) * cos(phis[j]), sin(thetas[j]) * sin(phis[j]), cos(thetas[j])};
		double th_hat[3] = {cos(thetas[j]) * cos(phis[j]), cos(thetas[j]) * sin(phis[j]), -sin(thetas[j])};
		double ph_hat[3] = { -sin(phis[j]), cos(phis[j]), 0.};

		for (int i = 0; i < Npol; i++)
		{
			double alpha = polAngle[3 * i], beta = polAngle[3 * i + 1], gamma = polAngle[3 * i + 2];
			Matrix3d rotMat = make_rotMat(alpha, beta, gamma);
			double khat_trans[3] = {0, 0, 0};
			double thhat_trans[3] = {0, 0, 0}, phhat_trans[3] = {0, 0, 0};
			for (int m = 0; m < 3; m++)
			{
				for (int n = 0; n < 3; n++)
				{
					khat_trans[m] -= rotMat(m, n) * khat[n];
					thhat_trans[m] += rotMat(m, n) * th_hat[n];
					phhat_trans[m] += rotMat(m, n) * ph_hat[n];
				}
			}
			double phi1 = atan2(khat_trans[1], khat_trans[0]), theta1 = acos(khat_trans[2]);
			double th_hat1[3] = {cos(theta1) * cos(phi1), cos(theta1) * sin(phi1), -sin(theta1)};
			double ph_hat1[3] = { -sin(phi1), cos(phi1), 0.};

			double th1_dot_th = th_hat1[0] * thhat_trans[0]
			                    + th_hat1[1] * thhat_trans[1] + th_hat1[2] * thhat_trans[2];
			double th1_dot_ph = th_hat1[0] * phhat_trans[0]
			                    + th_hat1[1] * phhat_trans[1] + th_hat1[2] * phhat_trans[2];
			double ph1_dot_th = ph_hat1[0] * thhat_trans[0]
			                    + ph_hat1[1] * thhat_trans[1] + ph_hat1[2] * thhat_trans[2];
			double ph1_dot_ph = ph_hat1[0] * phhat_trans[0]
			                    + ph_hat1[1] * phhat_trans[1] + ph_hat1[2] * phhat_trans[2];

			int index = i * Ndirection + j;
			complex<double> ffth1_tmp = ff_interpolation2D(Ntheta0, Nphi0, theta0, phi0, ff0_th,
			                            theta1, phi1);
			complex<double> ffph1_tmp = ff_interpolation2D(Ntheta0, Nphi0, theta0, phi0, ff0_ph,
			                            theta1, phi1);
			ff_th[index] = th1_dot_th * ffth1_tmp + ph1_dot_th * ffph1_tmp;
			ff_ph[index] = th1_dot_ph * ffth1_tmp + ph1_dot_ph * ffph1_tmp;
		}
	}
}
