#include <Eigen/Dense>
using namespace Eigen;
#include "interpolation.hpp"

double dilichlet(int n, double x)
{
	int nh = n / 2;
	double d, denom = sin(0.5 * x);
	if (abs(denom) < 1e-10)
	{
		d = double(2 * nh) + 1.;
	}
	else
	{
		d = sin((nh + 0.5) * x) / denom;
	}
	return d;
}

complex<double> interpolation1D(int Nx0, double x0[], complex<double> f0[], double x1)
{
	complex<double> f1 = 0;
	for (int j = 0; j < Nx0; j++)
	{
		double D = dilichlet(Nx0, x1 - x0[j]) / Nx0;
		f1 += D * f0[j];
	}
	return f1;
}

complex<double> ff_interpolation2D(int Ntheta0, int Nphi0, double theta0[], double phi0[],
                                   complex<double> ff0[],
                                   double theta1, double phi1)
{
	double* th0_round = new double[2 * Ntheta0];
	for (int j = 0; j < Ntheta0; j++)
	{
		th0_round[j] = theta0[j];
	}
	for (int j = Ntheta0; j < 2 * Ntheta0; j++)
	{
		th0_round[j] = 2 * M_PI - theta0[2 * Ntheta0 - j - 1];
	}
	complex<double>* fthn = new complex<double>[2 * Ntheta0];
	for (int j = 0; j < Ntheta0; j++)
	{
		int ind_start = j * Nphi0;
		fthn[j] = interpolation1D(Nphi0, phi0, &ff0[ind_start], phi1);
		fthn[2 * Ntheta0 - j - 1]
		    = -interpolation1D(Nphi0, phi0, &(ff0[ind_start]), phi1 + M_PI);
	}
	return interpolation1D(2 * Ntheta0, th0_round, fthn, theta1);
}
