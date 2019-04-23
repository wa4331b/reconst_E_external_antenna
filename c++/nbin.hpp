#include <iostream>
#include <Eigen/Dense>
// #include <random>
using namespace std;
using namespace Eigen;

inline ArrayXf
conjArray(const ArrayXf& x)
{
    return x;
}

inline ArrayXd
conjArray(const ArrayXd& x)
{
    return x;
}

template <class T>
inline Array<std::complex<T>, Dynamic, 1>
conjArray(const Array<std::complex<T>, Dynamic, 1>& x)
{
    Array<std::complex<T>, Dynamic, 1> y(x.size());
    y = conj(x);
    return y;
}

template <class T>
double norm2(const Array<T, Dynamic, 1>& x)
{
    return sqrt(abs((x * conjArray(x)).sum()));
}

template <class T>
double squareNorm2(const Array<T, Dynamic, 1>& x)
{
    return abs((x * conjArray(x)).sum());
}

void nbin(const ArrayXXcd& A,
          const int MAXITER,
          const double TOL,
          ArrayXd& dl,
          ArrayXd& dr);

void wrap_nbin(complex<double> A_tmp[],
               int Nrow, int Ncol, int MAXITER, double TOL,
               double dl_tmp[], double dr_tmp[]);