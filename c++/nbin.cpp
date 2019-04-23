#include "nbin.hpp"

void nbin(const ArrayXXcd& A,
          const int MAXITER,
          const double TOL,
          ArrayXd& dl,
          ArrayXd& dr)
{
    const int m = A.rows(), n = A.cols();
    const double r = sqrt((double) n / (double) m);
    double err = 1;
    ArrayXd dl0(m), dl1 = ArrayXd::Zero(m), dr0(n), dr1 = ArrayXd::Zero(n);
    ArrayXXcd Ak(m, n);

    dl = ArrayXd::Ones(m);
    dr = ArrayXd::Ones(n);
    Ak = A;

    for (int k = 0; k < MAXITER; k++)
    {
        dl0 = dl1;
        for (int i = 0; i < m; i++)
        {
            ArrayXcd A_tmp = Ak.row(i);
            dl1(i) = sqrt(1 / norm2(A_tmp));
            dl(i) *= dl1(i);
        }

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Ak(i, j) = dl(i) * A(i, j) * dr(j);
            }
        }

        dr0 = dr1;
        for (int j = 0; j < n; j++)
        {
            ArrayXcd A_tmp = Ak.col(j);
            dr1(j) = sqrt(r / norm2(A_tmp));
            dr(j) *= dr1(j);
        }

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < n; j++)
            {
                Ak(i, j) = dl(i) * A(i, j) * dr(j);
            }
        }

        ArrayXd resi_l(m), resi_r(n);
        resi_l = dl1 - dl0;
        resi_r = dr1 - dr0;
        err = norm2(resi_l) / norm2(dl1)
              + norm2(resi_r) / norm2(dr1);

        if (err <= TOL)
        {
            cout << "converge! " << "err = " << err << ", "
                 << "iter = " << k + 1 << endl;
            break;
        }
    }
}

void wrap_nbin(complex<double> A_tmp[],
               int Nrow, int Ncol, int MAXITER, double TOL,
               double dl_tmp[], double dr_tmp[])
{
    typedef Array<complex<double>, Dynamic, Dynamic, RowMajor> ArrayXXcd_RowMajor;
    ArrayXXcd A = Map<ArrayXXcd_RowMajor>(A_tmp, Nrow, Ncol);
    ArrayXd dl(Nrow), dr(Ncol);

    nbin(A, MAXITER, TOL, dl, dr);

    for (int i = 0; i < Nrow; i++) dl_tmp[i] = dl(i);
    for (int j = 0; j < Ncol; j++) dr_tmp[j] = dr(j);
}