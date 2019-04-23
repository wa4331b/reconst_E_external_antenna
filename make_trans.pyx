import numpy as np

cdef extern from "c++/make_trans.hpp":
    cdef void make_trans(int Nobs, double r_obs[], complex k,
                         int L, int Ndirection, double thetas[], double phis[], complex alpha[])
    cdef void make_ff(int Npol, double polAngle[],
                      int Ntheta0, int Nphi0, double theta0[], double phi0[],
                      complex ff0_th[], complex ff0_ph[],
                      int Ndirection, double thetas[], double phis[],
                      complex ff_th[], complex ff_ph[])

def Pywrap_makeTrans(double[:, ::1] r_obs, complex k, int L, double[::1] thetas, double[::1] phis):
    cdef int Nobs = r_obs.shape[0], Ndirection = thetas.size
    cdef complex[:, ::1] alpha = np.zeros([Nobs, Ndirection], dtype=complex)
    make_trans(Nobs, &r_obs[0, 0], k, L, Ndirection, &thetas[0], &phis[0], &alpha[0, 0])

    return np.asarray(alpha)

def Pywrap_makeFFforTrans(double[:, ::1] polAngle, double[::1] theta0, double[::1] phi0,
                          complex[::1] ff0_th, complex[::1] ff0_ph,
                          double[::1] thetas, double[::1] phis):

    cdef int Npol = polAngle.shape[0], Ntheta0 = theta0.size, Nphi0 = phi0.size
    cdef int Ndirection = thetas.size
    cdef complex[:, ::1] ff_th = np.zeros([Npol, Ndirection], dtype=complex)
    cdef complex[:, ::1] ff_ph = np.zeros([Npol, Ndirection], dtype=complex)
    make_ff(Npol, &polAngle[0, 0], Ntheta0, Nphi0, &theta0[0], &phi0[0], &ff0_th[0], &ff0_ph[0],
            Ndirection, &thetas[0], &phis[0], &ff_th[0, 0], &ff_ph[0, 0])

    return np.asarray(ff_th), np.asarray(ff_ph)