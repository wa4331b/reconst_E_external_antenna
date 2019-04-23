import numpy as np

cdef extern from "c++/nbin.hpp":
    cdef void wrap_nbin(complex A_tmp[],
                        int Nrow, int Ncol, int MAXITER, double TOL,
                        double dl_tmp[], double dr_tmp[])

def Pywrap_nbin(complex[:, ::1] A, int MAXITER, double TOL):

    cdef int Nrow = A.shape[0]
    cdef int Ncol = A.shape[1]
    cdef double[::1] dl = np.empty(Nrow, dtype="d")
    cdef double[::1] dr = np.empty(Ncol, dtype="d")

    wrap_nbin(&A[0, 0], Nrow, Ncol, MAXITER, TOL, &dl[0], &dr[0])

    return np.asarray(dl), np.asarray(dr)