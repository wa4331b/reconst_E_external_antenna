import numpy as np

cdef extern from "c++/compute_FFmat.hpp":
    cdef void compute_FFmat(complex FF_J_th[], complex FF_J_ph[],
                            complex FF_M_th[], complex FF_M_ph[],
                            int N_direction, double thetas[], double phis[],
                            int N_RWG,
                            double RWGNumber_trianglesCoord_tmp[],
                            int T,
                            int triangle_NRWG_tmp[],
                            int triangle_RWGNumber_tmp[],
                            int triangle_signInRWG_tmp[],
                            int triangle_surfaces_tmp[],
                            double w, complex eps_r,
                            complex mu_r)

def Pywrap_compute_FFmat(double[::1] thetas, double[::1] phis,
                         double[:, ::1] RWGNumber_trianglesCoord,
                         int[::1] triangle_NRWG, int[:, ::1] triangle_RWGNumber,
                         int[:, ::1] triangle_signInRWG, int[::1] triangle_surface,
                         double w, complex eps_r, complex mu_r):

    cdef int N_direction = thetas.size
    cdef int N_RWG = RWGNumber_trianglesCoord.shape[0]
    cdef int T = triangle_NRWG.size

    cdef complex[:, ::1] FF_J_th = np.zeros([N_direction, N_RWG], dtype=complex)
    cdef complex[:, ::1] FF_J_ph = np.zeros([N_direction, N_RWG], dtype=complex)
    cdef complex[:, ::1] FF_M_th = np.zeros([N_direction, N_RWG], dtype=complex)
    cdef complex[:, ::1] FF_M_ph = np.zeros([N_direction, N_RWG], dtype=complex)

    compute_FFmat(&FF_J_th[0, 0], &FF_J_ph[0, 0], &FF_M_th[0, 0], &FF_M_ph[0, 0],
                  N_direction, &thetas[0], &phis[0], N_RWG,
                  &RWGNumber_trianglesCoord[0, 0], T,
                  &triangle_NRWG[0], &triangle_RWGNumber[0, 0],
                  &triangle_signInRWG[0, 0], &triangle_surface[0],
                  w, eps_r, mu_r)

    return np.asarray(FF_J_th), np.asarray(FF_J_ph), np.asarray(FF_M_th), np.asarray(FF_M_ph)
