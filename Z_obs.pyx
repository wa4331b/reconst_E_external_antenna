import numpy as np

cdef extern from "c++/compute_Z_obs.hpp":
    
    cdef void wrap_compute_Z_obs(complex Z_EJ_obs_tmp[],
                                 complex Z_EM_obs_tmp[],
                                 complex Z_HJ_obs_tmp[],
                                 complex Z_HM_obs_tmp[],
                                 int N_obs,
                                 double r_obs_tmp[],
                                 int N_RWG,
                                 double RWGNumber_trianglesCoord_tmp[],
                                 int T,
                                 int triangle_NRWG_tmp[],
                                 int triangle_RWGNumber_tmp[],
                                 int triangle_signInRWG_tmp[],
                                 int triangle_surfaces_tmp[],
                                 double w, complex eps_r,
                                 complex mu_r,
                                 int N_target_surface,
                                 int target_surface_tmp[])

def Pywrap_compute_Z_obs(double[:, ::1] r_obs,
                         double[:, ::1] RWGNumber_trianglesCoord,
                         int[::1] triangle_NRWG, int[:, ::1] triangle_RWGNumber,
                         int[:, ::1] triangle_signInRWG, int[::1] triangle_surface,
                         double w, complex eps_r, complex mu_r, int[::1] target_surface):

    cdef int N_obs = r_obs.shape[0]
    cdef int N_RWG = RWGNumber_trianglesCoord.shape[0]
    cdef int T = triangle_NRWG.size

    cdef complex[:, :, ::1] Z_EJ_obs = np.empty([N_obs, 3, N_RWG], dtype=complex)
    cdef complex[:, :, ::1] Z_EM_obs = np.empty([N_obs, 3, N_RWG], dtype=complex)
    cdef complex[:, :, ::1] Z_HJ_obs = np.empty([N_obs, 3, N_RWG], dtype=complex)
    cdef complex[:, :, ::1] Z_HM_obs = np.empty([N_obs, 3, N_RWG], dtype=complex)

    wrap_compute_Z_obs(&Z_EJ_obs[0, 0, 0], &Z_EM_obs[0, 0, 0],
                       &Z_HJ_obs[0, 0, 0], &Z_HM_obs[0, 0, 0],
                       N_obs, &r_obs[0, 0], N_RWG,
                       &RWGNumber_trianglesCoord[0, 0], T, &triangle_NRWG[0],
                       &triangle_RWGNumber[0, 0], &triangle_signInRWG[0, 0],
                       &triangle_surface[0],
                       w, eps_r, mu_r, target_surface.size, &target_surface[0])

    return np.asarray(Z_EJ_obs), np.asarray(Z_EM_obs), \
           np.asarray(Z_HJ_obs), np.asarray(Z_HM_obs)
