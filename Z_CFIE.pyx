import numpy as np
from libcpp cimport bool

cdef extern from "c++/Z_EJ_Z_HJ.hpp":

    cdef void wrap_Z_CFIE_computation(complex Z_CFIE_EJ_tmp[], complex Z_CFIE_EM_tmp[],
                                      complex Z_CFIE_HJ_tmp[], complex Z_CFIE_HM_tmp[],
                                      complex CFIE[], int CFIE_bool_tmp[], int M_current_tmp,
                                      double signSurfObs, double signSurfSrc, 
                                      int N_RWG, int N_RWG_src, int N_RWG_test, int V,
                                      int numbers_RWG_src[], int numbers_RWG_test[],
                                      int srcRWGNumber_surface[],
                                      int RWGNumber_signedTriangles_tmp[],
                                      int RWGNumber_nodes_tmp[], double nodesCoord_tmp[],
                                      double w, complex eps_r, complex mu_r,
                                      int N_target_surface, int target_surface_tmp[])

    
    cdef void wrap_Z_CFIE_test_computation(complex Z_CFIE_tmp[],
                                           int N_RWG, int V,
                                           int srcRWGNumber_surface[],
                                           int RWGNumber_signedTriangles_tmp[],
                                           int RWGNumber_nodes_tmp[],
                                           double nodesCoord_tmp[],
                                           double w,
                                           complex eps_r,
                                           complex mu_r,
                                           int N_target_surface, int target_surface_tmp[])

def Pywrap_Z_CFIE_computation(complex[::1] CFIE, int[::1] CFIE_bool, int M_current,
                              double signSurfObs, double signSurfSrc,
                              int[::1] numbers_RWG_src,
                              int[::1] numbers_RWG_test,
                              int[::1] srcRWGNumber_surface,
                              int[:, ::1] RWGNumber_signedTriangles,
                              int[:, ::1] RWGNumber_nodes,
                              double[:, ::1] nodesCoord,
                              double w, complex eps_r, complex mu_r,
                              int[::1] target_surface):

    cdef int N_RWG = srcRWGNumber_surface.size
    cdef int N_RWG_src = numbers_RWG_src.size
    cdef int N_RWG_test = numbers_RWG_test.size
    cdef int V = nodesCoord.shape[0]
    cdef int N_target_surface = target_surface.size
    cdef complex[:, ::1] Z_CFIE_EJ = np.zeros([N_RWG_test, N_RWG_src], dtype=complex)
    cdef complex[:, ::1] Z_CFIE_EM = np.zeros([N_RWG_test, N_RWG_src], dtype=complex)
    cdef complex[:, ::1] Z_CFIE_HJ = np.zeros([N_RWG_test, N_RWG_src], dtype=complex)
    cdef complex[:, ::1] Z_CFIE_HM = np.zeros([N_RWG_test, N_RWG_src], dtype=complex)

    wrap_Z_CFIE_computation(&Z_CFIE_EJ[0, 0], &Z_CFIE_EM[0, 0],
                            &Z_CFIE_HJ[0, 0], &Z_CFIE_HM[0, 0],
                            &CFIE[0], &CFIE_bool[0], M_current,
                            signSurfObs, signSurfSrc, N_RWG, N_RWG_src, N_RWG_test, V,
                            &numbers_RWG_src[0], &numbers_RWG_test[0],
                            &srcRWGNumber_surface[0], &RWGNumber_signedTriangles[0, 0],
                            &RWGNumber_nodes[0, 0], &nodesCoord[0, 0],
                            w, eps_r, mu_r, N_target_surface, &target_surface[0])

    return np.asarray(Z_CFIE_EJ), np.asarray(Z_CFIE_EM), \
           np.asarray(Z_CFIE_HJ), np.asarray(Z_CFIE_HM)

def Pywrap_wrap_Z_CFIE_test_computation(int[::1] srcRWGNumber_surface,
                                        int[:, ::1] RWGNumber_signedTriangles,
                                        int[:, ::1] RWGNumber_nodes,
                                        double[:, ::1] nodesCoord,
                                        double w, complex eps_r, complex mu_r,
                                        int[::1] target_surface):

    cdef int N_RWG = srcRWGNumber_surface.size
    cdef int V = nodesCoord.shape[0]
    cdef int N_target_surface = target_surface.size
    cdef complex[:, ::1] Z_CFIE = np.zeros([N_RWG, N_RWG], dtype=complex)
    wrap_Z_CFIE_test_computation(&Z_CFIE[0, 0], N_RWG, V,
                                 &srcRWGNumber_surface[0], &RWGNumber_signedTriangles[0, 0],
                                 &RWGNumber_nodes[0, 0], &nodesCoord[0, 0],
                                 w, eps_r, mu_r, N_target_surface, &target_surface[0])

    return np.asarray(Z_CFIE)
