import numpy as np

cdef extern from "c++/poynting.hpp":
    cdef void poyntingMatrix(int N_RWG,
                             double RWGNumber_trianglesCoord_tmp[],
                             int T,
                             int triangle_NRWG[],
                             int triangle_RWGNumber_tmp[],
                             int triangle_signInRWG_tmp[],
                             int triangle_surfaces[],
                             int N_target_surface,
                             int target_surface_tmp[],
                             complex poyntingMatrix_tmp[])

def Pywrap_poynting(double[:, ::1] RWGNumber_trianglesCoord,
                    complex[::1] I_J, complex[::1] I_M,
                    int[::1] triangle_NRWG, int[:, ::1] triangle_RWGNumber,
                    int[:, ::1] triangle_signInRWG, int[::1] triangle_surface,
                    int[::1] target_surface):
    matrix = Pywrap_poyntingMatrix(RWGNumber_trianglesCoord, triangle_NRWG, triangle_RWGNumber,
                                   triangle_signInRWG, triangle_surface, target_surface)

    return np.vdot(I_J, np.dot(matrix, I_M))

def Pywrap_poyntingMatrix(double[:, ::1] RWGNumber_trianglesCoord,
                          int[::1] triangle_NRWG, int[:, ::1] triangle_RWGNumber,
                          int[:, ::1] triangle_signInRWG, int[::1] triangle_surface,
                          int[::1] target_surface):
    cdef int numberOfRWG = RWGNumber_trianglesCoord.shape[0]
    cdef complex[:, ::1] matrix = np.zeros((numberOfRWG, numberOfRWG), dtype=complex)
    poyntingMatrix(numberOfRWG, &RWGNumber_trianglesCoord[0, 0],
                   triangle_NRWG.shape[0], &triangle_NRWG[0],
                   &triangle_RWGNumber[0, 0], &triangle_signInRWG[0, 0],
                   &triangle_surface[0], target_surface.size, &target_surface[0], &matrix[0, 0])
    return np.asarray(matrix)
