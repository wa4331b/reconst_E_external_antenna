import sys, time, os

try:
    import cPickle
except ImportError:
    import pickle as cPickle
from meshClass import MeshClass
from MoMPostProcessing import *
from ReadWriteBlitzArray import *
from scipy import cos, sin, conj, log10, real, sum, dot, pi, sqrt, exp


def computeCurrentsVisualization(resultPath, params_simu, I_coeff,
                                 filename_norm = None, filename_vec = None,
                                 target_surface = np.asarray([0, 1])):

    target_mesh \
        = MeshClass(params_simu.targetName, 1., 0.)

    # target_mesh.constructFromGmshFile()
    target_mesh.constructFromSavedArrays("./mesh")
    N_RWG = readIntFromDisk("./mesh/N_RWG.txt")

    print("............Constructing visualisation of currents")

    J_centroids_triangles = JMCentroidsTriangles(I_coeff, target_mesh)

    norm_J_centroids_triangles \
        = normJMCentroidsTriangles(J_centroids_triangles)

    if filename_vec == None:
        filename = params_simu.targetName + '.J_centroids_triangles.pos'
    else:
        filename = filename_vec
    write_VectorFieldTrianglesCentroidsGMSH(resultPath + filename,
                                            real(J_centroids_triangles), target_mesh,
                                            target_surface)

    if filename_norm == None:
        filename = params_simu.targetName + '.norm_J_centroids_triangles.pos'
    else:
        filename = filename_norm
    write_ScalarFieldTrianglesCentroidsGMSH(resultPath + filename,
                                            norm_J_centroids_triangles, target_mesh,
                                            target_surface)

    print("............end of currents visualisation construction. "
          + "Open with gmsh the *.pos files in result directory.")