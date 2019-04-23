import sys, os, argparse

sys.path.append('./code')
try:
    import commands
except ImportError:
    import subprocess as commands
try:
    import cPickle
except ImportError:
    import pickle as cPickle
from scipy import take, sum, mean
from ReadWriteBlitzArray import writeScalarToDisk, writeBlitzArrayToDisk
from ReadWriteBlitzArray import readIntFromDisk, readFloatFromDisk, read1DBlitzArrayFromDisk
from ReadWriteBlitzArray import readBlitzArrayFromDisk
import numpy as np


def compute_RWG_CFIE_OK(triangles_surfaces, RWGNumber_signedTriangles,
                        IS_CLOSED_SURFACE):
    RWGNumber_CFIE_OK_tmp1 \
        = take(triangles_surfaces, RWGNumber_signedTriangles, axis=0)
    RWGNumber_CFIE_OK_tmp2 \
        = take(IS_CLOSED_SURFACE, RWGNumber_CFIE_OK_tmp1, axis=0)

    # following the Taskinen et al. paper in PIER
    # we can have a CFIE on a junction straddling a dielectric and metallic body
    RWGNumber_CFIE_OK \
        = ((sum(RWGNumber_CFIE_OK_tmp2, axis=1) >= 1) * 1).astype('i')
    # We cannot have M on a junction between a dielectric and metallic body!
    # The following expression is lacking the fact that a surface can be metallic
    # or dielectric.
    # If metallic, then there is no M current, even if surface is closed
    RWGNumber_M_CURRENT_OK \
        = ((sum(RWGNumber_CFIE_OK_tmp2, axis=1) > 1) * 1).astype('i')

    return RWGNumber_CFIE_OK, RWGNumber_M_CURRENT_OK


def cube_mesh(meshPath, a, vertexes_coord, RWGNumber_edgeCentroidCoord, n_level_max=None):
    from Cubes import cube_lower_coord_computation, RWGNumber_cubeNumber_computation
    from Cubes import cubeIndex_RWGNumbers_computation, findCubeNeighbors

    max_N_cubes_1D, N_levels, big_cube_lower_coord, big_cube_center_coord \
        = cube_lower_coord_computation(a, vertexes_coord, RWGNumber_edgeCentroidCoord)
    if n_level_max is not None:
        N_levels = n_level_max
    RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord \
        = RWGNumber_cubeNumber_computation(a, max_N_cubes_1D, big_cube_lower_coord,
                                           RWGNumber_edgeCentroidCoord)
    cubes_RWGsNumbers, cube_N_RWGs, cubes_centroids \
        = cubeIndex_RWGNumbers_computation(RWGNumber_cubeNumber, RWGNumber_cubeCentroidCoord)
    # cubes_neighborsIndexes, cube_N_neighbors \
    #     = findCubeNeighbors(max_N_cubes_1D, big_cube_lower_coord, cubes_centroids, a)
    C = cubes_centroids.shape[0]

    writeScalarToDisk(N_levels, meshPath + "/N_levels.txt")
    writeScalarToDisk(max_N_cubes_1D, meshPath + "max_N_cubes_1D.txt")
    writeScalarToDisk(C, meshPath + "C.txt")

    np.save(meshPath + 'big_cube_center_coord.npy', big_cube_center_coord)
    np.save(meshPath + 'big_cube_lower_coord.npy', big_cube_lower_coord)
    np.save(meshPath + 'cubes_centroids.npy', cubes_centroids)
    np.save(meshPath + 'cubes_RWGsNumbers.npy', cubes_RWGsNumbers)
    np.save(meshPath + 'cube_N_RWGs.npy', cube_N_RWGs)
    # np.save(meshPath + 'cubes_neighborsIndexes.npy', cubes_neighborsIndexes)
    # np.save(meshPath + 'cube_N_neighbors.npy', cube_N_neighbors)


def setup_mesh(params_simu, c, mesh_path):
    octtree_path = mesh_path + "octtree_data/"

    wave_length = c / params_simu.f
    average_RWG_length = readFloatFromDisk(os.path.join(mesh_path, 'average_RWG_length.txt'))
    print("average RWG length = " + str(average_RWG_length)
          + "m = lambda /" + str(wave_length / average_RWG_length))
    a = wave_length * params_simu.a_factor
    writeScalarToDisk(a, octtree_path + "leaf_side_length.txt")

    is_closed_surface \
        = read1DBlitzArrayFromDisk(mesh_path + "is_closed_surface.txt", "i")
    S = len(is_closed_surface)
    print("test of the closed surfaces : " + str(is_closed_surface))

    N_RWG = readIntFromDisk(os.path.join(mesh_path, "N_RWG.txt"))
    RWGNumber_edgeVertexes = readBlitzArrayFromDisk(mesh_path + "/RWGNumber_edgeVertexes.txt",
                                                    N_RWG, 2, 'i')
    V = readIntFromDisk(mesh_path + "/V.txt")
    vertexes_coord = readBlitzArrayFromDisk(mesh_path + "/vertexes_coord.txt", V, 3, "d")

    from mesh_functions_seb import compute_RWGNumber_edgeCentroidCoord
    RWGNumber_edgeCentroidCoord \
        = compute_RWGNumber_edgeCentroidCoord(vertexes_coord, RWGNumber_edgeVertexes)
    r_obs = np.load(mesh_path + "/r_obs.npy")

    print("Using good old weave!")
    file_path = octtree_path + "cube_RWG/"
    cube_mesh(file_path, a, vertexes_coord, RWGNumber_edgeCentroidCoord)
    average_NRWG_inCube = np.load(file_path + "cube_N_RWGs.npy").mean()
    cubes_centroids = np.load(file_path + "cubes_centroids.npy")

    n_level_rwg = readIntFromDisk(file_path + "N_levels.txt")
    file_path = octtree_path + "cube_obs/"
    cube_mesh(file_path, a, vertexes_coord, r_obs, n_level_rwg)
    average_Nobs_inCube = np.load(file_path + "cube_N_RWGs.npy").mean()
    cubes_centroids_obs = np.load(file_path + "cubes_centroids.npy")

    print("Searching obs neighbors starting...")
    from Cubes import findCubeNeighbors_obs
    cubes_neighborsIndexes, cube_N_neighbors \
        = findCubeNeighbors_obs(cubes_centroids, cubes_centroids_obs, a)

    n_total_neighbors = sum(cube_N_neighbors)
    if n_total_neighbors != 0:
        print("N_neighbors_obs = ", n_total_neighbors)
        print("Error: Neighbor obs cubes exists. Exit...")
        sys.exit(1)


if __name__ == "__main__":
    from inputParams.simulation_parameters import params_simu
    from inputParams import EM_constants

    meshPath = "./mesh/"
    setup_mesh(params_simu, EM_constants.c, meshPath)
