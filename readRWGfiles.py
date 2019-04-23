from ReadWriteBlitzArray import *
import numpy as np


class MeshData:
    def __init__(self, mesh_path):
        self.N_RWG = readIntFromDisk(mesh_path + "N_RWG.txt")

        RWGNumber_oppVertexes = readBlitzArrayFromDisk(mesh_path + "RWGNumber_oppVertexes.txt",
                                                       self.N_RWG, 2, "i")
        RWGNumber_edgeVertexes = readBlitzArrayFromDisk(mesh_path + "RWGNumber_edgeVertexes.txt",
                                                        self.N_RWG, 2, "i")

        self.RWGNumber_nodes = np.zeros([self.N_RWG, 4], dtype="i")
        self.RWGNumber_nodes[:, 0] = RWGNumber_oppVertexes[:, 0]
        self.RWGNumber_nodes[:, 1] = RWGNumber_edgeVertexes[:, 0]
        self.RWGNumber_nodes[:, 2] = RWGNumber_edgeVertexes[:, 1]
        self.RWGNumber_nodes[:, 3] = RWGNumber_oppVertexes[:, 1]

        self.V = readIntFromDisk(mesh_path + "V.txt")
        self.nodesCoord = readBlitzArrayFromDisk(mesh_path + "vertexes_coord.txt",
                                                 self.V, 3, float)

        self.RWGNumber_trianglesCoord = np.empty([self.N_RWG, 12], dtype=float)
        self.RWGNumber_trianglesCoord[:, :3] = self.nodesCoord[self.RWGNumber_nodes[:, 0], :]
        self.RWGNumber_trianglesCoord[:, 3:6] = self.nodesCoord[self.RWGNumber_nodes[:, 1], :]
        self.RWGNumber_trianglesCoord[:, 6:9] = self.nodesCoord[self.RWGNumber_nodes[:, 2], :]
        self.RWGNumber_trianglesCoord[:, 9:] = self.nodesCoord[self.RWGNumber_nodes[:, 3], :]

        self.RWGNumber_surface \
            = read1DBlitzArrayFromDisk(mesh_path + "RWGNumber_surfaces.txt", "i")

        self.numberOfSurfaces = np.max(self.RWGNumber_surface) + 1

        self.triangle_NRWG = np.load(mesh_path + "triangle_NRWG.npy").astype("i")
        self.triangle_RWGNumber = np.load(mesh_path + "triangle_RWGNumber.npy").astype("i")
        self.triangle_signInRWG = np.load(mesh_path + "triangle_signInRWG.npy").astype("i")
        self.triangle_surface = read1DBlitzArrayFromDisk(mesh_path + "triangles_surfaces.txt", "i")

        self.RWGNumber_signedTriangles \
            = readBlitzArrayFromDisk(mesh_path + "RWGNumber_signedTriangles.txt",
                                     self.N_RWG, 2, "i")

    def surface_rwg(self):
        surface_RWGnumbers = []
        surface_NRWG = []
        for i in range(self.numberOfSurfaces):
            surface_RWGnumbers.append(np.where(self.RWGNumber_surface == i)[0].astype("i"))
            surface_NRWG.append(surface_RWGnumbers[i].size)

        return surface_RWGnumbers, surface_NRWG
