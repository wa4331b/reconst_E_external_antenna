#######################################################################
##
## read_MLFMA_mesh_part1.py
##
## Copyright (C) 2014 Idesbald Van den Bosch
##
## This file is part of Puma-EM.
## 
## Puma-EM is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## Puma-EM is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Puma-EM.  If not, see <http://www.gnu.org/licenses/>.
##
## Suggestions/bugs : <vandenbosch.idesbald@gmail.com>
##
#######################################################################

import sys, os, argparse, time
from ReadWriteBlitzArray import writeScalarToDisk, writeBlitzArrayToDisk
from read_mesh import read_mesh_GMSH_2
import numpy as np

def setup_mesh(params_simu, simuDirName, inputDirName):
    """Sets up the mesh.
       params_simu is a class instance that contains the parameters for the simulation.
    """
    tmpDirName = simuDirName
    geoDirName = os.path.join(simuDirName, 'geo')
    # size of cube at finest level

    # reading the mesh
    path, name = geoDirName, params_simu.meshFileName
    meshPath = os.path.join(tmpDirName, "mesh")
    print("read_MLFMA_mesh_part1.py: reading" + os.path.join(path, name) + "...")
    vertexes_coord, triangle_vertexes, triangles_physicalSurface \
        = read_mesh_GMSH_2(os.path.join(path, name), 1., 0.)
    T = triangle_vertexes.shape[0]
    V = vertexes_coord.shape[0]
    print("number of triangles = " + str(T) + "\n")

    writeScalarToDisk(T, os.path.join(meshPath, "T.txt"))
    writeScalarToDisk(V, os.path.join(meshPath, "V.txt"))
    writeBlitzArrayToDisk(vertexes_coord, os.path.join(meshPath, 'vertexes_coord.txt'))
    writeBlitzArrayToDisk(triangle_vertexes, os.path.join(meshPath, 'triangle_vertexes.txt'))

    print("Reading input files in" + params_simu.inputFileDir)
    v = np.loadtxt(inputDirName + '/inputFiles/' + params_simu.inputFileDir + "/voltages.txt")
    v = v[:, 0] + 1j * v[:, 1]
    rObs = np.loadtxt(inputDirName + '/inputFiles/' + params_simu.inputFileDir + "/rObs.txt")
    polAngles = np.loadtxt(inputDirName + '/inputFiles/' + params_simu.inputFileDir + "/polAngles.txt")
    polNumbers = np.loadtxt(
        inputDirName + '/inputFiles/' + params_simu.inputFileDir + "/polNumbers.txt",
        dtype="i"
    )
    rObsNumbers = np.loadtxt(
        inputDirName + '/inputFiles/' + params_simu.inputFileDir + "/rObsNumbers.txt",
        dtype="i"
    )

    np.save(meshPath+"/vObs.npy", v)
    np.save(meshPath+"/polAngles.npy", polAngles)
    np.save(meshPath+"/r_obs.npy", rObs)
    np.save(meshPath+"/polNumbers.npy", polNumbers)
    np.save(meshPath+"/rObsNumbers.npy", rObsNumbers)

def run(simuDirName, inputDirName):
    from inputParams.simulation_parameters import params_simu
    setup_mesh(params_simu, simuDirName, inputDirName)


if __name__ == "__main__":
    simuDirName = "."
    inputDirName = "./inputParams"
    run(simuDirName, inputDirName)
