#######################################################################
##
## meshClass.py
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

import os.path, sys, time
try:
    import cPickle
except ImportError:
    import pickle as cPickle

from scipy import zeros, ones, arange, array, take, reshape, sort, argsort, put
from scipy import mean, sqrt, sum, compress, nonzero, prod
from ReadWriteBlitzArray import *
import copy

class MeshClass:
    
    def __init__(self, targetName, targetDimensions_scaling_factor, z_offset):

        self.targetName = targetName
        self.z_offset = z_offset
        self.targetDimensions_scaling_factor = targetDimensions_scaling_factor

    def constructFromSavedArrays(self, path):
        
        self.E = readIntFromDisk(os.path.join(path, "N_RWG.txt"))
        self.T = readIntFromDisk(os.path.join(path, "T.txt"))
        self.V = readIntFromDisk(os.path.join(path, "V.txt"))
        self.N_RWG = self.E
        # now the arrays
        self.vertexes_coord \
            = readBlitzArrayFromDisk(os.path.join(path, "vertexes_coord.txt"), self.V, 3, 'd')
        self.triangle_vertexes \
            = readBlitzArrayFromDisk(os.path.join(path, "triangle_vertexes.txt"), self.T, 3, 'i')
        self.triangles_surfaces \
            = read1DBlitzArrayFromDisk(os.path.join(path, "triangles_surfaces.txt"), 'i')
        self.RWGNumber_signedTriangles \
            = readBlitzArrayFromDisk(os.path.join(path, "RWGNumber_signedTriangles.txt"),
                                     self.E, 2, 'i')
        self.RWGNumber_edgeVertexes \
            = readBlitzArrayFromDisk(os.path.join(path, "RWGNumber_edgeVertexes.txt"),
                                     self.E, 2, 'i')
        self.RWGNumber_oppVertexes \
            = readBlitzArrayFromDisk(os.path.join(path, "RWGNumber_oppVertexes.txt"),
                                     self.E, 2, 'i')