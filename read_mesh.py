#######################################################################
##
## read_mesh.py
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

import os, sys, time
from scipy import zeros, size, compress, sort, take, put, array
from PyGmsh import executeGmsh, write_geo
#from scipy import weave
#from scipy.weave import converters

def preRead_mesh_GMSH_2(meshFile_name):
    """
       this function splits the meshFile_name.msh file into smaller entities.
       How do we do it?

       Each time we encounter a new 'entity' in the mesh file, such as 'Nodes' 
       for example, we create a new file for it, name it with the entity name,
       say 'meshFile_name.msh.Nodes', and place there all corresponding entities,
       i.e. Nodes in this example.
    """
    file = open(meshFile_name, 'r')
    content = {}
    # the special fields that can be encountered
    specialFields = ['PhysicalNames', 'Nodes', 'NOD', 'Coordinates', 'Elements', 'ELM']
    for line in file:
        # here we encounter a new entity, or field
        if line[0]=='$' and not ('End' in line or 'END' in line):
            if line.split('\n')[0] not in content:
                newFieldTmp = line.split('\n')[0]
                # we actualise newField
                newField = newFieldTmp.split('$')[1]
                content[newField] = []
                # we create the file corresponding to the entity
                fileForField = open(meshFile_name + '.' + newField, 'w')
                numberOfLines = 0
        # here we encounter the end of the new entity, or field
        elif newField in line and ('End' in line or 'END' in line):
            fileForField.close()
            content[newField] = [numberOfLines]
            pass
        else:
            if (len(content[newField])==0) and (newField in specialFields):
                content[newField] = [0]
                pass
            elif (newField in ['Elements', 'ELM']):
                elementTmp2 = line.split()
                if elementTmp2[1]=='2': # type 'planar triangle' has number 2
                  fileForField.write(line)
                  numberOfLines += 1
            else:
                fileForField.write(line)
                numberOfLines += 1
    file.close()
    return content

def read_mesh_GMSH_2(name, targetDimensions_scaling_factor, z_offset):
    """function that reads the mesh and puts it into nice arrays"""
    content = preRead_mesh_GMSH_2(name)
    print(content)
    if 'Nodes' in content:
        vertexes_key = 'Nodes'
    elif 'NOD' in content:
        vertexes_key = 'NOD'
    else:
        print("read_mesh.py: Error in the GMSH file format: not supported!!")
        sys.exit(1)
    V = int(content[vertexes_key][0]) # we get the number of vertexes
    vertexes_numbers = zeros(V, 'i')
    vertexes_coord = zeros( (V, 3), 'd' )

    index = 0
    file = open(name + '.' + vertexes_key, 'r')
    for line in file:
        tmp = line.split()
        vertexes_numbers[index] = int(tmp[0])
        vertexes_coord[index, :] = list(map(float, tmp[1:]))
        index += 1
    if not (targetDimensions_scaling_factor==1.0):
        vertexes_coord *= targetDimensions_scaling_factor
    vertexes_coord[:, -1] += z_offset
    file.close()

    # we check whether we have a source
    SOURCE = False
    PhysicalSurfaceNumberToName = {}
    if 'PhysicalNames' in content:
        keyPhysicalNames = 'PhysicalNames'
        NumberPhysicalNames = content[keyPhysicalNames][0]
        file = open(name + '.' + keyPhysicalNames, 'r')
        for line in file:
            elem2 = line.split()
            newKey = int(elem2[0])
            PhysicalSurfaceNumberToName[newKey] = elem2[1]
            if "source" in PhysicalSurfaceNumberToName[newKey]:
                SOURCE = True
        file.close()

    # we now extract the triangles and write them to a file
    if 'Elements' in content:
        elements_key = 'Elements'
    elif 'ELM' in content:
        elements_key = 'ELM'
    else:
        print("read_mesh.py: Error in the GMSH file format: not supported!!")
        sys.exit(1)
    T = content[elements_key][0] # N is the number of triangles
    del content

    g = open(name + "." + elements_key, 'r')
    triangles_nodes = zeros( (T, 3), 'i')
    triangles_physicalSurface = zeros(T, 'i')
    index = 0
    for line in g:
        tmp = list(map(int, line.split()))
        triangles_nodes[index, :] = tmp[-3:]
        triangles_physicalSurface[index] = tmp[3]
        index += 1
    g.close()
    
    # indexes of elements in Python/C++ arrays start at 0.
    # However, triangles_nodes don't necessarily start at 0.
    # So the following 4 lines correct that.
    vertex_number_max = max(vertexes_numbers)
    nodes_vertexes = zeros( vertex_number_max + 1, 'i' )
    put(nodes_vertexes, vertexes_numbers, range(V))
    triangles_vertexes = take(nodes_vertexes, triangles_nodes, axis=0).astype('i')
    # we will now eliminate the points that we don't need
#    del triangles_nodes # gain some memory, we now work only with triangles_vertexes...
#    max_encountered_vertexes = max(triangles_vertexes.flat)
#    encountered_vertexesTmp = (zeros(max_encountered_vertexes + 1, 'i') - 1).astype('i')
#    wrapping_code = """
#    for (int i=0 ; i<triangles_vertexes.extent(0) ; i++) {
#      for (int j=0 ; j<triangles_vertexes.extent(1) ; j++) {
#        const int vertex = triangles_vertexes(i, j);
#        encountered_vertexesTmp(vertex) = vertex;
#      }
#    }
#    """
#    weave.inline(wrapping_code,
#                 ['encountered_vertexesTmp', 'triangles_vertexes'],
#                 type_converters = converters.blitz,
#                 include_dirs = [],
#                 library_dirs = [],
#                 libraries = [],
#                 headers = ['<iostream>'],
#                 compiler = 'gcc',
#                 extra_compile_args = ['-O3', '-pthread', '-w'])
#    encountered_vertexes = compress(encountered_vertexesTmp>-1, encountered_vertexesTmp, axis=0)
#    del encountered_vertexesTmp
#    oldVertex_to_newVertex = zeros(max_encountered_vertexes+1, 'i')
#    oldVertex_to_newVertex[encountered_vertexes] = range(len(encountered_vertexes))
#    triangles_vertexes = take(oldVertex_to_newVertex, triangles_vertexes, axis=0)
#    vertexes_coord = take(vertexes_coord, encountered_vertexes, axis=0)
    return vertexes_coord.astype('d'), triangles_vertexes.astype('i'), triangles_physicalSurface.astype('i')

if __name__=="__main__":
    path = './geo'
    targetName = 'sphere'
    write_geo(path, targetName, 'lc', 0.05)
    write_geo(path, targetName, 'lx', 0.051)
    write_geo(path, targetName, 'ly', 0.1)
    write_geo(path, targetName, 'lz', 0.2)
    write_geo(path, targetName, 'frill_width', 0.02)
    executeGmsh(path, targetName, 0)
    z_offset = 0.0
    targetDimensions_scaling_factor = 1.0
    t0 = time.time()
    vertexes_coord_1, triangles_vertexes_1, triangles_physicalSurface_1 = read_mesh_GMSH_1(os.path.join(path, targetName) + '.msh', targetDimensions_scaling_factor, z_offset)
    print("time for classical *.msh file reading = " + str(time.time() - t0))
    t0 = time.time()
    vertexes_coord_2, triangles_vertexes_2, triangles_physicalSurface_2 = read_mesh_GMSH_2(os.path.join(path, targetName) + '.msh', targetDimensions_scaling_factor, z_offset)
    print("time for new *.msh file reading = " + str(time.time() - t0))
    
    print
    print("difference between python and C++ code. If results different than 0, there is a problem.")
    print(str(sum(abs(vertexes_coord_1 - vertexes_coord_2))))
    print(str(sum(abs(triangles_vertexes_1 - triangles_vertexes_2))))
    print(str(sum(abs(triangles_physicalSurface_1 - triangles_physicalSurface_2))))
    
    vertexes_coord_GiD, triangles_vertexes_GiD, triangles_physicalSurface_GiD = read_mesh_GiD(os.path.join(path, 'aaa1') + '.msh', targetDimensions_scaling_factor, z_offset)

    content = preRead_mesh_ANSYS('./geo')
    vertexes_coord, triangles_vertexes, triangles_physicalSurface = read_mesh_ANSYS(path, "whatever", targetDimensions_scaling_factor, z_offset)
    print(vertexes_coord.shape, triangles_vertexes.shape, triangles_physicalSurface.shape)

