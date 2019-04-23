#######################################################################
##
## mesh_functions_seb.py
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

import os.path, sys, argparse

sys.path.append('./code')
from scipy import zeros, ones, arange, array, take, reshape, sort, argsort
from scipy import sum, compress, nonzero, sqrt, dot, arccos
# from read_mesh import read_mesh_GMSH_1, read_mesh_GMSH_2
# from PyGmsh import executeGmsh, write_geo, findParameter, findParameterValue
import copy
import numpy as np
from ReadWriteBlitzArray import *


def triangles_centroids_computation(vertexes_coord, triangle_vertexes):
    """guess what? computes the centroids of the triangles"""
    triangles_centroids = take(vertexes_coord, triangle_vertexes[:, 0], axis=0)
    triangles_centroids += take(vertexes_coord, triangle_vertexes[:, 1], axis=0)
    triangles_centroids += take(vertexes_coord, triangle_vertexes[:, 2], axis=0)
    triangles_centroids /= 3.0
    return triangles_centroids


def edges_computation(triangle_vertexes, vertexes_coord):
    """This function builds the edges matrix from the triangles"""

    T = triangle_vertexes.shape[0]  # T is the number of triangles
    E = T * 3  # there are 3 edges per triangles
    edges_vertexes = zeros((E, 3), 'i')  # edges_vertexes[i, :] -> [v_start, v_end, v_opposite]
    # the edge kind is related to its number of triangles it belongs to.
    # 1 is a physical border, 2 is a normal RWG, and 3 or more is a junction

    print("    construction of edges_vertexes...")
    sys.stdout.flush()
    # we first construct a flattened view of edges_vertexes, such that
    # all edges corresponding to 1 triangle are on the same line of the array view
    flat_edges_vertexes = edges_vertexes.reshape((T, -1))

    # we then assign efficiently thanks to the flattened view:
    # edge0 = r0 -> r1, opposite = r2
    # edge1 = r1 -> r2, opposite = r0
    # edge2 = r2 -> r0, opposite = r1
    flat_edges_vertexes[:, 0] = flat_edges_vertexes[:, 5] \
        = flat_edges_vertexes[:, 7] = triangle_vertexes[:, 0]
    flat_edges_vertexes[:, 1] = flat_edges_vertexes[:, 3] \
        = flat_edges_vertexes[:, 8] = triangle_vertexes[:, 1]
    flat_edges_vertexes[:, 2] = flat_edges_vertexes[:, 4] \
        = flat_edges_vertexes[:, 6] = triangle_vertexes[:, 2]

    # For assigning a "number" and a "kind" to an edge (a set of two vertexes),
    # we need to find how many occurrences it has in "edges_vertexes". A costful
    # operation, because for each edge we have to compare the edge to all other
    # that appear in the "edges_vertexes" array.
    #
    # However, if we want efficiency, we can create an array of edges sorted as follows:
    # 1) sort the elements of the first two columns alongside dimension 1; 
    # 2) sort the pairs following their 1st element alongside dimension 0;
    # 3) sort the pairs following their 2nd element alongside dimension 0,
    #    but with keeping the order of their first element.
    # In such an array, all occurrences of an edge would be adjacent to each other. So, for
    # counting the occurrences of an edge, one should only count its similar
    # neighbors, thereby greatly reducing the computational cost of the algorithm

    # we now construct an array of sorted "edges_vertexes" following dimension 1 (the columns) 

    col_sorted_e_v = sort(edges_vertexes[:, :2], 1, kind='mergesort')
    # edges_opp_vertexes = edges_vertexes[:, 2]
    del edges_vertexes

    print("    construction of edgeNumber_triangles...")
    sys.stdout.flush()
    edgeNumber_triangles, edgeNumber_vertexes = compute_edgeNumber_triangles(col_sorted_e_v)
    # construction of triangle_adjacentTriangles matrix
    triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction \
        = compute_triangle_adjacentTriangles(T, edgeNumber_triangles)

    return edgeNumber_vertexes.astype(
        'i'), edgeNumber_triangles, triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction


def compute_edgeNumber_triangles(col_sorted_e_v):
    # Once the elements of the first 2 columns have been sorted alongside dimension 1,
    # we construct a 1-D array of real numbers, with:
    # 1) the entire part formed by the numbers of the first column
    # 2) the decimal part formed by the numbers of the second column
    max_decimal = max(col_sorted_e_v[:, 1])  # the maximum element of "col_sorted_e_v[:, 1]"
    X = 10
    while X < max_decimal:  # we look for smallest "X" such that "1eX > max_decimal"
        X *= 10
    decimal_e_v = col_sorted_e_v[:, 0] * X + col_sorted_e_v[:, 1]
    ind_sorted_e_v = argsort(decimal_e_v, kind='mergesort')
    sorted_decimal_e_v = take(decimal_e_v, ind_sorted_e_v, axis=0)

    diff = ones(sorted_decimal_e_v.shape[0], 'i')
    diff[1:] = abs(sorted_decimal_e_v[1:] - sorted_decimal_e_v[:-1])
    indexesEqualPreceding = compress(diff == 0, arange(len(diff)), axis=0)
    del diff, decimal_e_v, sorted_decimal_e_v

    print("        research of the same edges...")
    sys.stdout.flush()
    indexesEqualEdges = {}
    j = 0
    for i in indexesEqualPreceding:
        indexesEqualEdges[j] = [int(ind_sorted_e_v[i - 1]), int(ind_sorted_e_v[i])]
        if j > 0:
            if indexesEqualEdges[j][0] == indexesEqualEdges[j - 1][-1]:
                # in this case we have a junction: 3 or more triangles that
                # have a common edge. We then merge the two lists and
                # there is no need for increasing the index j
                indexesEqualEdges[j - 1] += indexesEqualEdges[j][1:]
                del indexesEqualEdges[j]
            else:
                j += 1
        else:
            j += 1

    # edge numbering
    print("        numbering of the edges...")
    sys.stdout.flush()
    edgeNumber_vertexes = ones((len(indexesEqualEdges), 2), 'i') * -1
    number = 0
    for key, value in indexesEqualEdges.items():
        #        All occurrences of an "edge" receive the same "edge_number".
        #        for i in value:
        #            edges_new_numbers[i] = number
        edgeNumber_vertexes[number] = col_sorted_e_v[value[0]]
        number += 1

    # construction of edgeNumber_triangles
    print("        construction of edgeNumber_triangles...")
    sys.stdout.flush()
    edgeNumber_triangles = indexesEqualEdges
    for key, value in edgeNumber_triangles.items():
        value_mod_3 = [int(x / 3) for x in value]
        edgeNumber_triangles[key] = value_mod_3

    return edgeNumber_triangles, edgeNumber_vertexes


def compute_triangle_adjacentTriangles(T, edgeNumber_triangles):
    print("    construction of triangle_adjacentTriangles...")
    sys.stdout.flush()
    triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction = {}, {}
    # initialisation of the dictionaries
    for i in range(T):
        triangle_adjacentTriangles[i] = []
    # filling of the dictionaries
    for adjacent_triangles in edgeNumber_triangles.values():
        N_adj_triangles = len(adjacent_triangles)
        IS_JUNCTION = (N_adj_triangles > 2)
        for tn in adjacent_triangles:
            listToAdd = [tj for tj in adjacent_triangles if tj != tn]
            triangle_adjacentTriangles[tn] += listToAdd
            if IS_JUNCTION:
                if tn in is_triangle_adjacentTriangles_via_junction:
                    is_triangle_adjacentTriangles_via_junction[
                        tn] += listToAdd  # [IS_JUNCTION] * (N_adj_triangles-1)
                else:
                    is_triangle_adjacentTriangles_via_junction[tn] = listToAdd
    return triangle_adjacentTriangles, is_triangle_adjacentTriangles_via_junction


def RWGNumber_signedTriangles_computation(edgeNumber_triangles, edgeNumber_vertexes,
                                          triangles_surfaces, is_closed_surface, triangle_vertexes,
                                          vertexes_coord):
    #     we now want to "intelligently" get rid of the junctions,
    #     that means, create new RWGs when needed.
    #     The following is correct for metal-metal junctions
    print("    computation of RWG to triangles relations...")
    sys.stdout.flush()
    N_edges = len(edgeNumber_triangles)
    #    RWGNumber_signedTrianglesTmp_1 is an array of fixed size N_edges.
    #    It will be equal to edgeNumber_triangles if there are no junctions.
    #    If there are junctions, RWGNumber_signedTrianglesTmp_2 will be non-empty.
    #    This is for limiting the memory requirements on this part of the code.
    RWGNumber_signedTrianglesTmp_1, RWGNumber_signedTrianglesTmp_2 = zeros((N_edges, 2), 'i'), {}
    RWGNumber_edgeNumber = range(N_edges)
    for index in range(N_edges):
        numberOfTriangles = len(edgeNumber_triangles[index])
        if numberOfTriangles == 2:
            RWGNumber_signedTrianglesTmp_1[index] = array(edgeNumber_triangles[index], 'i')
        else:  # then we have a junction
            # we compute the edge unit vector: zHatEdge
            n0, n1 = edgeNumber_vertexes[index, 0], edgeNumber_vertexes[index, 1]
            r0, r1 = vertexes_coord[n0], vertexes_coord[n1]
            r1_r0 = r1 - r0
            zHatEdge = r1_r0 / sqrt(dot(r1_r0, r1_r0))
            xHatEdge, yHatEdge = zeros(3, 'd'), zeros(3, 'd')
            triangles = edgeNumber_triangles[index]
            triangles_angles = zeros(numberOfTriangles, 'd')
            # construction of the edge local coordinate system: it is based on the first triangle
            # this is for sorting the triangles wrt their respective angles
            # because we cannot have a RWG which is bisected by another RWG
            for i in range(numberOfTriangles):
                tr = triangles[i]
                tr_nodes = triangle_vertexes[tr]
                for n2 in tr_nodes:
                    if (n2 != n0) and (n2 != n1):
                        break
                r2 = vertexes_coord[n2]
                r2_r0 = r2 - r0
                r2_r0_perpendicular = r2_r0 - dot(r2_r0, zHatEdge) * zHatEdge
                r2_r0_perpendicularHat = r2_r0_perpendicular / sqrt(
                    dot(r2_r0_perpendicular, r2_r0_perpendicular))
                if i == 0:
                    xHatEdge = r2_r0_perpendicularHat[
                               :]  # the first triangle will be the reference: will define the xHat vector
                    yHatEdge[0] = zHatEdge[1] * xHatEdge[2] - zHatEdge[2] * xHatEdge[1]
                    yHatEdge[1] = zHatEdge[2] * xHatEdge[0] - zHatEdge[0] * xHatEdge[2]
                    yHatEdge[2] = zHatEdge[0] * xHatEdge[1] - zHatEdge[1] * xHatEdge[0]
                    triangles_angles[0] = 0.0
                else:
                    Cos = dot(r2_r0_perpendicularHat, xHatEdge)
                    Sin = dot(r2_r0_perpendicularHat, yHatEdge)
                    if Sin >= 0:
                        triangles_angles[i] = arccos(Cos)
                    else:
                        triangles_angles[i] = 2 * np.pi - arccos(Cos)
            # we now sort the triangles wrt their respective position wrt xHatEdge
            ind_sortedTriangles = argsort(triangles_angles, kind='mergesort')
            sortedTriangles = take(triangles, ind_sortedTriangles)
            sortedTriangleSurfaces = triangles_surfaces[sortedTriangles]
            # we now form all the possible RWGs with the sorted triangles. Normally none of these RWGs can be bisected by a triangle now
            possibleTrianglesPairsForRWGs = []
            for i in range(1, numberOfTriangles):
                if (sortedTriangleSurfaces[i] != sortedTriangleSurfaces[i - 1]):
                    possibleTrianglesPairsForRWGs.append(
                        [sortedTriangles[i - 1], sortedTriangles[i]])
            # the last possibility of RWG: between the last and first triangle
            if len(
                    possibleTrianglesPairsForRWGs) < numberOfTriangles - 1:  # if our list is too small
                if (sortedTriangleSurfaces[numberOfTriangles - 1] != sortedTriangleSurfaces[
                    0]):  # if the triangles belong to a different surface
                    possibleTrianglesPairsForRWGs.append(
                        [sortedTriangles[numberOfTriangles - 1], sortedTriangles[0]])
            for i in range(1, len(
                    possibleTrianglesPairsForRWGs)):  # if more than one, we have a junction
                key = len(RWGNumber_signedTrianglesTmp_2)
                RWGNumber_signedTrianglesTmp_2[key] = possibleTrianglesPairsForRWGs[i]
                RWGNumber_edgeNumber.append(index)
            # we don't need the whole list of edges for current line of RWGNumber_signedTriangles
            RWGNumber_signedTrianglesTmp_1[index] = array(possibleTrianglesPairsForRWGs[0], 'i')
    # we can finally construct the final RWGNumber_signedTriangles array
    N_RWG = RWGNumber_signedTrianglesTmp_1.shape[0] + len(RWGNumber_signedTrianglesTmp_2)
    RWGNumber_signedTriangles = zeros((N_RWG, 2), 'i')
    RWGNumber_signedTriangles[:N_edges, :] = RWGNumber_signedTrianglesTmp_1
    index = N_edges
    for key, value in RWGNumber_signedTrianglesTmp_2.items():
        RWGNumber_signedTriangles[index] = array(RWGNumber_signedTrianglesTmp_2[key], 'i')
        index += 1
    if not (index == N_RWG):
        print("Error at the end of RWGNumber_signedTriangles_computation. Exiting")
        sys.exit(1)
        # computation of RWGNumber_edgeVertexes
    # RWGNumber_edgeVertexes = (take(edgeNumber_vertexes, array(RWGNumber_edgeNumber, 'i'), axis=0)).astype('i')
    RWGNumber_edgeVertexes = zeros((N_RWG, 2), 'i')
    for i in range(N_RWG):
        t0, t1 = RWGNumber_signedTriangles[i, 0], RWGNumber_signedTriangles[i, 1]
        s0, s1 = triangles_surfaces[t0], triangles_surfaces[t1]
        t = t0
        if not s0 == s1:
            if is_closed_surface[s0] == 1:
                t = t0
            elif is_closed_surface[s1] == 1:
                t = t1
            else:
                t = t0
        e0, e1 = edgeNumber_vertexes[RWGNumber_edgeNumber[i], 0], \
                 edgeNumber_vertexes[RWGNumber_edgeNumber[i], 1]
        n0, n1, n2 = triangle_vertexes[t, 0], triangle_vertexes[t, 1], triangle_vertexes[t, 2]
        if ((e0 == n0) and (e1 == n1)) or ((e0 == n1) and (e1 == n2)) or (
            (e0 == n2) and (e1 == n0)):
            RWGNumber_edgeVertexes[i, 0] = e0
            RWGNumber_edgeVertexes[i, 1] = e1
        else:
            RWGNumber_edgeVertexes[i, 0] = e1
            RWGNumber_edgeVertexes[i, 1] = e0

    return RWGNumber_signedTriangles.astype('i'), \
           RWGNumber_edgeVertexes.astype('i'), N_edges, N_RWG


def RWGNumber_oppVertexes_computation(RWGNumber_signedTriangles, RWGNumber_edgeVertexes,
                                      triangle_vertexes):
    print("    computation of RWG opposite vertexes...")
    sys.stdout.flush()

    N_RWG = RWGNumber_signedTriangles.shape[0]
    RWGNumber_oppVertexes = zeros((N_RWG, 2), 'i')
    for j in range(N_RWG):
        t0 = RWGNumber_signedTriangles[j, 0]
        t1 = abs(RWGNumber_signedTriangles[j, 1])
        edgeVertexes = RWGNumber_edgeVertexes[j]
        t0Vertexes = triangle_vertexes[t0]
        t1Vertexes = triangle_vertexes[t1]
        for i in t0Vertexes:
            if (i != edgeVertexes[0]) and (i != edgeVertexes[1]):
                break
        RWGNumber_oppVertexes[j, 0] = i
        for i in t1Vertexes:
            if (i != edgeVertexes[0]) and (i != edgeVertexes[1]):
                break
        RWGNumber_oppVertexes[j, 1] = i

    return RWGNumber_oppVertexes.astype('i')


def compute_RWGNumber_edgeCentroidCoord(vertexes_coord, RWGNumber_edgeVertexes):
    edgeCentroidCoord = take(vertexes_coord, RWGNumber_edgeVertexes[:, 0], axis=0)
    edgeCentroidCoord += take(vertexes_coord, RWGNumber_edgeVertexes[:, 1], axis=0)
    return edgeCentroidCoord / 2.0


def compute_RWGNumber_edgeLength(vertexes_coord, RWGNumber_edgeVertexes):
    RWGNumber_length = take(vertexes_coord, RWGNumber_edgeVertexes[:, 0], axis=0)
    RWGNumber_length -= take(vertexes_coord, RWGNumber_edgeVertexes[:, 1], axis=0)
    RWGNumber_length = sqrt(sum((RWGNumber_length) ** 2, axis=1))
    return RWGNumber_length.astype('d')


def compute_RWG_meanEdgeLength(vertexes_coord, RWGNumber_edgeVertexes, stride):
    r0 = take(vertexes_coord, RWGNumber_edgeVertexes[::stride, 0], axis=0)
    r1 = take(vertexes_coord, RWGNumber_edgeVertexes[::stride, 1], axis=0)
    RWGNumber_length = sqrt(sum((r0 - r1) ** 2, axis=1))
    return sum(RWGNumber_length) / RWGNumber_length.shape[0]


def edgeNumber_triangles_indexes(list_of_edges_numbers, RWGNumber_signedTriangles):
    """This function returns a 1-D array of the indexes of the triangles corresponding
    to a 1-D array of edges_numbers. This function is important for creating lists of triangles
    that will participate to the MoM, given a particular criterium concerning the edges."""
    indexes_of_triangles_tmp1 = take(RWGNumber_signedTriangles, list_of_edges_numbers, axis=0).flat
    indexes_of_triangles_tmp2 = sort(indexes_of_triangles_tmp1, kind='mergesort')
    indexes_of_triangles_to_take = ones(indexes_of_triangles_tmp2.shape[0], 'i')
    indexes_of_triangles_to_take[1:] = indexes_of_triangles_tmp2[1:] - indexes_of_triangles_tmp2[
                                                                       :-1]
    indexes_of_triangles = compress(indexes_of_triangles_to_take != 0, indexes_of_triangles_tmp2)
    return indexes_of_triangles.astype('i')


def change_triangle_circulation(t0, t1, triangle_vertexes):
    """This function reorders the nodes of triangle t1 such that the common edge
    between reference triangle t0 and t1 is parcouru following one direction and
    then following the opposite. In this way we will have compatible normals"""
    # nodes of reference triangle
    node00 = triangle_vertexes[t0, 0]
    node01 = triangle_vertexes[t0, 1]
    node02 = triangle_vertexes[t0, 2]
    # nodes of triangle to change
    node10 = triangle_vertexes[t1, 0]
    node11 = triangle_vertexes[t1, 1]
    node12 = triangle_vertexes[t1, 2]
    # circulation on these triangles is a succession of paths between nodes: 0 -> 1, 1 -> 2, 2 -> 0
    t0_circ = [[node00, node01], [node01, node02], [node02, node00]]
    t1_circ = [[node10, node11], [node11, node12], [node12, node10]]
    # "coincide_circ": a variable that will tell if the edge is parcouru alike by each triangle.
    # In this case, one will have to reorder the nodes of t1.
    coincide_circ = 0;
    for k in range(3):
        if coincide_circ == 1:
            break
        for l in range(3):
            if t0_circ[k] == t1_circ[l]:
                coincide_circ = 1
                break
    if coincide_circ == 1:
        triangle_vertexes[t1, :] = array([node10, node12, node11])


def reorder_triangle_vertexes(triangle_adjacentTriangles,
                              is_triangle_adjacentTriangles_via_junction,
                              triangle_vertexes, vertexes_coord):
    """This function is necessary to orientate all the normals of the triangles
    coherently for a given surface. For this purpose it gives a good order of
    appearance of the triangles vertexes. This function acts directly on
    the array 'triangle_vertexes'.
    
    This function also returns triangles_surfaces"""

    print("      reordering triangles for normals coherency...")
    sys.stdout.flush()
    T = len(triangle_adjacentTriangles)
    is_triangle_reordered = [
                                0] * T  # create a list of length "T" that tells if the triangle has been reordered
    triangles_surfaces = array([-1] * T, 'i')
    surf_number = -1
    while 0 in is_triangle_reordered:  # as long as there is a triangle that has not been reordered
        surf_number += 1
        t_start = is_triangle_reordered.index(0)  # the triangle from which we start the growing
        is_triangle_reordered[t_start] = 1
        triangles_surfaces[t_start] = surf_number
        # index = 0
        list_t_to_reorder = []
        for tn in triangle_adjacentTriangles[t_start]:
            is_triangle_not_adjacent_via_junction = True
            if t_start in is_triangle_adjacentTriangles_via_junction:
                if tn in is_triangle_adjacentTriangles_via_junction[t_start]:
                    is_triangle_not_adjacent_via_junction = False
            if (is_triangle_reordered[int(tn)] == 0) and (is_triangle_not_adjacent_via_junction):
                list_t_to_reorder.append(tn)
                # index += 1
        # list_t_to_reorder = [tn for tn in triangle_adjacentTriangles[t_start] if is_triangle_reordered[int(tn)]==0] # initialization of the list to process
        list_calling_t = [t_start] * len(list_t_to_reorder)
        while list_t_to_reorder:
            t = list_t_to_reorder.pop()
            calling_t = list_calling_t.pop()
            # we change circulation of t according to calling_t
            change_triangle_circulation(calling_t, t, triangle_vertexes)
            is_triangle_reordered[int(t)] = 1
            triangles_surfaces[int(t)] = surf_number
            # index = 0
            t_adjacent_triangles = []
            for tn in triangle_adjacentTriangles[t]:
                is_triangle_not_adjacent_via_junction = True
                if t in is_triangle_adjacentTriangles_via_junction:
                    if tn in is_triangle_adjacentTriangles_via_junction[t]:
                        is_triangle_not_adjacent_via_junction = False
                if (is_triangle_reordered[int(tn)] == 0) and (
                is_triangle_not_adjacent_via_junction):
                    t_adjacent_triangles.append(tn)
                    # index += 1
            # we extend the list to reorder with the triangles adjacent to t
            list_t_to_reorder.extend(t_adjacent_triangles)
            # we extend the list of "calling" triangles with t
            list_calling_t.extend([int(t)] * len(t_adjacent_triangles))
    # we loop on the surfaces, because all normals are coherent but maybe not directed outwards closed surfaces...
    print("      redirecting the normals outward...")
    sys.stdout.flush()
    S = surf_number
    triangles_centroids_z = triangles_centroids_computation(vertexes_coord, triangle_vertexes)[:, 2]
    for s in range(S + 1):
        ind_t_on_s = nonzero(triangles_surfaces == s)[0]
        #        max_height_centroids_s = np.max(take(triangles_centroids_z, ind_t_on_s, axis=0))
        ind_max_height_centroids_s_tmp \
            = np.argmax(take(triangles_centroids_z, ind_t_on_s, axis=0))
        #        ind_max_height_centroids_s_tmp2 = triangles_centroids_z.tolist().index(max_height_centroids_s)

        ind_max_height_centroids_s = ind_t_on_s[ind_max_height_centroids_s_tmp]
        r0 = vertexes_coord[triangle_vertexes[ind_max_height_centroids_s, 0]]
        r1 = vertexes_coord[triangle_vertexes[ind_max_height_centroids_s, 1]]
        r2 = vertexes_coord[triangle_vertexes[ind_max_height_centroids_s, 2]]
        r1_r0 = r1 - r0
        r2_r0 = r2 - r0
        triangle_normal = zeros(3, 'd')
        triangle_normal[0] = r1_r0[1] * r2_r0[2] - r1_r0[2] * r2_r0[1]
        triangle_normal[1] = r1_r0[2] * r2_r0[0] - r1_r0[0] * r2_r0[2]
        triangle_normal[2] = r1_r0[0] * r2_r0[1] - r1_r0[1] * r2_r0[0]
        # we now test the normal associated to ind_max_height_centroids_s
        z_hat = array([0, 0, 1], 'd')
        if sum(triangle_normal * z_hat) < 0:
            # we also swap the columns of triangle_vertexes
            for index in ind_t_on_s:  # this part does not work with put...
                vertexes = copy.copy(triangle_vertexes[index, :])
                triangle_vertexes[index, 1] = vertexes[2]
                triangle_vertexes[index, 2] = vertexes[1]
    return triangles_surfaces.astype('i')


def is_surface_closed(triangles_surfaces, edgeNumber_triangles):
    """this function determines if a surface is open or closed.
       it also provides relationships between surfaces (linked or not)"""
    S = max(triangles_surfaces) + 1
    NUMBER_TRIANGLES_IN_SURFACE = zeros(S, 'i')
    for s in triangles_surfaces:
        NUMBER_TRIANGLES_IN_SURFACE[s] += 1
    connected_surfaces = {}
    # we now count the number of INNER edges for each surface.
    # the edges that are junctions receive a special treatment:
    # only if the edge has two triangles on the given surface, 
    # can it be counted as an inner edge, which will then be 
    # counted in NUMBER_EDGES_IN_SURFACE
    NUMBER_EDGES_IN_SURFACE = zeros(S, 'i')
    for edge_number, triangles_tmp in edgeNumber_triangles.items():
        surfaces_appeared_already = zeros(S, 'i')
        if len(triangles_tmp) > 2:  # we have a junction here
            for t in triangles_tmp:
                surface = triangles_surfaces[t]
                if surfaces_appeared_already[surface] == 0:
                    surfaces_appeared_already[surface] = 1
                else:
                    NUMBER_EDGES_IN_SURFACE[surface] += 1
            surfaces_present = compress(surfaces_appeared_already > 0, arange(S))
            if len(surfaces_present) == 2:
                s0, s1 = min(surfaces_present), max(surfaces_present)
                if (s0, s1) in connected_surfaces:
                    connected_surfaces[(s0, s1)].append(edge_number)
                else:
                    connected_surfaces[(s1, s0)] = [edge_number]
            else:
                for index1 in arange(len(surfaces_present)):
                    for index2 in arange(index1 + 1, len(surfaces_present)):
                        s1 = min(surfaces_present[index1], surfaces_present[index2])
                        s2 = max(surfaces_present[index1], surfaces_present[index2])
                        if (s1, s2) in connected_surfaces:
                            connected_surfaces[(s1, s2)].append(edge_number)
                        else:
                            connected_surfaces[(s1, s2)] = [edge_number]
        else:
            surface = triangles_surfaces[triangles_tmp[0]]
            NUMBER_EDGES_IN_SURFACE[surface] += 1

    is_closed_surface = ((NUMBER_EDGES_IN_SURFACE * 2) == (NUMBER_TRIANGLES_IN_SURFACE * 3))
    # we now check for potential closed surfaces: surfaces which can be closed
    # and on which we can therefore apply the CFIE
    potential_closed_surfaces = {}
    for key, item in connected_surfaces.items():
        s0, s1 = key[0], key[1]
        numberEdges0, numberEdges1 = NUMBER_EDGES_IN_SURFACE[s0], NUMBER_EDGES_IN_SURFACE[s1]
        numberTriangles0, numberTriangles1 = NUMBER_TRIANGLES_IN_SURFACE[s0], \
                                             NUMBER_TRIANGLES_IN_SURFACE[s1]
        if (numberEdges0 + numberEdges1 + len(item)) * 2 == 3 * (
            numberTriangles0 + numberTriangles1):
            potential_closed_surfaces[key] = item
    return is_closed_surface * 1, connected_surfaces, potential_closed_surfaces


def triangles_unnormalized_normals_computation(vertexes_coord, triangle_vertexes, t):
    """This function returns the non-normalized normals of each triangle.
    t is an array of the t indexes to be considered"""
    triangles_normals = zeros((t.shape[0], 3), 'd')
    stride = 10000
    startIndex, stopIndex = 0, min(stride, t.shape[0])
    # it is coded this way for memory optimization: this function is really a memory hog!
    while startIndex < triangles_normals.shape[0]:
        indexes = t[startIndex:stopIndex]
        v0 = take(triangle_vertexes[:, 0], indexes, axis=0)
        v1 = take(triangle_vertexes[:, 1], indexes, axis=0)
        v2 = take(triangle_vertexes[:, 2], indexes, axis=0)
        r0 = take(vertexes_coord, v0, axis=0)  # first vertexes of all triangles
        r1_r0 = take(vertexes_coord, v1, axis=0) - r0
        r2_r0 = take(vertexes_coord, v2, axis=0) - r0
        triangles_normals[indexes, 0] = r1_r0[:, 1] * r2_r0[:, 2] - r1_r0[:, 2] * r2_r0[:, 1]
        triangles_normals[indexes, 1] = r1_r0[:, 2] * r2_r0[:, 0] - r1_r0[:, 0] * r2_r0[:, 2]
        triangles_normals[indexes, 2] = r1_r0[:, 0] * r2_r0[:, 1] - r1_r0[:, 1] * r2_r0[:, 0]
        startIndex = stopIndex
        stopIndex = min(stopIndex + stride, t.shape[0])
    return triangles_normals.astype('d')


def triangles_areas_normals_computation(vertexes_coord, triangle_vertexes, triangles_surfaces):
    """this function"""
    T = triangle_vertexes.shape[0]
    #    S = max(triangles_surfaces)
    triangles_normals = triangles_unnormalized_normals_computation(vertexes_coord,
                                                                   triangle_vertexes,
                                                                   arange(T).astype('i')).astype(
        'd')
    norm_triangles_normals = reshape(sqrt(sum(triangles_normals ** 2, 1)), (T, 1))
    triangles_areas = norm_triangles_normals / 2.0
    triangles_normals /= norm_triangles_normals
    return triangles_areas.astype('d'), triangles_normals.astype('d')


def write_normals(name, triangles_centroids, triangles_normals, triangles_surfaces, surface):
    """function that writes the normals of a given surface to a file readable by GMSH for viewing."""
    f = open(name, 'w')
    f.write('View "normals of surfaces xxx" {\n')
    T = triangles_surfaces.shape[0]
    for k in range(T):
        write_condition = (triangles_surfaces[k] == surface) or (
        surface == -1)  # if surface == -1, we write all the normals
        if write_condition:
            string_to_write = 'VP(' + str(triangles_centroids[k, :].tolist())[1:-1] + ')'
            string_to_write += '{' + str(triangles_normals[k, :].tolist())[1:-1] + '};\n'
            f.write(string_to_write)
    f.write('};\n')
    f.close()


# Functions for computing the barycentric mesh and the P (M for me) matrix from barycentric to original mesh
def divide_triangles(RWGNumber_signedTriangles, RWGNumber_edgeVertexes, reordered_triangle_vertexes,
                     vertexes_coord):
    N_RWG = RWGNumber_signedTriangles.shape[0]
    V = vertexes_coord.shape[0]
    T = reordered_triangle_vertexes.shape[0]
    divided_triangles_vertexes = zeros((T, 7), 'int32')
    # in divided_triangles_vertexes, the column indexes:
    # 0, 2, 4 correspond to the nodes 0, 1, 2 of the triangles_vertexes
    # 1, 3, 5 correspond to the midpoints of the edges n01, n12, n20 of the triangles
    # 6 corresponds to the r_grav of the triangle 
    divided_triangles_vertexes[:, 0] = reordered_triangle_vertexes[:, 0]
    divided_triangles_vertexes[:, 2] = reordered_triangle_vertexes[:, 1]
    divided_triangles_vertexes[:, 4] = reordered_triangle_vertexes[:, 2]
    divided_triangles_vertexes[:, 6] = arange(V, V + T)
    # first we take care of the points located on RWG edges
    number_of_edge_centroid = V + T
    for i in range(N_RWG):
        e0 = RWGNumber_edgeVertexes[i, 0]
        e1 = RWGNumber_edgeVertexes[i, 1]
        edge = [e0, e1]
        edge.sort()
        index_in_t = [0, 0]
        edge_centroid_number = [0, 0]
        for j in range(2):
            t = RWGNumber_signedTriangles[i, j]
            n0 = reordered_triangle_vertexes[t, 0]
            n1 = reordered_triangle_vertexes[t, 1]
            n2 = reordered_triangle_vertexes[t, 2]
            edges_triangle = [[n0, n1], [n1, n2], [n2, n0]]
            for e in edges_triangle:
                e.sort()
            index_in_t[j] = edges_triangle.index(edge) * 2 + 1
            edge_centroid_number[j] = divided_triangles_vertexes[t, index_in_t[j]]
        # we now check if one of the edges has already a number. If not we assign a number.
        if edge_centroid_number[0] != 0:
            number = edge_centroid_number[0]
        elif edge_centroid_number[1] != 0:
            number = edge_centroid_number[1]
        else:
            number = number_of_edge_centroid
            number_of_edge_centroid += 1
        for j in range(2):
            t = RWGNumber_signedTriangles[i, j]
            divided_triangles_vertexes[t, index_in_t[j]] = number

    # then we take care of the points located on triangle edges that are not RWGs
    for t in range(T):
        for j in range(3):
            if (divided_triangles_vertexes[t, j * 2 + 1] == 0):
                divided_triangles_vertexes[t, j * 2 + 1] = number_of_edge_centroid
                number_of_edge_centroid += 1

    # that's all for the barycentric division of the triangles.
    return divided_triangles_vertexes, number_of_edge_centroid

def compute_triangle_arrays(RWGNumber_signedTriangles, T):
    N_RWG = RWGNumber_signedTriangles.shape[0]

    triangle_RWGNumber = -np.ones([T, 3], dtype="i")
    triangle_signInRWG = np.zeros([T, 3], dtype="i")
    triangle_NRWG = np.zeros(T, dtype="i")
    triangleNumbers = RWGNumber_signedTriangles.reshape(2 * N_RWG)
    signInRWG = np.c_[np.ones(N_RWG, dtype="i"), -np.ones(N_RWG, dtype="i")]
    signInRWG = signInRWG.reshape(2 * N_RWG)

    ind = np.argsort(triangleNumbers)
    RWGNumbers = ind/2
    triangleNumbers = triangleNumbers[ind]
    signInRWG = signInRWG[ind]

    tr_number, j = -1, 0
    for i in range(N_RWG*2):

        if (tr_number != triangleNumbers[i]):
            j, tr_number = 0, triangleNumbers[i]
        else:
            j += 1

        triangle_RWGNumber[tr_number, j] = RWGNumbers[i]
        triangle_signInRWG[tr_number, j] = signInRWG[i]
        triangle_NRWG[tr_number] += 1

    return triangle_NRWG, triangle_RWGNumber, triangle_signInRWG


def run(simuDirName):
    mesh_data_path = os.path.join(simuDirName, "mesh")

    T = readIntFromDisk(mesh_data_path + "/T.txt")
    V = readIntFromDisk(mesh_data_path + "/V.txt")

    triangle_vertexes \
        = readBlitzArrayFromDisk(mesh_data_path + "/triangle_vertexes.txt",
                                 T, 3, "i")
    vertexes_coord \
        = readBlitzArrayFromDisk(mesh_data_path + "/vertexes_coord.txt",
                                 V, 3, "d")

    edgeNumber_vertexes, edgeNumber_triangles, triangle_adjacentTriangles, \
    is_triangle_adjacentTriangles_via_junction \
        = edges_computation(triangle_vertexes, vertexes_coord)

    triangles_surfaces \
        = reorder_triangle_vertexes(triangle_adjacentTriangles,
                                    is_triangle_adjacentTriangles_via_junction,
                                    triangle_vertexes, vertexes_coord)

    print("    checking open and closed surfaces...")
    is_closed_surface, connected_surfaces, potential_closed_surfaces \
        = is_surface_closed(triangles_surfaces, edgeNumber_triangles)
    print("    is_closed_surface = " + str(is_closed_surface * 1))
    print("    connected surfaces = " + str(connected_surfaces))
    print("    potential closed surfaces = " + str(potential_closed_surfaces))

    # if 0 in is_closed_surface:
    #     print "error: some surfaces is not closed"
    #     sys.exit()

    RWGNumber_signedTriangles, RWGNumber_edgeVertexes, N_edges, N_RWG \
        = RWGNumber_signedTriangles_computation(edgeNumber_triangles, edgeNumber_vertexes,
                                                triangles_surfaces, is_closed_surface,
                                                triangle_vertexes, vertexes_coord)
    RWGNumber_surfaces = triangles_surfaces.take(RWGNumber_signedTriangles[:, 0])

    # print RWGNumber_signedTriangles
    RWGNumber_oppVertexes \
        = RWGNumber_oppVertexes_computation(RWGNumber_signedTriangles,
                                            RWGNumber_edgeVertexes, triangle_vertexes)
    print("    Number of edges = " + str(N_edges))
    print("    Number of RWG = " + str(N_RWG))

    triangle_NRWG, triangle_RWGNumber, triangle_signInRWG \
        = compute_triangle_arrays(RWGNumber_signedTriangles, T)

    node0, node1 = RWGNumber_edgeVertexes[:, 0], RWGNumber_edgeVertexes[:, 1]
    RWG_edgeVector = vertexes_coord[node0, :] - vertexes_coord[node1, :]
    RWG_length = np.linalg.norm(RWG_edgeVector, axis=1)
    mean_RWG_length = np.sum(RWG_length) / N_RWG
    S = is_closed_surface.size
    if S != np.max(RWGNumber_surfaces) + 1:
        print("error: S is not coincidence...")

    writeScalarToDisk(S, mesh_data_path + "/S.txt")
    writeScalarToDisk(N_edges, mesh_data_path + "/N_edges.txt")
    writeScalarToDisk(N_RWG, mesh_data_path + "/N_RWG.txt")
    writeScalarToDisk(mean_RWG_length, mesh_data_path + "/average_RWG_length.txt")
    writeBlitzArrayToDisk(triangle_vertexes,
                          mesh_data_path + "/triangle_vertexes.txt")
    writeBlitzArrayToDisk(RWGNumber_signedTriangles,
                          mesh_data_path + "/RWGNumber_signedTriangles.txt")
    writeBlitzArrayToDisk(RWGNumber_edgeVertexes,
                          mesh_data_path + "/RWGNumber_edgeVertexes.txt")
    writeBlitzArrayToDisk(RWGNumber_oppVertexes,
                          mesh_data_path + "/RWGNumber_oppVertexes.txt")
    writeBlitzArrayToDisk(triangles_surfaces,
                          mesh_data_path + "/triangles_surfaces.txt")
    writeBlitzArrayToDisk(RWGNumber_surfaces,
                          mesh_data_path + "/RWGNumber_surfaces.txt")
    writeBlitzArrayToDisk(is_closed_surface.astype("i"),
                          mesh_data_path + "/is_closed_surface.txt")

    np.save(mesh_data_path + "/triangle_NRWG.npy", triangle_NRWG)
    np.save(mesh_data_path + "/triangle_RWGNumber", triangle_RWGNumber)
    np.save(mesh_data_path + "/triangle_signInRWG", triangle_signInRWG)

if __name__ == "__main__":
    run("./")
