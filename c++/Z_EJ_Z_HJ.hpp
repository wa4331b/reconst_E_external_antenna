/**********************************************************************
 *
 * Z_EJ_Z_HJ.h
 *
 * Copyright (C) 2014 Idesbald Van den Bosch
 *
 * This file is part of Puma-EM.
 *
 * Puma-EM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Puma-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Puma-EM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Suggestions/bugs : <vandenbosch.idesbald@gmail.com>
 *
 **********************************************************************/

#ifndef Z_EJ_Z_HJ_H
#define Z_EJ_Z_HJ_H

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

void wrap_Z_CFIE_computation(complex<double> Z_CFIE_EJ_tmp[],
                             complex<double> Z_CFIE_EM_tmp[],
                             complex<double> Z_CFIE_HJ_tmp[],
                             complex<double> Z_CFIE_HM_tmp[],
                             complex<double> CFIE[],
                             int CFIE_bool[],
                             int M_current_tmp,
                             double signSurfObs,
                             double signSurfSrc,
                             int N_RWG, int N_RWG_src, int N_RWG_test, int V,
                             int numbers_RWG_src_tmp[],
                             int numbers_RWG_test_tmp[],
                             int srcRWGNumber_surface[],
                             int RWGNumber_signedTriangles_tmp[],
                             int RWGNumber_nodes_tmp[],
                             double nodesCoord_tmp[],
                             double w,
                             complex<double> eps_r,
                             complex<double> mu_r,
                             int N_target_surface, int target_surface_tmp[]);

void wrap_Z_CFIE_test_computation(complex<double> Z_CFIE_tmp[],
                                  int N_RWG, int V,
                                  int srcRWGNumber_surface[],
                                  int RWGNumber_signedTriangles_tmp[],
                                  int RWGNumber_nodes_tmp[],
                                  double nodesCoord_tmp[],
                                  double w,
                                  complex<double> eps_r,
                                  complex<double> mu_r,
                                  int N_target_surface, int target_surface_tmp[]);

#endif
