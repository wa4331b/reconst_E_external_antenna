/**********************************************************************
 *
 * octtree.h
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

#ifndef OCTTREE_H
#define OCTTREE_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <algorithm>

using namespace std;

//#include "mesh.hpp"
#include "level.hpp"

class Octtree
{
	int L;
	int numberOfUpdates;
public:
	double w;
	int VERBOSE;
	int N_GaussOnTriangle;
	int N_levels, ALLOW_CEILING_LEVEL;
	complex<double> k, k_d;
	complex<double> eps_r, mu_r;
	vector<Level> levels;
	double big_cube_lower_coord[3];
	double big_cube_center_coord[3];
	ArrayXcd CFIE;
	ArrayXi target_surface;
	string octtreeDataPath;

	// constructors
	Octtree(void)
	{
	}

	Octtree(const string octtree_data_path,
	        const string media_path,
	        int C, int C_obs,
	        double cubes_centroids_tmp[], double cubes_centroids_obs_tmp[],
	        int N_target_surface, int target_surface_tmp[],
	        int octtreeNthetas_tmp[], double octtreeXthetas_tmp[], double octtreeWthetas_tmp[],
	        int octtreeNphis_tmp[], double octtreeXphis_tmp[], double octtreeWphis_tmp[],
	        int LExpansion_tmp[], complex<double> CFIEcoeffs_tmp[],
	        double bigCubeLowerCoord_tmp[], double bigCubeCenterCoord_tmp[]);

	Octtree(const string octtree_data_path,
	        const string media_path,
	        const ArrayXXd& cubes_centroids,
	        const ArrayXXd& cubes_centroids_obs,
	        const ArrayXi& target_surface,
	        ArrayXi& octtreeNthetas, ArrayXXd& octtreeXthetas, ArrayXXd& octtreeWthetas,
	        ArrayXi& octtreeNphis, ArrayXXd& octtreeXphis, ArrayXXd& octtreeWphis,
	        ArrayXi& LExpansion, ArrayXcd& CFIEcoeffs,
	        ArrayXd& bigCubeLowerCoord, ArrayXd& bigCubeCenterCoord);

	void copyOcttree(const Octtree&);
	Octtree(const Octtree&); // copy constructor
	Octtree& operator=(const Octtree&); // copy assignment operator
	void constructArrays(void);
	void computeGaussLocatedArguments(
		const ArrayXi& local_cubes_NRWG,
	 const ArrayXi& local_RWG_numbers,
	 const ArrayXi& local_RWG_numbers_surf,
	 const ArrayXXi& local_RWGNumbers_signedTriangles,
	 const ArrayXXd& local_RWGNumbers_trianglesCoord
	 	);
	void computeObsLocatedArguments(const ArrayXi& cubes_Nobs,
	                                const ArrayXi& obs_numbers,
	                                const ArrayXXd& obsNumber_obsCoord);
	void RWGs_renumbering(void);
	void computeIndexesOfCubesInOriginalMesh(ArrayXi& oldIndexesOfCubes,
	                                         ArrayXi& oldIndexesOfCubes_obs);
	~Octtree();

	// // functions
	// int getNumberOfUpdates (void) const {return numberOfUpdates;}
	// void setNumberOfUpdates (const int n) {numberOfUpdates = n;}
	// int getL(void) const {return L;}
	// int getLevelsSize (void) const {return levels.size();}
	// complex<double> getK(void) const {return k;}
	// complex<float> getEps_r(void) const {return eps_r;}
	// complex<float> getMu_r(void) const {return mu_r;}
	// float getW(void) const {return w;}
	Level getLevel(const int l) const { return levels[l]; }

	Cube getCubeLevel(const int i, const int l)
	const { return levels[l].getCube(i); }

	// void assignCubesToProcessors(const int /*CUBES_DISTRIBUTION*/);
	// void writeAssignedLeafCubesToDisk(const string /*path*/,
	//                                   const string /*filename*/);
	void updateSup(const ArrayXcd&,
	               const bool);
	void updateSup_obs(const ArrayXXcd&,
	                   const bool, const bool);
	void alphaTranslations(const bool);
	void alphaTranslations_obs();
	void alphaTranslations_obs_transpose();
	// void exchangeSupsIndividually
	// (blitz::Array< blitz::Array<complex<float>, 2>, 1>& /*SupThisLevel*/,
	//  const int /*l*/, const vector<int> & /*localCubesIndexes*/);
	// void exchangeSupsInBlocks
	// (blitz::Array< blitz::Array<complex<float>, 2>, 1>& /*SupThisLevel*/,
	//  const int /*l*/, const vector<int> & /*localCubesIndexes*/);
	void ZIFarComputation(ArrayXcd& ZI,
	                      const ArrayXcd& I_PQ,
	                      const bool, const bool);
	void ZINearComputation(ArrayXXcd& E_obs,
	                       const ArrayXcd& I_PQ,
	                       const bool M_current);
	void ZINearComputation_transpose(ArrayXXcd& E_obs,
	                                 ArrayXcd& I_PQ,
	                                 const bool M_current);
	// blitz::Array<complex<float>, 1> getCFIE(void) const {return CFIE;}
	vector<int> getNeighborsSonsIndexes(const int, const int) const;
	vector<int> getNeighborsSonsIndexes_obs(const int, const int) const;
	void findAlphaTransParticipantsIndexes(const int l);
	void SupAlphaMultiplication
	(ArrayXXcd& SupAlpha,
	 const ArrayXXcd& Sup,
	 const ArrayXcd& alphaTranslation,
	 const ArrayXi& alphaTranslationIndexesNonZeros,
	 const ArrayXi& alphaTranslationIndexes,
	 const int alphaCartesianCoord[3]);
	// void SupAlphaMultiplicationDirections
	// (blitz::Array<complex<float>, 2>& /*SupAlpha*/,
	//  const blitz::Array<complex<float>, 2>& /*Sup*/,
	//  const blitz::Array<complex<float>, 1>& /*alphaTranslation*/,
	//  const blitz::Array<int, 1>& /*alphaTranslationIndexesNonZeros*/,
	//  const int alphaCartesianCoord[3]);
	void alphaTranslationsToCube
	(ArrayXXcd& /*S_tmp*/,
	 const Array<ArrayXXcd, Dynamic, 1>& /*LevelSup*/,
	 const int /*l*/,
	 const int /*cubeIndex*/,
	 const vector<int>&, /*indexesAlphaParticipants*/
	 const bool, const bool);
	void alphaTranslationsToCube_obs
	(ArrayXXcd& S_tmp,
	 const Array<ArrayXXcd, Dynamic, 1>& LevelSup,
	 const int l,
	 const int cubeIndex,
	 const vector<int>& indexesAlphaParticipants);
	void shiftExp(ArrayXXcd& S,
	              const ArrayXcd& shiftingArray);
	void computeFarField(complex<double> e_theta_far_tmp[],
	                     complex<double> e_phi_far_tmp[],
	                     double r_phase_center_tmp[],
	                     double octtreeXthetas_coarsest_tmp[],
	                     double octtreeXphis_coarsest_tmp[],
	                     int octtreeNthetas_coarsest,
	                     int octtreeNphis_coarsest,
	                     complex<double> I_J_tmp[],
	                     complex<double> I_M_tmp[],
	                     int N_RWG,
	                     string octtree_data_path);
	void computeSupLastLevel(ArrayXXcd& SupLastLevel,
	                         const ArrayXd& r_phase_center,
	                         const ArrayXd& octtreeXthetas_coarsest,
	                         const ArrayXd& octtreeXphis_coarsest,
	                         const ArrayXcd& I_PQ,
	                         const string octtree_data_path);
	// void computeSourceFarField
	// (blitz::Array<complex<float>, 2>& e_theta_far,
	//  blitz::Array<complex<float>, 2>& e_phi_far,
	//  const blitz::Array<float, 1>& r_phase_center,
	//  const blitz::Array<float, 1>& octtreeXthetas_coarsest,
	//  const blitz::Array<float, 1>& octtreeXphis_coarsest,
	//  const int J_DIPOLES_EXCITATION,
	//  const blitz::Array<complex<float>, 2>& J_dip,
	//  const blitz::Array<float, 2>& r_J_dip,
	//  const int M_DIPOLES_EXCITATION,
	//  const blitz::Array<complex<float>, 2>& M_dip,
	//  const blitz::Array<float, 2>& r_M_dip);
	// void computeDipoleSup(blitz::Array<complex<float>, 2> & Sup,
	//                       const complex<float> J_dipole[3],
	//                       const int IS_J_CURRENT,
	//                       const float r_dipole[3],
	//                       const float rCenter[3],
	//                       const blitz::Array<float, 1>& thetas,
	//                       const blitz::Array<float, 1>& phis);
	// void resizeSdownLevelsToZero(void) {for (unsigned int i = 0 ; i < levels.size() ; ++i) levels[i].Sdown.resize(0);}
};
#endif
