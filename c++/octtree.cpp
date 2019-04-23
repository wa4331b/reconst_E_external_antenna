/**********************************************************************
 *
 * octtree.cpp
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

#include <fstream>
#include <iostream>
#include <string>
#include <complex>
#include <cmath>
#include <vector>

using namespace std;

#include "octtree.hpp"
#include "readWriteBlitzArrayFromFile.hpp"
#include "interpolation.hpp"

/****************************************************************************/
/********************************** Octtree *********************************/
/****************************************************************************/
Octtree::Octtree(const string octtree_data_path,
                 const string media_path,
                 int C, int C_obs,
                 double cubes_centroids_tmp[], double cubes_centroids_obs_tmp[],
                 int N_target_surface, int target_surface_tmp[],
                 int octtreeNthetas_tmp[], double octtreeXthetas_tmp[], double octtreeWthetas_tmp[],
                 int octtreeNphis_tmp[], double octtreeXphis_tmp[], double octtreeWphis_tmp[],
                 int LExpansion_tmp[], complex<double> CFIEcoeffs_tmp[],
                 double bigCubeLowerCoord_tmp[], double bigCubeCenterCoord_tmp[])
{
	typedef Array<double, Dynamic, Dynamic, RowMajor> ArrayXXd_RowMajor;

	ArrayXXd cubes_centroids = Map<ArrayXXd_RowMajor>(cubes_centroids_tmp, C, 3);
	ArrayXXd cubes_centroids_obs = Map<ArrayXXd_RowMajor>(cubes_centroids_obs_tmp, C_obs, 3);
	ArrayXi target_surface = Map<ArrayXi>(target_surface_tmp, N_target_surface);

	readIntFromASCIIFile(octtree_data_path + "N_active_levels.txt", N_levels);

	ArrayXi octtreeNthetas = Map<ArrayXi>(octtreeNthetas_tmp, N_levels);
	ArrayXi octtreeNphis = Map<ArrayXi>(octtreeNphis_tmp, N_levels);

	int N_Xthetas = octtreeNthetas(N_levels - 1), N_Xphis = octtreeNphis(N_levels - 1);
	ArrayXXd octtreeXthetas = Map<ArrayXXd_RowMajor>(octtreeXthetas_tmp, N_levels, N_Xthetas);
	ArrayXXd octtreeWthetas = Map<ArrayXXd_RowMajor>(octtreeWthetas_tmp, N_levels, N_Xthetas);
	ArrayXXd octtreeXphis = Map<ArrayXXd_RowMajor>(octtreeXphis_tmp, N_levels, N_Xphis);
	ArrayXXd octtreeWphis = Map<ArrayXXd_RowMajor>(octtreeWphis_tmp, N_levels, N_Xphis);

	ArrayXi LExpansion = Map<ArrayXi>(LExpansion_tmp, N_levels);
	ArrayXcd CFIEcoeffs = Map<ArrayXcd>(CFIEcoeffs_tmp, 4);
	ArrayXd bigCubeLowerCoord = Map<ArrayXd>(bigCubeLowerCoord_tmp, 3);
	ArrayXd bigCubeCenterCoord = Map<ArrayXd>(bigCubeCenterCoord_tmp, 3);

	this->target_surface = target_surface;

	octtreeDataPath = octtree_data_path;
	// verbose or not?
	readIntFromASCIIFile(octtree_data_path + "VERBOSE.txt", VERBOSE);

	if (VERBOSE == 1) cout << "creating the tree " << " from disk data" << endl;
	numberOfUpdates = 0;
	readIntFromASCIIFile(octtree_data_path + "N_active_levels.txt", N_levels);
	if (VERBOSE == 1) cout << "N active levels = " << N_levels << endl;
	readIntFromASCIIFile
		(octtree_data_path + "ALLOW_CEILING_LEVEL.txt", ALLOW_CEILING_LEVEL);
	if (VERBOSE == 1)
		cout << "ALLOW_CEILING_LEVEL = " << ALLOW_CEILING_LEVEL << endl;
	readIntFromASCIIFile
		(octtree_data_path + "N_GaussOnTriangle.txt", N_GaussOnTriangle);
	if (VERBOSE == 1) cout << "N_GaussOnTriangle = " << N_GaussOnTriangle << endl;

	// theta data
	int INCLUDED_THETA_BOUNDARIES, PERIODIC_Theta, CYCLIC_Theta, NOrderInterpTheta;
	readIntFromASCIIFile(octtree_data_path + "INCLUDED_THETA_BOUNDARIES.txt",
	                     INCLUDED_THETA_BOUNDARIES);
	readIntFromASCIIFile(octtree_data_path + "PERIODIC_Theta.txt", PERIODIC_Theta);
	readIntFromASCIIFile(octtree_data_path + "CYCLIC_Theta.txt", CYCLIC_Theta);
	readIntFromASCIIFile(octtree_data_path + "NOrderInterpTheta.txt",
	                     NOrderInterpTheta);

	// phi data
	int INCLUDED_PHI_BOUNDARIES, PERIODIC_Phi, CYCLIC_Phi, NOrderInterpPhi;
	readIntFromASCIIFile(octtree_data_path + "INCLUDED_PHI_BOUNDARIES.txt",
	                     INCLUDED_PHI_BOUNDARIES);
	readIntFromASCIIFile(octtree_data_path + "PERIODIC_Phi.txt", PERIODIC_Phi);
	readIntFromASCIIFile(octtree_data_path + "CYCLIC_Phi.txt", CYCLIC_Phi);
	readIntFromASCIIFile(octtree_data_path + "NOrderInterpPhi.txt",
	                     NOrderInterpPhi);

	double leaf_side_length;
	readDoubleFromASCIIFile
		(octtree_data_path + "leaf_side_length.txt", leaf_side_length);

	if (VERBOSE == 1) cout << "L expansions = " << LExpansion << endl;
	if (VERBOSE == 1) cout << "CFIE coeffs = " << CFIEcoeffs << endl;

	for (int i = 0; i < 3; ++i)
	{
		big_cube_lower_coord[i] = bigCubeLowerCoord(i);
		big_cube_center_coord[i] = bigCubeCenterCoord(i);
	}

	readDoubleFromASCIIFile(octtree_data_path + "w.txt", this->w);
	readComplexDoubleFromASCIIFile(octtree_data_path + media_path + "k.txt", this->k);
	readComplexDoubleFromASCIIFile(octtree_data_path + media_path + "eps_r.txt", this->eps_r);
	readComplexDoubleFromASCIIFile(octtree_data_path + media_path + "mu_r.txt", this->mu_r);

	cout << "w = " << w << ", mu_rel = " << mu_r
		<< ", eps_rel = " << eps_r << endl;

	CFIE.resize(CFIEcoeffs.size());
	CFIE = CFIEcoeffs;
	levels.reserve(N_levels); // reserve the memory
	double minCubesArraysSize;
	int indMinCubesArraysSize;
	{
		// construction of the first level
		Array<double, 1, Dynamic> XthetasThisLevel, WthetasThisLevel, XphisThisLevel, WphisThisLevel;
		XthetasThisLevel = octtreeXthetas.block(0, 0, 1, octtreeNthetas(0));
		WthetasThisLevel = octtreeWthetas.block(0, 0, 1, octtreeNthetas(0));
		XphisThisLevel = octtreeXphis.block(0, 0, 1, octtreeNphis(0));
		WphisThisLevel = octtreeWphis.block(0, 0, 1, octtreeNphis(0));

		Array<double, 1, Dynamic> XthetasNextLevel, XphisNextLevel;
		if (octtreeXthetas.rows() == 1)
		{
			XthetasNextLevel = octtreeXthetas.block(0, 0, 1, octtreeNthetas(0));
			XphisNextLevel = octtreeXphis.block(0, 0, 1, octtreeNphis(0));
		}
		else
		{
			XthetasNextLevel = octtreeXthetas.block(1, 0, 1, octtreeNthetas(1));
			XphisNextLevel = octtreeXphis.block(1, 0, 1, octtreeNphis(1));
		}
		levels.push_back(Level(
			N_levels + 1, LExpansion(0), leaf_side_length, this->big_cube_lower_coord,
			cubes_centroids, cubes_centroids_obs, k, XthetasThisLevel,
			WthetasThisLevel, INCLUDED_THETA_BOUNDARIES, octtreeNthetas(0),
			PERIODIC_Theta, CYCLIC_Theta, NOrderInterpTheta,
			XphisThisLevel, WphisThisLevel, INCLUDED_PHI_BOUNDARIES,
			octtreeNphis(0), PERIODIC_Phi, CYCLIC_Phi, NOrderInterpPhi,
			XthetasNextLevel, XphisNextLevel, VERBOSE
		));
		levels[0].sortCubesByParents();
		minCubesArraysSize = levels[0].getCubesSizeMB();
		indMinCubesArraysSize = 0;
	} // end of first level construction...using braces
	// construction of the coarser levels
	for (int j = 1; j < N_levels; ++j)
	{
		Array<double, 1, Dynamic> XthetasThisLevel, WthetasThisLevel, XphisThisLevel, WphisThisLevel;
		XthetasThisLevel = octtreeXthetas.block(j, 0, 1, octtreeNthetas(j));
		WthetasThisLevel = octtreeWthetas.block(j, 0, 1, octtreeNthetas(j));
		XphisThisLevel = octtreeXphis.block(j, 0, 1, octtreeNphis(j));
		WphisThisLevel = octtreeWphis.block(j, 0, 1, octtreeNphis(j));

		Array<double, 1, Dynamic> XthetasNextLevel, XphisNextLevel;
		if (j < N_levels - 1)
		{
			XthetasNextLevel = octtreeXthetas.block(j + 1, 0, 1, octtreeNthetas(j + 1));
			XphisNextLevel = octtreeXphis.block(j + 1, 0, 1, octtreeNphis(j + 1));
		}
		else
		{
			XthetasNextLevel = octtreeXthetas.block(j, 0, 1, octtreeNthetas(j));
			XphisNextLevel = octtreeXphis.block(j, 0, 1, octtreeNphis(j));
		}
		levels.push_back(Level(
			levels[j - 1], LExpansion(j), this->big_cube_lower_coord,
			XthetasThisLevel, WthetasThisLevel, INCLUDED_THETA_BOUNDARIES,
			octtreeNthetas(j), PERIODIC_Theta, CYCLIC_Theta, NOrderInterpTheta,
			XphisThisLevel, WphisThisLevel, INCLUDED_PHI_BOUNDARIES,
			octtreeNphis(j), PERIODIC_Phi, CYCLIC_Phi, NOrderInterpPhi,
			XthetasNextLevel, XphisNextLevel, VERBOSE
		));
		//		cout << octtreeXthetas.block(j, 0, 1, octtreeNthetas(j)).rows()
		//			<< " " << octtreeXthetas.block(j, 0, 1, octtreeNthetas(j)).cols() << endl;
		levels[j].sortCubesByParents();
		levels[j - 1].updateFatherIndexes(levels[j]);
		if (levels[j].getCubesSizeMB() < minCubesArraysSize)
		{
			minCubesArraysSize = levels[j].getCubesSizeMB();
			indMinCubesArraysSize = j;
		}
		if (j == N_levels - 1) levels[j].setCeiling(1); // we have reached the top level
	}

	if (ALLOW_CEILING_LEVEL == 1)
	{ //we remove the non-necessary levels
		int j = levels.size() - 1;
		while (j > indMinCubesArraysSize)
		{
			levels.pop_back();
			j--;
		}
		levels[indMinCubesArraysSize].setCeiling(1);
		vector<Level>(levels).swap(levels); // trick for trimming exceeding capacity
		N_levels = levels.size();
	}

	for (int j = 1; j < N_levels; j++) levels[j].searchCubesNeighborsIndexes();
	// alpha translations
	if (VERBOSE == 1)
	{
		cout << "Searching the indexes of the possible cubes "
			<< "for alpha translations.........." << endl;
	}
	cout << "N_levels = " << N_levels << endl;
	for (int l = 0; l < N_levels; l++) this->findAlphaTransParticipantsIndexes(l);

	N_levels = levels.size();
	if (VERBOSE == 1) cout << "Tree construction terminated." << endl;
}

Octtree::Octtree
(const string octtree_data_path,
 const string media_path,
 const ArrayXXd& cubes_centroids,
 const ArrayXXd& cubes_centroids_obs,
 const ArrayXi& target_surface,
 ArrayXi& octtreeNthetas, ArrayXXd& octtreeXthetas, ArrayXXd& octtreeWthetas,
 ArrayXi& octtreeNphis, ArrayXXd& octtreeXphis, ArrayXXd& octtreeWphis,
 ArrayXi& LExpansion, ArrayXcd& CFIEcoeffs,
 ArrayXd& bigCubeLowerCoord, ArrayXd& bigCubeCenterCoord)
{
	this->target_surface = target_surface;

	//	cout << octtreeXthetas << endl;
	//	cout << octtreeXthetas.block(0, 0, 1, octtreeNthetas(0)) << endl;

	octtreeDataPath = octtree_data_path;
	// verbose or not?
	readIntFromASCIIFile(octtree_data_path + "VERBOSE.txt", VERBOSE);

	if (VERBOSE == 1) cout << "creating the tree " << " from disk data" << endl;
	numberOfUpdates = 0;
	readIntFromASCIIFile(octtree_data_path + "N_active_levels.txt", N_levels);
	if (VERBOSE == 1) cout << "N active levels = " << N_levels << endl;
	readIntFromASCIIFile
		(octtree_data_path + "ALLOW_CEILING_LEVEL.txt", ALLOW_CEILING_LEVEL);
	if (VERBOSE == 1)
		cout << "ALLOW_CEILING_LEVEL = " << ALLOW_CEILING_LEVEL << endl;
	readIntFromASCIIFile
		(octtree_data_path + "N_GaussOnTriangle.txt", N_GaussOnTriangle);
	if (VERBOSE == 1) cout << "N_GaussOnTriangle = " << N_GaussOnTriangle << endl;

	// theta data
	int INCLUDED_THETA_BOUNDARIES, PERIODIC_Theta, CYCLIC_Theta, NOrderInterpTheta;
	readIntFromASCIIFile(octtree_data_path + "INCLUDED_THETA_BOUNDARIES.txt",
	                     INCLUDED_THETA_BOUNDARIES);
	readIntFromASCIIFile(octtree_data_path + "PERIODIC_Theta.txt", PERIODIC_Theta);
	readIntFromASCIIFile(octtree_data_path + "CYCLIC_Theta.txt", CYCLIC_Theta);
	readIntFromASCIIFile(octtree_data_path + "NOrderInterpTheta.txt",
	                     NOrderInterpTheta);

	// phi data
	int INCLUDED_PHI_BOUNDARIES, PERIODIC_Phi, CYCLIC_Phi, NOrderInterpPhi;
	readIntFromASCIIFile(octtree_data_path + "INCLUDED_PHI_BOUNDARIES.txt",
	                     INCLUDED_PHI_BOUNDARIES);
	readIntFromASCIIFile(octtree_data_path + "PERIODIC_Phi.txt", PERIODIC_Phi);
	readIntFromASCIIFile(octtree_data_path + "CYCLIC_Phi.txt", CYCLIC_Phi);
	readIntFromASCIIFile(octtree_data_path + "NOrderInterpPhi.txt",
	                     NOrderInterpPhi);

	double leaf_side_length;
	readDoubleFromASCIIFile
		(octtree_data_path + "leaf_side_length.txt", leaf_side_length);

	if (VERBOSE == 1) cout << "L expansions = " << LExpansion << endl;
	if (VERBOSE == 1) cout << "CFIE coeffs = " << CFIEcoeffs << endl;

	for (int i = 0; i < 3; ++i)
	{
		big_cube_lower_coord[i] = bigCubeLowerCoord(i);
		big_cube_center_coord[i] = bigCubeCenterCoord(i);
	}

	readDoubleFromASCIIFile(octtree_data_path + "w.txt", this->w);
	readComplexDoubleFromASCIIFile(octtree_data_path + media_path + "k.txt", this->k);
	readComplexDoubleFromASCIIFile(octtree_data_path + media_path + "eps_r.txt", this->eps_r);
	readComplexDoubleFromASCIIFile(octtree_data_path + media_path + "mu_r.txt", this->mu_r);

	cout << "w = " << w << ", mu_rel = " << mu_r
		<< ", eps_rel = " << eps_r << endl;

	CFIE.resize(CFIEcoeffs.size());
	CFIE = CFIEcoeffs;
	levels.reserve(N_levels); // reserve the memory
	double minCubesArraysSize;
	int indMinCubesArraysSize;
	{
		// construction of the first level
		Array<double, 1, Dynamic> XthetasThisLevel, WthetasThisLevel, XphisThisLevel, WphisThisLevel;
		XthetasThisLevel = octtreeXthetas.block(0, 0, 1, octtreeNthetas(0));
		WthetasThisLevel = octtreeWthetas.block(0, 0, 1, octtreeNthetas(0));
		XphisThisLevel = octtreeXphis.block(0, 0, 1, octtreeNphis(0));
		WphisThisLevel = octtreeWphis.block(0, 0, 1, octtreeNphis(0));

		Array<double, 1, Dynamic> XthetasNextLevel, XphisNextLevel;
		if (octtreeXthetas.rows() == 1)
		{
			XthetasNextLevel = octtreeXthetas.block(0, 0, 1, octtreeNthetas(0));
			XphisNextLevel = octtreeXphis.block(0, 0, 1, octtreeNphis(0));
		}
		else
		{
			XthetasNextLevel = octtreeXthetas.block(1, 0, 1, octtreeNthetas(1));
			XphisNextLevel = octtreeXphis.block(1, 0, 1, octtreeNphis(1));
		}
		levels.push_back(Level(N_levels + 1,
		                       LExpansion(0),
		                       leaf_side_length,
		                       this->big_cube_lower_coord,
		                       cubes_centroids,
		                       cubes_centroids_obs,
		                       k,
		                       XthetasThisLevel,
		                       WthetasThisLevel,
		                       INCLUDED_THETA_BOUNDARIES,
		                       octtreeNthetas(0),
		                       PERIODIC_Theta,
		                       CYCLIC_Theta,
		                       NOrderInterpTheta,
		                       XphisThisLevel,
		                       WphisThisLevel,
		                       INCLUDED_PHI_BOUNDARIES,
		                       octtreeNphis(0),
		                       PERIODIC_Phi,
		                       CYCLIC_Phi,
		                       NOrderInterpPhi,
		                       XthetasNextLevel,
		                       XphisNextLevel,
		                       VERBOSE));
		//		return;
		levels[0].sortCubesByParents();
		minCubesArraysSize = levels[0].getCubesSizeMB();
		indMinCubesArraysSize = 0;
	} // end of first level construction...using braces
	// construction of the coarser levels
	for (int j = 1; j < N_levels; ++j)
	{
		Array<double, 1, Dynamic> XthetasThisLevel, WthetasThisLevel, XphisThisLevel, WphisThisLevel;
		XthetasThisLevel = octtreeXthetas.block(j, 0, 1, octtreeNthetas(j));
		WthetasThisLevel = octtreeWthetas.block(j, 0, 1, octtreeNthetas(j));
		XphisThisLevel = octtreeXphis.block(j, 0, 1, octtreeNphis(j));
		WphisThisLevel = octtreeWphis.block(j, 0, 1, octtreeNphis(j));

		Array<double, 1, Dynamic> XthetasNextLevel, XphisNextLevel;
		if (j < N_levels - 1)
		{
			XthetasNextLevel = octtreeXthetas.block(j + 1, 0, 1, octtreeNthetas(j + 1));
			XphisNextLevel = octtreeXphis.block(j + 1, 0, 1, octtreeNphis(j + 1));
		}
		else
		{
			XthetasNextLevel = octtreeXthetas.block(j, 0, 1, octtreeNthetas(j));
			XphisNextLevel = octtreeXphis.block(j, 0, 1, octtreeNphis(j));
		}
		levels.push_back(Level(
			levels[j - 1], LExpansion(j), this->big_cube_lower_coord,
			XthetasThisLevel, WthetasThisLevel, INCLUDED_THETA_BOUNDARIES,
			octtreeNthetas(j), PERIODIC_Theta, CYCLIC_Theta, NOrderInterpTheta,
			XphisThisLevel, WphisThisLevel, INCLUDED_PHI_BOUNDARIES,
			octtreeNphis(j), PERIODIC_Phi, CYCLIC_Phi, NOrderInterpPhi,
			XthetasNextLevel, XphisNextLevel, VERBOSE
		));
		//		cout << octtreeXthetas.block(j, 0, 1, octtreeNthetas(j)).rows()
		//			<< " " << octtreeXthetas.block(j, 0, 1, octtreeNthetas(j)).cols() << endl;
		levels[j].sortCubesByParents();
		levels[j - 1].updateFatherIndexes(levels[j]);
		if (levels[j].getCubesSizeMB() < minCubesArraysSize)
		{
			minCubesArraysSize = levels[j].getCubesSizeMB();
			indMinCubesArraysSize = j;
		}
		if (j == N_levels - 1) levels[j].setCeiling(1); // we have reached the top level
	}

	if (ALLOW_CEILING_LEVEL == 1)
	{ //we remove the non-necessary levels
		int j = levels.size() - 1;
		while (j > indMinCubesArraysSize)
		{
			levels.pop_back();
			j--;
		}
		levels[indMinCubesArraysSize].setCeiling(1);
		vector<Level>(levels).swap(levels); // trick for trimming exceeding capacity
		N_levels = levels.size();
	}

	for (int j = 1; j < N_levels; j++) levels[j].searchCubesNeighborsIndexes();
	// alpha translations
	if (VERBOSE == 1)
	{
		cout << "Searching the indexes of the possible cubes "
			<< "for alpha translations.........." << endl;
	}
	cout << "N_levels = " << N_levels << endl;
	for (int l = 0; l < N_levels; l++) this->findAlphaTransParticipantsIndexes(l);

	N_levels = levels.size();
	if (VERBOSE == 1) cout << "Tree construction terminated." << endl;
}

// void Octtree::writeAssignedLeafCubesToDisk
// (const string path, const string filename)
// {
// 	const int NC = levels[0].cubes.size();
// 	Array<int, 2> cubesNumberToProcessNumber(NC, 3);
// 	for (int i = 0 ; i < NC ; ++i) {
// 		cubesNumberToProcessNumber(i, 0) = levels[0].cubes[i].getOldIndex();
// 		cubesNumberToProcessNumber(i, 1) = levels[0].cubes[i].getNumber();
// 		cubesNumberToProcessNumber(i, 2) = 0;
// 	}
// 	writeIntBlitzArray2DToASCIIFile(path + filename, cubesNumberToProcessNumber);
// }

void Octtree::computeIndexesOfCubesInOriginalMesh
(ArrayXi& oldIndexesOfCubes,
 ArrayXi& oldIndexesOfCubes_obs)
{
	levels[0].computeOldIndexesOfCubes(oldIndexesOfCubes, oldIndexesOfCubes_obs);
}

void Octtree::computeGaussLocatedArguments
(const ArrayXi& local_cubes_NRWG,
 const ArrayXi& local_RWG_numbers,
 const ArrayXi& local_RWG_Numbers_surf,
 const ArrayXXi& local_RWGNumbers_signedTriangles,
 const ArrayXXd& local_RWGNumbers_trianglesCoord)
{
	if (VERBOSE == 1)
		cout << "computing the leaf level Gauss located arguments.........."
			<< endl;
	levels[0].computeGaussLocatedArguments
	(local_cubes_NRWG, local_RWG_numbers,
	 local_RWG_Numbers_surf, local_RWGNumbers_signedTriangles,
	 local_RWGNumbers_trianglesCoord, this->N_GaussOnTriangle);
}

void Octtree::computeObsLocatedArguments
(const ArrayXi& cubes_Nobs,
 const ArrayXi& obs_numbers,
 const ArrayXXd& obsNumber_obsCoord)
{
	if (VERBOSE == 1)
		cout << "computing the leaf level obsevation located arguments....."
			<< endl;
	levels[0].computeObsLocatedArguments(cubes_Nobs, obs_numbers, obsNumber_obsCoord);
}

void Octtree::RWGs_renumbering(void)
{
	levels[0].RWGs_renumbering();
}

void Octtree::constructArrays(void)
{
	const int N_levels = levels.size();
	cout << "computing the alpha translations and shifting arrays.........."
		<< endl;

	double alphaTranslation_smoothing_factor;
	double alphaTranslation_thresholdRelValueMax;
	double alphaTranslation_RelativeCountAboveThreshold;

	readDoubleFromASCIIFile(this->octtreeDataPath
	                        + "alphaTranslation_smoothing_factor.txt",
	                        alphaTranslation_smoothing_factor);
	readDoubleFromASCIIFile(this->octtreeDataPath
	                        + "alphaTranslation_thresholdRelValueMax.txt",
	                        alphaTranslation_thresholdRelValueMax);
	readDoubleFromASCIIFile(this->octtreeDataPath
	                        + "alphaTranslation_RelativeCountAboveThreshold.txt",
	                        alphaTranslation_RelativeCountAboveThreshold);

	if (alphaTranslation_smoothing_factor < 1.0)
		alphaTranslation_smoothing_factor = 1.0;
	if (alphaTranslation_smoothing_factor > 2.0)
		alphaTranslation_smoothing_factor = 2.0;
	if (alphaTranslation_RelativeCountAboveThreshold > 1.0)
		alphaTranslation_RelativeCountAboveThreshold = 1.0;
	if (alphaTranslation_RelativeCountAboveThreshold < 0.0)
		alphaTranslation_RelativeCountAboveThreshold = 0.0;

	for (int j = 0; j < N_levels; ++j)
	{
		levels[j].NCubesXYZComputation();
		levels[j].alphaTranslationsComputation(alphaTranslation_smoothing_factor,
		                                       alphaTranslation_thresholdRelValueMax,
		                                       alphaTranslation_RelativeCountAboveThreshold);
		levels[j].shiftingArraysComputation();
	}

	for (int j = 0; j < N_levels; ++j)
	{
		cout << "for level " << levels[j].getLevel()
			<< ", sizes of arrays are:" << endl;
		cout << "  size of alphaTranslations = "
			<< levels[j].getAlphaTranslationsSizeMB() << " MB;" << endl;
		cout << "  size of alphaTranslationsIndexes = "
			<< levels[j].getAlphaTranslationsIndexesSizeMB() << " MB;" << endl;
		cout << "  size of cubes = "
			<< levels[j].getCubesSizeMB() << " MB, number of cubes = "
			<< levels[j].cubes.size() << endl;
		cout << "  size of indexes for alphaTranslations candidates = "
			<< levels[j].getSizeMBOfAlphaTransParticipantsIndexes()
			<< " MB;" << endl;
		cout << "  size of shiftingArrays = "
			<< levels[j].getShiftingArraysSizeMB() << " MB;" << endl;
	}
}

void Octtree::copyOcttree(const Octtree& octtreeTocopy) /// copy constructor
{
	cout << "\ncopying octtree... " << endl;
	target_surface.resize(octtreeTocopy.target_surface.size());
	target_surface = octtreeTocopy.target_surface;
	numberOfUpdates = octtreeTocopy.numberOfUpdates;
	L = octtreeTocopy.L;
	k = octtreeTocopy.k;
	N_GaussOnTriangle = octtreeTocopy.N_GaussOnTriangle;
	eps_r = octtreeTocopy.eps_r;
	mu_r = octtreeTocopy.mu_r;
	w = octtreeTocopy.w;
	octtreeDataPath = octtreeTocopy.octtreeDataPath;
	int NLevels = octtreeTocopy.N_levels;
	cout << "The number of Levels is " << NLevels << endl;
	levels.resize(NLevels);
	for (int j = 0; j < NLevels; j++)
		levels[j].copyLevel(octtreeTocopy.getLevel(j));
	for (int i = 0; i < 3; i++)
		big_cube_lower_coord[i] = octtreeTocopy.big_cube_lower_coord[i];
	for (int i = 0; i < 3; i++)
		big_cube_center_coord[i] = octtreeTocopy.big_cube_center_coord[i];
	CFIE.resize(octtreeTocopy.CFIE.size());
	CFIE = octtreeTocopy.CFIE;
	N_levels = octtreeTocopy.N_levels;
	ALLOW_CEILING_LEVEL = octtreeTocopy.ALLOW_CEILING_LEVEL;
	VERBOSE = octtreeTocopy.VERBOSE;
	cout << "end of copying octtree... " << endl;
}

Octtree::Octtree(const Octtree& octtreeToCopy) /// copy constructor
{
	copyOcttree(octtreeToCopy);
}

Octtree& Octtree::operator=(const Octtree& octtreeToCopy) /// copy assignment
{
	copyOcttree(octtreeToCopy);
	return *this;
}

Octtree::~Octtree()
{
	levels.clear();
	//	CFIE.free();
}

vector<int> Octtree::getNeighborsSonsIndexes(const int index, const int l) const
{
	vector<int> sonsOfNeighbors, sonsOfNeighborsTmp;
	const vector<int> neighborsIndexes(getCubeLevel(index, l).neighborsIndexes);
	for (unsigned int m = 0; m < neighborsIndexes.size(); ++m)
	{
		sonsOfNeighborsTmp = getCubeLevel(neighborsIndexes[m], l).sonsIndexes;
		for (unsigned int n = 0; n < sonsOfNeighborsTmp.size(); ++n) sonsOfNeighbors.push_back(sonsOfNeighborsTmp[n]);
	}
	return sonsOfNeighbors;
}

vector<int> Octtree::getNeighborsSonsIndexes_obs(const int index, const int l) const
{
	vector<int> sonsOfNeighbors, sonsOfNeighborsTmp;
	const vector<int> neighborsIndexes(levels[l].cubes_obs[index].neighborsIndexes);
	for (unsigned int m = 0; m < neighborsIndexes.size(); ++m)
	{
		sonsOfNeighborsTmp = getCubeLevel(neighborsIndexes[m], l).sonsIndexes;
		for (unsigned int n = 0; n < sonsOfNeighborsTmp.size(); ++n) sonsOfNeighbors.push_back(sonsOfNeighborsTmp[n]);
	}
	return sonsOfNeighbors;
}

void Octtree::findAlphaTransParticipantsIndexes(const int l)
{
	const int N_cubes = levels[l].getLevelSize();
	// we treat the ceiling level differently than the regular levels
	// hereafter is the code for the ceiling level
	if (l == N_levels - 1)
	{ // if we are at the ceiling level

		for (int i = 0; i < N_cubes; ++i)
		{ // loop on the local cubes

			vector<int> AlphaTransParticipantsIndexes;

			for (int j = 0; j < N_cubes; ++j)
			{ // loop on all the cubes (because ceiling level)

				const double* diffAbsCartCoord_1(levels[l].cubes[i].absoluteCartesianCoord);
				const double* diffAbsCartCoord_2(levels[l].cubes[j].absoluteCartesianCoord);
				const double diffAbsCartCoord[3]
					= {diffAbsCartCoord_1[0] - diffAbsCartCoord_2[0],
						diffAbsCartCoord_1[1] - diffAbsCartCoord_2[1],
						diffAbsCartCoord_1[2] - diffAbsCartCoord_2[2]
					};

				bool condition = false;
				// condition = true if cubes are not touching

				for (int mm = 0; mm < 3; ++mm)
					condition = (condition || (abs(diffAbsCartCoord[mm]) > 1.0));
				if (condition) AlphaTransParticipantsIndexes.push_back(j);
			} // end for

			// we now trim the excess capacity of alphaTransParticipantsIndexes
			levels[l].cubes[i].AlphaTransParticipantsIndexes.resize
			                  (AlphaTransParticipantsIndexes.size());

			for (unsigned int j = 0; j < AlphaTransParticipantsIndexes.size(); ++j)
				levels[l].cubes[i].AlphaTransParticipantsIndexes[j]
					= AlphaTransParticipantsIndexes[j];
		}
	}

	else
	{ // for the NON-CEILING level

		for (int i = 0; i < N_cubes; ++i)
		{
			vector<int> AlphaTransParticipantsIndexes;
			const vector<int> possibleIndexes
				(getNeighborsSonsIndexes(levels[l].cubes[i].getFatherIndex(), l + 1));

			for (unsigned int j = 0; j < possibleIndexes.size(); j++)
			{
				// possible indexes of the alpha trans participants
				const int possibleIndex = possibleIndexes[j];
				const double* diffAbsCartCoord_1(levels[l].cubes[i].absoluteCartesianCoord);
				const double* diffAbsCartCoord_2(levels[l].cubes[possibleIndex].absoluteCartesianCoord);
				const double diffAbsCartCoord[3] = {diffAbsCartCoord_1[0] - diffAbsCartCoord_2[0],
					diffAbsCartCoord_1[1] - diffAbsCartCoord_2[1],
					diffAbsCartCoord_1[2] - diffAbsCartCoord_2[2]
				};

				bool condition = false;

				// condition = true if cubes are not touching
				for (int mm = 0; mm < 3; ++mm)
					condition = (condition || (abs(diffAbsCartCoord[mm]) > 1.0));
				if (condition)
				{
					AlphaTransParticipantsIndexes.push_back(possibleIndex);
				}
			}

			// we now trim the excess capacity of alphaTransParticipantsIndexes
			levels[l].cubes[i].AlphaTransParticipantsIndexes.resize
			                  (AlphaTransParticipantsIndexes.size());

			for (unsigned int j = 0; j < AlphaTransParticipantsIndexes.size(); ++j)
				levels[l].cubes[i].AlphaTransParticipantsIndexes[j]
					= AlphaTransParticipantsIndexes[j];
		}
	}

	const int N_cubes_obs = levels[l].cubes_obs.size();
	if (l == N_levels - 1)
	{ // if we are at the ceiling level

		for (int i = 0; i < N_cubes_obs; ++i)
		{ // loop on the local cubes

			vector<int> AlphaTransParticipantsIndexes;

			for (int j = 0; j < N_cubes; ++j)
			{ // loop on all the cubes (because ceiling level)

				const double* diffAbsCartCoord_1(levels[l].cubes_obs[i].absoluteCartesianCoord);
				const double* diffAbsCartCoord_2(levels[l].cubes[j].absoluteCartesianCoord);
				const double diffAbsCartCoord[3]
					= {diffAbsCartCoord_1[0] - diffAbsCartCoord_2[0],
						diffAbsCartCoord_1[1] - diffAbsCartCoord_2[1],
						diffAbsCartCoord_1[2] - diffAbsCartCoord_2[2]
					};

				bool condition = false;
				// condition = true if cubes are not touching
				for (int mm = 0; mm < 3; ++mm)
					condition = (condition || (abs(diffAbsCartCoord[mm]) > 1.0));
				if (condition)
				{
					AlphaTransParticipantsIndexes.push_back(j);
					levels[l].cubes[j].AlphaTransParticipantsIndexes_obs.push_back(i);
				}
			} // end for

			// we now trim the excess capacity of alphaTransParticipantsIndexes
			levels[l].cubes_obs[i].AlphaTransParticipantsIndexes.resize
			                      (AlphaTransParticipantsIndexes.size());
			for (unsigned int j = 0; j < AlphaTransParticipantsIndexes.size(); ++j)
				levels[l].cubes_obs[i].AlphaTransParticipantsIndexes[j]
					= AlphaTransParticipantsIndexes[j];
		}
	}

	else
	{ // for the NON-CEILING level

		for (int i = 0; i < N_cubes_obs; ++i)
		{
			vector<int> AlphaTransParticipantsIndexes;
			const vector<int> possibleIndexes
				(getNeighborsSonsIndexes_obs(levels[l].cubes_obs[i].getFatherIndex(), l + 1));

			for (unsigned int j = 0; j < possibleIndexes.size(); j++)
			{
				// possible indexes of the alpha trans participants
				const int possibleIndex = possibleIndexes[j];
				const double* diffAbsCartCoord_1(levels[l].cubes_obs[i].absoluteCartesianCoord);
				const double* diffAbsCartCoord_2(levels[l].cubes[possibleIndex].absoluteCartesianCoord);
				const double diffAbsCartCoord[3]
					= {diffAbsCartCoord_1[0] - diffAbsCartCoord_2[0],
						diffAbsCartCoord_1[1] - diffAbsCartCoord_2[1],
						diffAbsCartCoord_1[2] - diffAbsCartCoord_2[2]
					};

				bool condition = false;
				// condition = true if cubes are not touching
				for (int mm = 0; mm < 3; ++mm)
					condition = (condition || (abs(diffAbsCartCoord[mm]) > 1.0));
				if (condition)
				{
					AlphaTransParticipantsIndexes.push_back(possibleIndex);
					levels[l].cubes[possibleIndex].AlphaTransParticipantsIndexes_obs.push_back(i);
				}
			}

			// we now trim the excess capacity of alphaTransParticipantsIndexes
			levels[l].cubes_obs[i].AlphaTransParticipantsIndexes.resize
			                      (AlphaTransParticipantsIndexes.size());
			for (unsigned int j = 0; j < AlphaTransParticipantsIndexes.size(); ++j)
				levels[l].cubes_obs[i].AlphaTransParticipantsIndexes[j]
					= AlphaTransParticipantsIndexes[j];
		}
	}
}

void Octtree::shiftExp(ArrayXXcd& S,
                       const ArrayXcd& shiftingArray)
{
	const int N = shiftingArray.size();
	for (int i = 0; i < N; i++)
	{
		S(0, i) *= shiftingArray(i);
		S(1, i) *= shiftingArray(i);
	}
}

void Octtree::updateSup(const ArrayXcd& I_PQ,
                        const bool is_transpose)
{
	const int N_levels = levels.size();
	// first we update the Sups...which are stored (temporarily) in the Sdowns!!
	for (int l = 0; l < N_levels; ++l)
	{
		const int N_cubes = levels[l].getLevelSize();
		const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
		const int N_directions = N_theta * N_phi;
		if (levels[l].Sdown.size() == 0) levels[l].Sdown.resize(N_cubes);
		// we first compute the Sups of all the cubes at the given level
		if (l == 0)
		{ // we use computeSup only for the leaf cubes
			for (int i = 0; i < N_cubes; ++i)
			{
				if (levels[l].Sdown(i).size() == 0)
					levels[l].Sdown(i).resize(2, N_theta * N_phi);
				levels[l].computeSup(levels[l].Sdown(i), k, I_PQ, levels[l].cubes[i],
				                     levels[l].thetas, levels[l].phis, target_surface);
			}
		}
		else
		{ // the Sups are obtained by interpolation from the sonsCubes Sups
			ArrayXXcd S_tmp(2, N_directions);
			for (int i = 0; i < N_cubes; ++i)
			{
				levels[l].Sdown(i) = ArrayXXcd::Zero(2, N_directions);
				const int NSons = levels[l].cubes[i].sonsIndexes.size();
				for (int j = 0; j < NSons; ++j)
				{
					const int sonIndex = levels[l].cubes[i].sonsIndexes[j];
					// interpolation
					ArrayXcd S_tmp_row(S_tmp.cols());
					interpolate2Dlfi(S_tmp_row,
					                 levels[l - 1].Sdown(sonIndex).row(0),
					                 levels[l - 1].lfi2D);
					S_tmp.row(0) = S_tmp_row;
					interpolate2Dlfi(S_tmp_row,
					                 levels[l - 1].Sdown(sonIndex).row(1),
					                 levels[l - 1].lfi2D);
					S_tmp.row(1) = S_tmp_row;
					// shifting
					const double* rc_1(levels[l].cubes[i].rCenter);
					const double* rc_2(levels[l - 1].cubes[sonIndex].rCenter);
					double DRcenters[3] = {rc_1[0] - rc_2[0],
						rc_1[1] - rc_2[1],
						rc_1[2] - rc_2[2]
					};
					shiftExp(S_tmp, levels[l].getShiftingArray(DRcenters[0],
					                                           DRcenters[1],
					                                           DRcenters[2]));
					levels[l].Sdown(i) += S_tmp;
				}
			}
		}
	}
	//	{
	//		int l = 0, i = 0;
	//		for (int j = 0; j < 3; j++) cout << levels[l].cubes[i].rCenter[j] << " ";
	//		cout << endl;
	//		cout << levels[l].Sdown(i) << endl;
	//	}
}

void Octtree::updateSup_obs(const ArrayXXcd& E_obs,
                            const bool is_transpose, const bool M_current)
{
	const int N_levels = levels.size();

	for (int l = 0; l < N_levels; ++l)
	{
		const int N_cubes = levels[l].cubes_obs.size();
		const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
		const int N_directions = N_theta * N_phi;

		if (levels[l].Sdown_obs.size() == 0) levels[l].Sdown_obs.resize(N_cubes);

		// we first compute the Sups of all the cubes at the given level
		if (l == 0)
		{ // we use computeSup only for the leaf cubes
			for (int i = 0; i < N_cubes; ++i)
			{
				if (levels[l].Sdown_obs(i).size() == 0)
					levels[l].Sdown_obs(i).resize(2, N_theta * N_phi);
				levels[l].computeSup_obs(levels[l].Sdown_obs(i), k, E_obs,
				                         levels[l].cubes_obs[i], levels[l].thetas,
				                         levels[l].phis, M_current);
			}
		}

		else
		{ // the Sups are obtained by interpolation from the sonsCubes Sups
			ArrayXXcd S_tmp(2, N_directions);

			for (int i = 0; i < N_cubes; ++i)
			{
				if (levels[l].Sdown_obs(i).size() == 0)
					levels[l].Sdown_obs(i).resize(2, N_theta * N_phi);

				levels[l].Sdown_obs(i) = ArrayXXcd::Zero(levels[l].Sdown_obs(i).rows(),
				                                         levels[l].Sdown_obs(i).cols());
				const int NSons = levels[l].cubes_obs[i].sonsIndexes.size();

				for (int j = 0; j < NSons; ++j)
				{
					const int sonIndex = levels[l].cubes_obs[i].sonsIndexes[j];
					// interpolation
					ArrayXcd S_tmp_row(S_tmp.cols());
					interpolate2Dlfi(S_tmp_row,
					                 levels[l - 1].Sdown_obs(sonIndex).row(0),
					                 levels[l - 1].lfi2D);
					S_tmp.row(0) = S_tmp_row;
					interpolate2Dlfi(S_tmp_row,
					                 levels[l - 1].Sdown_obs(sonIndex).row(1),
					                 levels[l - 1].lfi2D);
					S_tmp.row(1) = S_tmp_row;

					// shifting
					const double* rc_1(levels[l].cubes_obs[i].rCenter);
					const double* rc_2(levels[l - 1].cubes_obs[sonIndex].rCenter);
					double DRcenters[3] = {rc_1[0] - rc_2[0],
						rc_1[1] - rc_2[1], rc_1[2] - rc_2[2]
					};

					shiftExp(S_tmp, levels[l].getShiftingArray(DRcenters[0],
					                                           DRcenters[1],
					                                           DRcenters[2]));

					levels[l].Sdown_obs(i) += S_tmp;
				}
			}
		}
	}
}

void Octtree::SupAlphaMultiplication
(ArrayXXcd& SupAlpha,
 const ArrayXXcd& Sup,
 const ArrayXcd& alphaTranslation,
 const ArrayXi& alphaTranslationIndexesNonZeros,
 const ArrayXi& alphaTranslationIndexes,
 const int alphaCartesianCoord[3])
{
	if ((abs(alphaCartesianCoord[0]) > 1)
		|| (abs(alphaCartesianCoord[1]) > 1)
		|| (abs(alphaCartesianCoord[2]) > 1))
	{
		if (alphaTranslationIndexesNonZeros.size() == 0)
		{
			const int N_alpha(alphaTranslationIndexes.size());
			for (int i = 0; i < N_alpha; ++i)
			{
				const int newIndex(alphaTranslationIndexes(i));
				SupAlpha(0, i) += Sup(0, i) * alphaTranslation(newIndex);
				SupAlpha(1, i) += Sup(1, i) * alphaTranslation(newIndex);
			}
		}
		else
		{
			const int N_alpha(alphaTranslationIndexesNonZeros.size());
			for (int i = 0; i < N_alpha; ++i)
			{
				const int oldIndex = alphaTranslationIndexesNonZeros(i);
				const int newIndex = alphaTranslationIndexes(oldIndex);
				SupAlpha(0, newIndex) += Sup(0, newIndex) * alphaTranslation(i);
				SupAlpha(1, newIndex) += Sup(1, newIndex) * alphaTranslation(i);
			}
		}
	}
}

void Octtree::alphaTranslationsToCube
(ArrayXXcd& S_tmp,
 const Array<ArrayXXcd, Dynamic, 1>& LevelSup,
 const int l,
 const int cubeIndex,
 const vector<int>& indexesAlphaParticipants,
 const bool is_transpose,
 const bool is_obs)
{
	double cartCoord_1[3];
	if (is_obs) for (int j = 0; j < 3; ++j) cartCoord_1[j] = levels[l].cubes_obs[cubeIndex].absoluteCartesianCoord[j];
	else for (int j = 0; j < 3; ++j) cartCoord_1[j] = levels[l].cubes[cubeIndex].absoluteCartesianCoord[j];
	S_tmp = ArrayXXcd::Zero(S_tmp.rows(), S_tmp.cols());
	const int N_part = indexesAlphaParticipants.size();
	for (int j = 0; j < N_part; ++j)
	{
		const int indexParticipant = indexesAlphaParticipants[j];
		const double* cartCoord_2(levels[l].cubes[indexParticipant].absoluteCartesianCoord);
		double DRcenters[3] = {cartCoord_1[0] - cartCoord_2[0], cartCoord_1[1] - cartCoord_2[1], cartCoord_1[2] - cartCoord_2[2]};
		if (is_transpose) for (int i = 0; i < 3; ++i) DRcenters[i] = -DRcenters[i];
		const int alphaCartesianCoord[3] = {static_cast<int>(round(DRcenters[0])), static_cast<int>(round(DRcenters[1])), static_cast<int>(round(DRcenters[2]))};
		const int X = 1 * (alphaCartesianCoord[0] >= 0), Y = 1 * (alphaCartesianCoord[1] >= 0), Z = 1 * (alphaCartesianCoord[2] >= 0);
		const int m = abs(alphaCartesianCoord[0]), n = abs(alphaCartesianCoord[1]), p = abs(alphaCartesianCoord[2]);
		ArrayXcd alphaTranslation(levels[l].alphaTranslations[m][n][p].size());
		if (is_transpose) alphaTranslation = conj(levels[l].alphaTranslations[m][n][p]);
		else alphaTranslation = levels[l].alphaTranslations[m][n][p];
		SupAlphaMultiplication(S_tmp, LevelSup(indexParticipant), alphaTranslation,
		                       levels[l].alphaTranslationsIndexesNonZeros[m][n][p],
		                       levels[l].alphaTranslationsIndexes[X][Y][Z], alphaCartesianCoord);
	}
}

void Octtree::alphaTranslationsToCube_obs
(ArrayXXcd& S_tmp,
 const Array<ArrayXXcd, Dynamic, 1>& LevelSup,
 const int l,
 const int cubeIndex,
 const vector<int>& indexesAlphaParticipants)
{
	double cartCoord_1[3];
	for (int j = 0; j < 3; ++j) cartCoord_1[j] = levels[l].cubes[cubeIndex].absoluteCartesianCoord[j];
	S_tmp = ArrayXXcd::Zero(S_tmp.rows(), S_tmp.cols());
	const int N_part = indexesAlphaParticipants.size();
	//	cout << N_part << endl;
	for (int j = 0; j < N_part; ++j)
	{
		const int indexParticipant = indexesAlphaParticipants[j];
		const double* cartCoord_2(levels[l].cubes_obs[indexParticipant].absoluteCartesianCoord);
		double DRcenters[3] = {cartCoord_1[0] - cartCoord_2[0], cartCoord_1[1] - cartCoord_2[1], cartCoord_1[2] - cartCoord_2[2]};
		for (int i = 0; i < 3; ++i) DRcenters[i] = -DRcenters[i];
		const int alphaCartesianCoord[3] = {static_cast<int>(round(DRcenters[0])), static_cast<int>(round(DRcenters[1])), static_cast<int>(round(DRcenters[2]))};
		const int X = 1 * (alphaCartesianCoord[0] >= 0), Y = 1 * (alphaCartesianCoord[1] >= 0), Z = 1 * (alphaCartesianCoord[2] >= 0);
		const int m = abs(alphaCartesianCoord[0]), n = abs(alphaCartesianCoord[1]), p = abs(alphaCartesianCoord[2]);
		ArrayXcd alphaTranslation(levels[l].alphaTranslations[m][n][p].size());
		alphaTranslation = conj(levels[l].alphaTranslations[m][n][p]);
		SupAlphaMultiplication(S_tmp, LevelSup(indexParticipant), alphaTranslation, levels[l].alphaTranslationsIndexesNonZeros[m][n][p],
		                       levels[l].alphaTranslationsIndexes[X][Y][Z], alphaCartesianCoord);
		//		if (j == 0) cout << "test " << alphaTranslation << endl;
	}
}

void Octtree::alphaTranslations(const bool is_transpose)
{
	// now the translation stage!!!
	// cout << "alpha Translations..."; flush(cout);
	for (int l = 0; l < N_levels; ++l)
	{
		const int N_cubes = levels[l].getLevelSize();
		const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
		const int N_directions = N_theta * N_phi;
		// we must define a SupThisLevel
		Array<ArrayXXcd, Dynamic, 1> SupThisLevel(N_cubes);
		for (int i = 0; i < N_cubes; ++i)
		{
			SupThisLevel(i).resize(2, N_directions);
			SupThisLevel(i) = levels[l].Sdown(i);
			levels[l].Sdown(i) = ArrayXXcd::Zero(levels[l].Sdown(i).rows(),
			                                     levels[l].Sdown(i).cols());
		}
		for (int i = 0; i < N_cubes; ++i)
		{
			alphaTranslationsToCube(levels[l].Sdown(i), SupThisLevel, l, i,
			                        levels[l].cubes[i].AlphaTransParticipantsIndexes,
			                        is_transpose, false);
		}
	}
}

void Octtree::ZIFarComputation
(ArrayXcd& ZI,
 const ArrayXcd& I_PQ,
 const bool is_transpose,
 const bool M_current)
{
	// update of all the Sup of the tree
	this->updateSup(I_PQ, is_transpose);
	this->alphaTranslations(is_transpose);
	// cout << "tree descent"; flush(cout);
	const int N_levels = levels.size(), L = N_levels - 1;
	int thisLevel = L;

	// now the descent towards the bottom
	while (thisLevel > 0)
	{
		const int sonLevel = thisLevel - 1;
		const int N_cubes = levels[thisLevel].getLevelSize();
		ArrayXXcd Stmp2(2, levels[thisLevel].thetas.size() * levels[thisLevel].phis.size());
		ArrayXXcd Stmp3(2, levels[sonLevel].thetas.size() * levels[sonLevel].phis.size());
		for (int i = 0; i < N_cubes; ++i)
		{
			vector<int> sonsIndexes = levels[thisLevel].cubes[i].sonsIndexes;
			for (unsigned int j = 0; j < sonsIndexes.size(); ++j)
			{
				const int sonIndex = sonsIndexes[j];
				// shifting
				const double* rc_1(levels[sonLevel].cubes[sonIndex].rCenter);
				const double* rc_2(levels[thisLevel].cubes[i].rCenter);
				double DRcenters[3] = {rc_1[0] - rc_2[0], rc_1[1] - rc_2[1], rc_1[2] - rc_2[2]};
				Stmp2 = levels[thisLevel].Sdown(i);
				shiftExp(Stmp2, levels[thisLevel].getShiftingArray(DRcenters[0], DRcenters[1], DRcenters[2]));
				// anterpolate
				ArrayXcd Stmp3_row(Stmp3.cols());
				for (int m = 0; m < 2; m++)
				{
					anterpolate2Dlfi(Stmp3_row, Stmp2.row(m), levels[sonLevel].lfi2D);
					Stmp3.row(m) = Stmp3_row;
				}
				levels[sonLevel].Sdown(sonIndex) += Stmp3;
			}
		}
		thisLevel--;
	}
	// and finally the integration
	const int N_cubes = levels[thisLevel].getLevelSize();
	for (int i = 0; i < N_cubes; ++i)
	{
		levels[thisLevel].sphericalIntegration
		(ZI, levels[thisLevel].Sdown(i), levels[thisLevel].cubes[i],
		 levels[thisLevel].thetas, levels[thisLevel].phis, w, mu_r, k,
		 CFIE, M_current, target_surface);
	}
}

void Octtree::alphaTranslations_obs()
{
	// now the translation stage!!!
	// cout << "observation alpha Translations..."; flush(cout);
	for (int l = 0; l < N_levels; ++l)
	{
		const int N_cubes = levels[l].cubes.size();
		const int N_cubes_obs = levels[l].cubes_obs.size();
		const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
		const int N_directions = N_theta * N_phi;
		// we must define a SupThisLevel
		Array<ArrayXXcd, Dynamic, 1> SupThisLevel(N_cubes);
		for (int i = 0; i < N_cubes; ++i)
		{
			SupThisLevel(i).resize(2, N_directions);
			SupThisLevel(i) = levels[l].Sdown(i);
			levels[l].Sdown(i) = ArrayXXcd::Zero(levels[l].Sdown(i).rows(),
			                                     levels[l].Sdown(i).cols());
		}
		for (int i = 0; i < N_cubes_obs; ++i)
		{
			alphaTranslationsToCube(levels[l].Sdown_obs(i), SupThisLevel, l, i, levels[l].cubes_obs[i].AlphaTransParticipantsIndexes, false, true);
		}
	}
}

void Octtree::alphaTranslations_obs_transpose()
{
	for (int l = 0; l < N_levels; ++l)
	{
		const int N_cubes = levels[l].cubes.size();
		const int N_cubes_obs = levels[l].cubes_obs.size();
		const int N_theta = levels[l].thetas.size();
		const int N_phi = levels[l].phis.size();
		const int N_directions = N_theta * N_phi;

		// we must define a SupThisLevel
		Array<ArrayXXcd, Dynamic, 1> SupThisLevel(N_cubes_obs);

		for (int i = 0; i < N_cubes_obs; ++i)
		{
			SupThisLevel(i).resize(2, N_directions);
			SupThisLevel(i) = levels[l].Sdown_obs(i);
			levels[l].Sdown_obs(i) = ArrayXXcd::Zero(levels[l].Sdown_obs(i).rows(),
			                                         levels[l].Sdown_obs(i).cols());
			//			if (i == 0 && l == 0)
			//			{
			//				cout << SupThisLevel(i) << endl;
			//			}
		}

		for (int i = 0; i < N_cubes; ++i)
		{
			levels[l].Sdown(i) = ArrayXXcd::Zero(levels[l].Sdown(i).rows(),
			                                     levels[l].Sdown(i).cols());
			alphaTranslationsToCube_obs
			(levels[l].Sdown(i), SupThisLevel, l, i,
			 levels[l].cubes[i].AlphaTransParticipantsIndexes_obs);
			//			if (i == 0 && l == 0)
			//			{
			//				cout << levels[l].Sdown(i) << endl;
			//			}
		}
	}
}

void Octtree::ZINearComputation(ArrayXXcd& E_obs,
                                const ArrayXcd& I_PQ,
                                const bool M_current)
{
	for (int l = 0; l < N_levels; ++l)
	{
		const int N_cubes = levels[l].cubes_obs.size();
		const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
		const int N_directions = N_theta * N_phi;
		levels[l].Sdown_obs.resize(N_cubes);
		for (int i = 0; i < N_cubes; ++i)
		{
			levels[l].Sdown_obs(i).resize(2, N_directions);
		}
	}

	this->updateSup(I_PQ, false);
	this->alphaTranslations_obs();

	const int N_levels = levels.size(), L = N_levels - 1;
	int thisLevel = L;

	// now the descent towards the bottom
	while (thisLevel > 0)
	{
		const int sonLevel = thisLevel - 1;
		const int N_cubes_obs = levels[thisLevel].cubes_obs.size();
		ArrayXXcd Stmp2(2, levels[thisLevel].thetas.size() * levels[thisLevel].phis.size());
		ArrayXXcd Stmp3(2, levels[sonLevel].thetas.size() * levels[sonLevel].phis.size());
		for (int i = 0; i < N_cubes_obs; ++i)
		{
			vector<int> sonsIndexes = levels[thisLevel].cubes_obs[i].sonsIndexes;
			for (unsigned int j = 0; j < sonsIndexes.size(); ++j)
			{
				const int sonIndex = sonsIndexes[j];
				// shifting
				const double* rc_1(levels[sonLevel].cubes_obs[sonIndex].rCenter);
				const double* rc_2(levels[thisLevel].cubes_obs[i].rCenter);
				double DRcenters[3]
					= {rc_1[0] - rc_2[0], rc_1[1] - rc_2[1], rc_1[2] - rc_2[2]};
				Stmp2 = levels[thisLevel].Sdown_obs(i);
				shiftExp(Stmp2, levels[thisLevel].getShiftingArray
				         (DRcenters[0], DRcenters[1], DRcenters[2]));
				// anterpolate
				ArrayXcd Stmp3_row(Stmp3.cols());
				for (int m = 0; m < 2; m++)
				{
					anterpolate2Dlfi(Stmp3_row, Stmp2.row(m), levels[sonLevel].lfi2D);
					Stmp3.row(m) = Stmp3_row;
				}
				levels[sonLevel].Sdown_obs(sonIndex) += Stmp3;
			}
		}
		thisLevel--;
	}
	// and finally the integration
	E_obs = ArrayXXcd::Zero(E_obs.rows(), E_obs.cols());
	const int N_cubes_obs = levels[thisLevel].cubes_obs.size();
	for (int i = 0; i < N_cubes_obs; ++i)
	{
		levels[thisLevel].sphericalIntegration_obs
		(E_obs, levels[thisLevel].Sdown_obs(i), levels[thisLevel].cubes_obs[i],
		 levels[thisLevel].thetas, levels[thisLevel].phis, w, mu_r, k, CFIE,
		 M_current);
	}
}

void Octtree::ZINearComputation_transpose(ArrayXXcd& E_obs,
                                          ArrayXcd& I_PQ, const bool M_current)
{
	for (int l = 0; l < N_levels; ++l)
	{
		const int N_cubes = levels[l].cubes.size();
		const int N_theta = levels[l].thetas.size();
		const int N_phi = levels[l].phis.size();
		const int N_directions = N_theta * N_phi;
		levels[l].Sdown.resize(N_cubes);
		for (int i = 0; i < N_cubes; ++i)
		{
			levels[l].Sdown(i).resize(2, N_directions);
		}
	}

	updateSup_obs(E_obs, true, M_current);
	alphaTranslations_obs_transpose();
	const int N_levels = levels.size(), L = N_levels - 1;
	int thisLevel = L;

	// now the descent towards the bottom
	while (thisLevel > 0)
	{
		const int sonLevel = thisLevel - 1;
		const int N_cubes = levels[thisLevel].cubes.size();

		ArrayXXcd Stmp2(2, levels[thisLevel].thetas.size() * levels[thisLevel].phis.size());
		ArrayXXcd Stmp3(2, levels[sonLevel].thetas.size() * levels[sonLevel].phis.size());

		for (int i = 0; i < N_cubes; i++)
		{
			//			if (thisLevel == L && i == 0) cout << Stmp2 << endl;
			vector<int> sonsIndexes = levels[thisLevel].cubes[i].sonsIndexes;

			for (unsigned int j = 0; j < sonsIndexes.size(); ++j)
			{
				const int sonIndex = sonsIndexes[j];
				// shifting
				const double* rc_1(levels[sonLevel].cubes[sonIndex].rCenter);
				const double* rc_2(levels[thisLevel].cubes[i].rCenter);
				double DRcenters[3] = {rc_1[0] - rc_2[0],
					rc_1[1] - rc_2[1], rc_1[2] - rc_2[2]
				};
				Stmp2 = levels[thisLevel].Sdown(i);
				shiftExp(Stmp2, levels[thisLevel].getShiftingArray(DRcenters[0],
				                                                   DRcenters[1],
				                                                   DRcenters[2]));
				// anterpolate
				ArrayXcd Stmp3_row(Stmp3.cols());
				for (int m = 0; m < 2; m++)
				{
					anterpolate2Dlfi(Stmp3_row, Stmp2.row(m), levels[sonLevel].lfi2D);
					Stmp3.row(m) = Stmp3_row;
				}

				levels[sonLevel].Sdown(sonIndex) += Stmp3;
			}
		}
		thisLevel--;
	}
	// and finally the integration
	E_obs = ArrayXXcd::Zero(E_obs.rows(), E_obs.cols());
	const int N_cubes = levels[thisLevel].cubes.size();
	for (int i = 0; i < N_cubes; i++)
	{
		levels[thisLevel].sphericalIntegration
		(I_PQ, levels[thisLevel].Sdown(i), levels[thisLevel].cubes[i],
		 levels[thisLevel].thetas, levels[thisLevel].phis, w, mu_r, k,
		 CFIE, M_current, target_surface);
	}
}

void Octtree::computeSupLastLevel(ArrayXXcd& SupLastLevel,
                                  const ArrayXd& r_phase_center,
                                  const ArrayXd& octtreeXthetas_coarsest,
                                  const ArrayXd& octtreeXphis_coarsest,
                                  const ArrayXcd& I_PQ,
                                  const string octtree_data_path)
{
	const int N_levels = levels.size();
	int stopLevel;
	for (int l = 0; l < N_levels; ++l)
	{
		// cout << levels[l].getLevel() << "..."; flush(cout);
		const int N_cubes = levels[l].getLevelSize();
		const int N_theta = levels[l].thetas.size(), N_phi = levels[l].phis.size();
		const int N_directions = N_theta * N_phi;

		if (levels[l].Sdown.size() == 0) levels[l].Sdown.resize(N_cubes);
		// we first compute the Sups of all the cubes at the given level
		if (l == 0)
		{ // we use computeSup only for the leaf cubes
			for (int i = 0; i < N_cubes; ++i)
			{
				if (levels[l].Sdown(i).size() == 0)
					levels[l].Sdown(i).resize(2, N_theta * N_phi);
				levels[l].computeSup(levels[l].Sdown(i), k, I_PQ,
				                     levels[l].cubes[i], levels[l].thetas,
				                     levels[l].phis, target_surface);
			}
		}
		else if (l > 0)
		{
			// the Sups are obtained by interpolation
			// and shifting from the sonsCubes Sups
			ArrayXXcd S_tmp(2, N_theta * N_phi);
			ArrayXXcd S_tmp2(2, N_theta * N_phi);
			for (int i = 0; i < N_cubes; ++i)
			{
				if (levels[l].Sdown(i).size() == 0)
					levels[l].Sdown(i).resize(2, N_directions);
				S_tmp2 = ArrayXXcd::Zero(2, N_theta * N_phi);
				const int NSons = levels[l].cubes[i].sonsIndexes.size();
				for (int j = 0; j < NSons; ++j)
				{
					const int sonIndex = levels[l].cubes[i].sonsIndexes[j];
					// interpolation
					ArrayXcd S_tmp_row(N_theta * N_phi);
					interpolate2Dlfi(S_tmp_row,
					                 levels[l - 1].Sdown(sonIndex).row(0),
					                 levels[l - 1].lfi2D);
					S_tmp.row(0) = S_tmp_row;
					interpolate2Dlfi(S_tmp_row,
					                 levels[l - 1].Sdown(sonIndex).row(1),
					                 levels[l - 1].lfi2D);
					S_tmp.row(1) = S_tmp_row;
					// shifting
					const double* rc_1(levels[l].cubes[i].rCenter);
					const double* rc_2(levels[l - 1].cubes[sonIndex].rCenter);
					const double DRcenters[3] = {rc_1[0] - rc_2[0],
						rc_1[1] - rc_2[1], rc_1[2] - rc_2[2]
					};

					shiftExp(S_tmp, levels[l].getShiftingArray(DRcenters[0],
					                                           DRcenters[1],
					                                           DRcenters[2]));
					S_tmp2 += S_tmp;
				}
				levels[l].Sdown(i) = S_tmp2;
			}
		}
		stopLevel = l;
	}

	// we now define an interpolator from the last level
	// we went through to the coarsest level.
	const int N_thetaCoarseLevel(octtreeXthetas_coarsest.size());
	const int N_phiCoarseLevel(octtreeXphis_coarsest.size());
	const int N_directions = N_thetaCoarseLevel * N_phiCoarseLevel;
	int INCLUDED_THETA_BOUNDARIES, PERIODIC_Theta;
	int CYCLIC_Theta, NOrderInterpTheta;

	readIntFromASCIIFile(octtree_data_path + "INCLUDED_THETA_BOUNDARIES.txt",
	                     INCLUDED_THETA_BOUNDARIES);
	readIntFromASCIIFile(octtree_data_path + "PERIODIC_Theta.txt", PERIODIC_Theta);
	readIntFromASCIIFile(octtree_data_path + "CYCLIC_Theta.txt", CYCLIC_Theta);
	readIntFromASCIIFile(octtree_data_path + "NOrderInterpTheta.txt",
	                     NOrderInterpTheta);

	// phi data
	int INCLUDED_PHI_BOUNDARIES, PERIODIC_Phi, CYCLIC_Phi, NOrderInterpPhi;
	readIntFromASCIIFile(octtree_data_path + "INCLUDED_PHI_BOUNDARIES.txt",
	                     INCLUDED_PHI_BOUNDARIES);
	readIntFromASCIIFile(octtree_data_path + "PERIODIC_Phi.txt", PERIODIC_Phi);
	readIntFromASCIIFile(octtree_data_path + "CYCLIC_Phi.txt", CYCLIC_Phi);
	readIntFromASCIIFile(octtree_data_path + "NOrderInterpPhi.txt",
	                     NOrderInterpPhi);

	const double THETA_MIN = 0., THETA_MAX = M_PI;
	const double PHI_MIN = 0., PHI_MAX = 2.0 * M_PI;
	LagrangeFastInterpolator2D FarFieldInterpolator(octtreeXthetas_coarsest, levels[stopLevel].thetas,
	                                                THETA_MIN, THETA_MAX, INCLUDED_THETA_BOUNDARIES,
	                                                NOrderInterpTheta, PERIODIC_Theta, CYCLIC_Theta,
	                                                octtreeXphis_coarsest, levels[stopLevel].phis,
	                                                PHI_MIN, PHI_MAX, INCLUDED_PHI_BOUNDARIES,
	                                                NOrderInterpPhi, PERIODIC_Phi, CYCLIC_Phi);

	// actual radiation pattern computation
	ArrayXXcd SupLastLevelTmp = ArrayXXcd::Zero(2, N_directions);
	const int N_cubes = levels[stopLevel].getLevelSize();

	for (int i = 0; i < N_cubes; ++i)
	{
		SupLastLevelTmp = ArrayXXcd::Zero(2, N_directions);
		ArrayXcd Sup_tmp(N_directions);
		interpolate2Dlfi(Sup_tmp,
		                 levels[stopLevel].Sdown(i).row(0), FarFieldInterpolator);
		SupLastLevelTmp.row(0) = Sup_tmp;
		interpolate2Dlfi(Sup_tmp,
		                 levels[stopLevel].Sdown(i).row(1), FarFieldInterpolator);
		SupLastLevelTmp.row(1) = Sup_tmp;

		// shifting
		double DRcenters[3];
		for (int j = 0; j < 3; j++)
			DRcenters[j] = r_phase_center(j) - levels[stopLevel].cubes[i].rCenter[j];

		ArrayXcd shiftingArray(N_directions);

		for (int m = 0; m < N_thetaCoarseLevel; ++m)
		{
			const double sinTheta = sin(octtreeXthetas_coarsest(m));
			const double cosTheta = cos(octtreeXthetas_coarsest(m));
			for (int n = 0; n < N_phiCoarseLevel; ++n)
			{
				const double k_hat[3] = {sinTheta * cos(octtreeXphis_coarsest(n)),
					sinTheta * sin(octtreeXphis_coarsest(n)),
					cosTheta
				};

				shiftingArray(m + n * N_thetaCoarseLevel)
					= exp(-I * this->k * (k_hat[0] * DRcenters[0]
						+ k_hat[1] * DRcenters[1]
						+ k_hat[2] * DRcenters[2]));
			}
		}
		shiftExp(SupLastLevelTmp, shiftingArray);
		SupLastLevel += SupLastLevelTmp / (4.0 * M_PI);
	}
}

void Octtree::computeFarField(complex<double> e_theta_far_tmp[],
                              complex<double> e_phi_far_tmp[],
                              double r_phase_center_tmp[],
                              double octtreeXthetas_coarsest_tmp[],
                              double octtreeXphis_coarsest_tmp[],
                              int octtreeNthetas_coarsest,
                              int octtreeNphis_coarsest,
                              complex<double> I_J_tmp[],
                              complex<double> I_M_tmp[],
                              int N_RWG,
                              string octtree_data_path)
{
	typedef Array<complex<double>, Dynamic, Dynamic, RowMajor> ArrayXXcd_RowMajor;

	ArrayXd r_phase_center = Map<ArrayXd>(r_phase_center_tmp, 3);
	ArrayXd octtreeXthetas_coarsest = Map<ArrayXd>(octtreeXthetas_coarsest_tmp, octtreeNthetas_coarsest);
	ArrayXd octtreeXphis_coarsest = Map<ArrayXd>(octtreeXphis_coarsest_tmp, octtreeNphis_coarsest);
	ArrayXcd I_J = Map<ArrayXcd>(I_J_tmp, N_RWG), I_M = Map<ArrayXcd>(I_M_tmp, N_RWG);

	cout << "\nFar field computation" << ": level ";
	const complex<double> EJ_factor(-I * mu_0 * w * mu_r);
	const complex<double> EM_factor(-I * k);

	int N_directions = octtreeNthetas_coarsest * octtreeNphis_coarsest;
	ArrayXXcd SupLastLevel_J = ArrayXXcd::Zero(2, N_directions);
	ArrayXXcd SupLastLevel_M = ArrayXXcd::Zero(2, N_directions);

	computeSupLastLevel(SupLastLevel_J, r_phase_center, octtreeXthetas_coarsest,
	                    octtreeXphis_coarsest, I_J, octtree_data_path);
	computeSupLastLevel(SupLastLevel_M, r_phase_center, octtreeXthetas_coarsest,
	                    octtreeXphis_coarsest, I_M, octtree_data_path);

	const int N_thetaCoarseLevel(octtreeXthetas_coarsest.size());
	const int N_phiCoarseLevel(octtreeXphis_coarsest.size());
	// we now gather all
	ArrayXXcd e_theta_far(N_thetaCoarseLevel, N_phiCoarseLevel);
	ArrayXXcd e_phi_far(N_thetaCoarseLevel, N_phiCoarseLevel);
	// e_theta
	int index = 0;
	for (int m = 0; m < N_thetaCoarseLevel; ++m)
	{
		for (int n = 0; n < N_phiCoarseLevel; ++n)
		{
			e_theta_far_tmp[index]
				= EJ_factor * SupLastLevel_J(0, m + n * N_thetaCoarseLevel)
				+ EM_factor * SupLastLevel_M(1, m + n * N_thetaCoarseLevel);
			index++;
		}
	}
	// e_phi
	index = 0;
	for (int m = 0; m < N_thetaCoarseLevel; ++m)
	{
		for (int n = 0; n < N_phiCoarseLevel; ++n)
		{
			e_phi_far_tmp[index]
				= EJ_factor * SupLastLevel_J(1, m + n * N_thetaCoarseLevel)
				- EM_factor * SupLastLevel_M(0, m + n * N_thetaCoarseLevel);
			index++;
		}
	}
	cout << "finished!" << endl;
	flush(cout);
}

// void Octtree::computeSourceFarField
// (Array<complex<float>, 2>& e_theta_far,
//  Array<complex<float>, 2>& e_phi_far,
//  const Array<float, 1>& r_phase_center,
//  const Array<float, 1>& octtreeXthetas_coarsest,
//  const Array<float, 1>& octtreeXphis_coarsest,
//  const int J_DIPOLES_EXCITATION,
//  const Array<complex<float>, 2>& J_dip,
//  const Array<float, 2>& r_J_dip,
//  const int M_DIPOLES_EXCITATION,
//  const Array<complex<float>, 2>& M_dip,
//  const Array<float, 2>& r_M_dip)
// {
// 	cout << "\nSource Far field computation...";
// 	const int N_thetas = octtreeXthetas_coarsest.size(), N_phis = octtreeXphis_coarsest.size();

// 	Array<complex<float>, 2> SupTmp(2, N_thetas * N_phis), Sup(2, N_thetas * N_phis);
// 	Sup = 0.0;
// 	if (J_DIPOLES_EXCITATION == 1) {
// 		const int N_J_dipoles = J_dip.rows();
// 		for (int i = 0; i < N_J_dipoles; i++) {
// 			const int IS_J_CURRENT = 1;
// 			const complex<float> J[3] = {J_dip(i, 0), J_dip(i, 1), J_dip(i, 2)};
// 			const float r_dip[3] = {r_J_dip(i, 0), r_J_dip(i, 1), r_J_dip(i, 2)};
// 			const float rCenter[3] = {r_phase_center(0), r_phase_center(1), r_phase_center(2)};
// 			computeDipoleSup(SupTmp, J, IS_J_CURRENT, r_dip, rCenter, octtreeXthetas_coarsest, octtreeXphis_coarsest);
// 			Sup += (static_cast<complex<float> >(-I * mu_0) * w / static_cast<float>(4.0 * M_PI) * mu_r) * SupTmp;
// 		}
// 	}
// 	if (M_DIPOLES_EXCITATION == 1) {
// 		const int N_M_dipoles = M_dip.rows();
// 		for (int i = 0; i < N_M_dipoles; i++) {
// 			const int IS_J_CURRENT = 0;
// 			const complex<float> M[3] = {M_dip(i, 0), M_dip(i, 1), M_dip(i, 2)};
// 			const float r_dip[3] = {r_M_dip(i, 0), r_M_dip(i, 1), r_M_dip(i, 2)};
// 			const float rCenter[3] = {r_phase_center(0), r_phase_center(1), r_phase_center(2)};
// 			computeDipoleSup(SupTmp, M, IS_J_CURRENT, r_dip, rCenter,
// 			                 octtreeXthetas_coarsest, octtreeXphis_coarsest);
// 			Sup += static_cast<complex<float> >((I * k) / (4.0 * M_PI)) * SupTmp;
// 		}
// 	}
// 	e_theta_far.resize(N_thetas, N_phis);
// 	e_phi_far.resize(N_thetas, N_phis);
// 	for (int m = 0 ; m < N_thetas ; ++m) {
// 		for (int n = 0 ; n < N_phis ; ++n) {
// 			e_theta_far(m, n) = Sup(0, m + n * N_thetas);
// 			e_phi_far(m, n) = Sup(1, m + n * N_thetas);
// 		}
// 	}
// 	cout << "finished!" << endl;
// 	flush(cout);
// }

// void Octtree::computeDipoleSup(Array<complex<float>, 2> & Sup,
//                                const complex<float> J_dipole[3],
//                                const int IS_J_CURRENT,
//                                const float r_dipole[3],
//                                const float rCenter[3],
//                                const Array<float, 1>& thetas,
//                                const Array<float, 1>& phis)
// {
// 	const int NThetas = thetas.size(), NPhis = phis.size();

// 	Array< float, 2> kHats(NThetas * NPhis, 3),
// 	       thetaHats(NThetas * NPhis, 3), phiHats(NThetas * NPhis, 3);
// 	Array< complex<float>, 2> FC3Components(NThetas * NPhis, 3);
// 	vector<float> sin_thetas, cos_thetas;
// 	sin_thetas.resize(NThetas);
// 	cos_thetas.resize(NThetas);
// 	for (int p = 0 ; p < NThetas ; ++p) sincosf(thetas(p), &sin_thetas[p], &cos_thetas[p]);
// 	// initialisation of arrays
// 	for (int q = 0 ; q < NPhis ; ++q) {
// 		const float cos_phi = cos(phis(q)), sin_phi = sin(phis(q));
// 		for (int p = 0 ; p < NThetas ; ++p) {
// 			int index = p + q * NThetas;
// 			const float sin_theta = sin_thetas[p], cos_theta = cos_thetas[p];
// 			kHats(index, 0) = sin_theta * cos_phi;
// 			kHats(index, 1) = sin_theta * sin_phi;
// 			kHats(index, 2) = cos_theta;
// 			thetaHats(index, 0) = cos_theta * cos_phi;
// 			thetaHats(index, 1) = cos_theta * sin_phi;
// 			thetaHats(index, 2) = -sin_theta;
// 			phiHats(index, 0) = -sin_phi;
// 			phiHats(index, 1) = cos_phi;
// 			phiHats(index, 2) = 0.0;
// 			for (int i = 0 ; i < 3 ; i++) FC3Components(index, i) = 0.0;
// 		}
// 	}
// 	// computation of FC3Components array
// 	const complex<float> I_k(static_cast<complex<float> >(I * this->k));
// 	const complex<float> fj[3] = {J_dipole[0], J_dipole[1], J_dipole[2]};
// 	const float expArg[3] = {r_dipole[0] - rCenter[0], r_dipole[1] - rCenter[1], r_dipole[2] - rCenter[2]};
// 	for (int q = 0 ; q < NPhis ; q++) {
// 		const int index_1 = q * NThetas;
// 		for (int p = 0 ; p < NThetas ; p++) {
// 			const int index = p + index_1;
// 			const complex<float>
// 			a(I_k * (expArg[0]*kHats(index, 0)
// 			         + expArg[1]*kHats(index, 1) + expArg[2]*kHats(index, 2)));
// 			float c, s, e;
// 			e = (a.real() == 0.0) ? 1.0 : exp(a.real());
// 			sincosf(a.imag(), &s, &c);
// 			const complex<float> EXP(e * c, e * s);
// 			if (IS_J_CURRENT == 1) { // we have electric dipole
// 				FC3Components(index, 0) += fj[0] * EXP;
// 				FC3Components(index, 1) += fj[1] * EXP;
// 				FC3Components(index, 2) += fj[2] * EXP;
// 			} else { // we have magnetic dipole
// 				FC3Components(index, 0) += (kHats(index, 1) * fj[2] - kHats(index, 2) * fj[1]) * EXP;
// 				FC3Components(index, 1) += (kHats(index, 2) * fj[0] - kHats(index, 0) * fj[2]) * EXP;
// 				FC3Components(index, 2) += (kHats(index, 0) * fj[1] - kHats(index, 1) * fj[0]) * EXP;
// 			}
// 		}
// 	} // end q loop
// 	// transformation from cartesian to spherical coordinates and assignation to Sup
// 	for (int i = 0 ; i < Sup.extent(1) ; ++i) {
// 		Sup(0, i) = thetaHats(i, 0) * FC3Components(i, 0) + thetaHats(i, 1) * FC3Components(i, 1) + thetaHats(i, 2) * FC3Components(i, 2);
// 		Sup(1, i) = phiHats(i, 0) * FC3Components(i, 0) + phiHats(i, 1) * FC3Components(i, 1) + phiHats(i, 2) * FC3Components(i, 2);
// 	}
// }
