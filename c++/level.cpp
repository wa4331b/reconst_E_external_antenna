/**********************************************************************
 *
 * level.cpp
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
#include <algorithm>

using namespace std;

#include "level.hpp"
#include "readWriteBlitzArrayFromFile.hpp"
#include "alpha_computation.hpp"


/****************************************************************************/
/********************************* Level ************************************/
/****************************************************************************/
Level::Level()
{
};

Level::Level(const int l,
             const double leaf_side_length,
             const double big_cube_lower_coord[3],
             const ArrayXXd& cubes_centroids)
{
	numberTimesCopied = 0;
	level = l;
	cubeSideLength = leaf_side_length;
	maxNumberCubes1D = static_cast<int>(pow(2.0, level));
	int N_cubes_level_L = cubes_centroids.rows();
	cubes.reserve(N_cubes_level_L);
	for (int j = 0; j < N_cubes_level_L; ++j)
	{
		const double r_c[3] = {cubes_centroids(j, 0),
			cubes_centroids(j, 1), cubes_centroids(j, 2)
		};
		addNode(Cube(level, leaf_side_length, big_cube_lower_coord, r_c));
		cubes[j].setOldIndex(j); // the original index, from the python mesh
	}
}

Level::Level(const int l,
             const int N_expansion,
             const double leaf_side_length,
             const double big_cube_lower_coord[3],
             const ArrayXXd& cubes_centroids,
             const ArrayXXd& cubes_centroids_obs,
             const complex<double>& waveNumber,
             const ArrayXd& Xthetas,
             const ArrayXd& Wthetas,
             const int INCLUDED_THETA_BOUNDARIES,
             const int N_theta,
             const int PERIODIC_Theta,
             const int CYCLIC_Theta,
             const int NOrderInterpolatorTheta,
             const ArrayXd& Xphis,
             const ArrayXd& Wphis,
             const int INCLUDED_PHI_BOUNDARIES,
             const int N_phi,
             const int PERIODIC_Phi,
             const int CYCLIC_Phi,
             const int NOrderInterpolatorPhi,
             const ArrayXd& XthetasNextLevel,
             const ArrayXd& XphisNextLevel,
             const int VERBOSE) // leaf level constructor: no sons needed
{
	numberTimesCopied = 0;
	k = waveNumber;
	level = l;
	leaf = true;
	ceiling = 0;
	N = N_expansion;
	cubeSideLength = leaf_side_length;
	maxNumberCubes1D = static_cast<int>(pow(2.0, level));

	if (VERBOSE == 1) cout << "construction of the leaf (finest) level" << endl;
	int j, N_cubes_level_L = cubes_centroids.rows();
	int N_cubes_obs_level_L = cubes_centroids_obs.rows();

	cubes.resize(N_cubes_level_L);
	for (j = 0; j < N_cubes_level_L; ++j)
	{
		const double r_c[3] = {cubes_centroids(j, 0), cubes_centroids(j, 1),
			cubes_centroids(j, 2)
		};
		cubes[j] = Cube(level, leaf_side_length, big_cube_lower_coord, r_c);
		cubes[j].setOldIndex(j); // the original index, from the python mesh
	}

	cubes_obs.resize(N_cubes_obs_level_L);
	for (j = 0; j < N_cubes_obs_level_L; ++j)
	{
		const double r_c[3] = {cubes_centroids_obs(j, 0),
			cubes_centroids_obs(j, 1),
			cubes_centroids_obs(j, 2)
		};
		cubes_obs[j] = Cube(level, leaf_side_length, big_cube_lower_coord, r_c);
		cubes_obs[j].setOldIndex(j); // the original index, from the python mesh
	}

	if (VERBOSE == 1)
		cout << "cubes.size() = " << cubes.size()
			<< ", cubes.capacity() = " << cubes.capacity() << endl;
	if (VERBOSE == 1)
		cout << "cubes_obs.size() = " << cubes_obs.size()
			<< ", cubes_obs.capacity() = " << cubes_obs.capacity() << endl;

	thetas.resize(N_theta);
	thetas = Xthetas;
	weightsThetas.resize(N_theta);
	weightsThetas = Wthetas;
	phis.resize(N_phi);
	phis = Xphis;
	weightsPhis.resize(N_phi);
	weightsPhis = Wphis;

	// Lagrange Fast Interpolator
	const double THETA_MIN = 0., THETA_MAX = M_PI;
	const double PHI_MIN = 0., PHI_MAX = 2.0 * M_PI;
	lfi2D.setLfi2D(LagrangeFastInterpolator2D(XthetasNextLevel, Xthetas,
	                                          THETA_MIN, THETA_MAX,
	                                          INCLUDED_THETA_BOUNDARIES,
	                                          NOrderInterpolatorTheta, PERIODIC_Theta,
	                                          CYCLIC_Theta,
	                                          XphisNextLevel, Xphis,
	                                          PHI_MIN, PHI_MAX,
	                                          INCLUDED_PHI_BOUNDARIES,
	                                          NOrderInterpolatorPhi, PERIODIC_Phi,
	                                          CYCLIC_Phi));

	cout << "The total number of directions at level " << getLevel();
	cout << " is N_theta*N_phi = " << N_theta;
	cout << "*" << N_phi << " = " << N_theta * N_phi << endl << endl;
}

void Level::copyLevel(const Level& levelToCopy) // copy constructor
{
	numberTimesCopied = levelToCopy.getNumberTimesCopied();
	incrementNumberTimesCopied();
	leaf = levelToCopy.getLeaf();
	ceiling = levelToCopy.getCeiling();
	level = levelToCopy.getLevel();

	if (numberTimesCopied > 15)
		cout << "    Warning: level " << level
			<< "has been copied more than 15 times. numberTimesCopied = "
			<< numberTimesCopied << endl;

	N = levelToCopy.getN();
	cubeSideLength = levelToCopy.getCubeSideLength();
	maxNumberCubes1D = levelToCopy.getMaxNumberCubes1D();
	NCubesX = levelToCopy.getNCubesX();
	NCubesY = levelToCopy.getNCubesY();
	NCubesZ = levelToCopy.getNCubesZ();
	offsetAlphaIndexX = levelToCopy.getOffsetAlphaIndexX();
	offsetAlphaIndexY = levelToCopy.getOffsetAlphaIndexY();
	offsetAlphaIndexZ = levelToCopy.getOffsetAlphaIndexZ();
	k = levelToCopy.getK();
	cubes.resize(levelToCopy.cubes.size());

	for (unsigned int i = 0; i < levelToCopy.cubes.size(); ++i)
		cubes[i] = levelToCopy.cubes[i];
	cubes_obs.resize(levelToCopy.cubes_obs.size());
	for (unsigned int i = 0; i < levelToCopy.cubes_obs.size(); ++i)
		cubes_obs[i] = levelToCopy.cubes_obs[i];
	numbersToIndexes = levelToCopy.getNumbersToIndexes();

	const int N_theta = levelToCopy.getNThetas();
	const int N_phi = levelToCopy.getNPhis();
	thetas.resize(N_theta);
	phis.resize(N_phi);
	thetas = levelToCopy.getThetas();
	phis = levelToCopy.getPhis();
	Array3i alphaTranslationsExtents;
	alphaTranslationsExtents = levelToCopy.getAlphaTranslationsExtents();
	const int Nx = alphaTranslationsExtents(0);
	const int Ny = alphaTranslationsExtents(1), Nz = alphaTranslationsExtents(2);

	alphaTranslations = levelToCopy.alphaTranslations;
	alphaTranslationsIndexesNonZeros = levelToCopy.alphaTranslationsIndexesNonZeros;
	alphaTranslationsIndexes = levelToCopy.alphaTranslationsIndexes;

	shiftingArrays.resize(levelToCopy.getShiftingArrays().rows(),
	                      levelToCopy.getShiftingArrays().cols());
	shiftingArrays = levelToCopy.getShiftingArrays();
	weightsThetas = levelToCopy.weightsThetas;
	weightsPhis = levelToCopy.weightsPhis;
	lfi2D = levelToCopy.lfi2D;
	Sdown = levelToCopy.Sdown;
	Sdown_obs = levelToCopy.Sdown_obs;
}


Level::Level(const Level& levelToCopy) // copy constructor
{
	copyLevel(levelToCopy);
}

Level& Level::operator=(const Level& levelToCopy) // copy assignment
{
	copyLevel(levelToCopy);
	return *this;
}

Level::Level(const Level& sonLevel,
             const int N_expansion,
             const double big_cube_lower_coord[3],
             const ArrayXd& Xthetas,
             const ArrayXd& Wthetas,
             const int INCLUDED_THETA_BOUNDARIES,
             const int N_theta,
             const int PERIODIC_Theta,
             const int CYCLIC_Theta,
             const int NOrderInterpolatorTheta,
             const ArrayXd& Xphis,
             const ArrayXd& Wphis,
             const int INCLUDED_PHI_BOUNDARIES,
             const int N_phi,
             const int PERIODIC_Phi,
             const int CYCLIC_Phi,
             const int NOrderInterpolatorPhi,
             const ArrayXd& XthetasNextLevel,
             const ArrayXd& XphisNextLevel,
             const int VERBOSE)
{
	numberTimesCopied = 0;
	N = N_expansion;
	k = sonLevel.getK();
	int j;
	level = sonLevel.getLevel() - 1;
	leaf = false; // because this level is created from a lower (finer) level...
	ceiling = 0;
	cubeSideLength = 2.0 * sonLevel.getCubeSideLength();
	maxNumberCubes1D = sonLevel.getMaxNumberCubes1D() / 2;
	if (VERBOSE == 1) cout << "construction of level " << level << endl;
	flush(cout);
	addNode(Cube(sonLevel.getCube(0), level, big_cube_lower_coord, cubeSideLength)); // initialization
	for (j = 1; j < sonLevel.getLevelSize(); ++j)
	{ // we walk through the sons list
		if (sonLevel.getCube(j).getFatherNumber() == cubes.back().getNumber())
			cubes.back().addSon(sonLevel.getCube(j));
		else
			addNode(Cube(sonLevel.getCube(j), level, big_cube_lower_coord, cubeSideLength));
	}
	cubes_obs.push_back(Cube(sonLevel.cubes_obs[0], level, big_cube_lower_coord, cubeSideLength)); // initialization
	for (j = 1; j < int(sonLevel.cubes_obs.size()); ++j)
	// we walk through the sons list
	{
		if (sonLevel.cubes_obs[j].getFatherNumber() == cubes_obs.back().getNumber())
			cubes_obs.back().addSon(sonLevel.cubes_obs[j]);
		else
			cubes_obs.push_back(Cube(sonLevel.cubes_obs[j], level, big_cube_lower_coord, cubeSideLength));
	}
	vector<Cube>(cubes).swap(cubes);
	// swap trick for trimming exceeding capacity
	vector<Cube>(cubes_obs).swap(cubes_obs);
	// swap trick for trimming exceeding capacity
	if (VERBOSE == 1)
	{
		cout << "cubes.size() = " << cubes.size();
		cout << ", cubes.capacity() = " << cubes.capacity() << endl;
		cout << "cubes_obs.size() = " << cubes_obs.size();
		cout << ", cubes_obs.capacity() = " << cubes_obs.capacity() << endl;
	}
	thetas.resize(N_theta);
	thetas = Xthetas;
	weightsThetas.resize(N_theta);
	weightsThetas = Wthetas;
	phis.resize(N_phi);
	phis = Xphis;
	weightsPhis.resize(N_phi);
	weightsPhis = Wphis;
	// Lagrange Fast Interpolator. The coarsest level does not need one...
	const double THETA_MIN = 0., THETA_MAX = M_PI;
	const double PHI_MIN = 0., PHI_MAX = 2.0 * M_PI;
	if (level > 2)
		lfi2D.setLfi2D(LagrangeFastInterpolator2D(
			XthetasNextLevel, Xthetas, THETA_MIN, THETA_MAX, INCLUDED_THETA_BOUNDARIES,
			NOrderInterpolatorTheta, PERIODIC_Theta, CYCLIC_Theta,
			XphisNextLevel, Xphis, PHI_MIN, PHI_MAX, INCLUDED_PHI_BOUNDARIES,
			NOrderInterpolatorPhi, PERIODIC_Phi, CYCLIC_Phi
		));

	cout << "The total number of directions at level " << getLevel();
	cout << " is N_theta*N_phi = " << N_theta << "*" << N_phi << " = " << N_theta * N_phi << endl << endl;
}

Level::~Level()
{
	cubes.clear();
	cubes_obs.clear();
	numbersToIndexes.clear();
	listOfFcToBeReceived.clear();
	listOfFcToBeSent.clear();
	//	thetas.free();
	//	phis.free();
	//	weightsThetas.free();
	//	weightsPhis.free();
	//	alphaTranslations.free();
	//	alphaTranslationsIndexesNonZeros.free();
	//	alphaTranslationsIndexes.free();
	//	shiftingArrays.free();
	//	Sdown.free();
	//	Sdown_obs.free();
}

void Level::NCubesXYZComputation()
{
	int NxMax, NyMax, NzMax, NxMin, NyMin, NzMin;
	NxMax = NyMax = NzMax = 0;
	NxMin = NyMin = NzMin = maxNumberCubes1D;
	for (unsigned int j = 0; j < cubes.size(); ++j)
	{
		const double* absoluteCartesianCoordTmp(cubes[j].absoluteCartesianCoord);
		if (static_cast<int>(absoluteCartesianCoordTmp[0]) > NxMax)
			NxMax = static_cast<int>(absoluteCartesianCoordTmp[0]);
		if (static_cast<int>(absoluteCartesianCoordTmp[1]) > NyMax)
			NyMax = static_cast<int>(absoluteCartesianCoordTmp[1]);
		if (static_cast<int>(absoluteCartesianCoordTmp[2]) > NzMax)
			NzMax = static_cast<int>(absoluteCartesianCoordTmp[2]);
		if (static_cast<int>(absoluteCartesianCoordTmp[0]) < NxMin)
			NxMin = static_cast<int>(absoluteCartesianCoordTmp[0]);
		if (static_cast<int>(absoluteCartesianCoordTmp[1]) < NyMin)
			NyMin = static_cast<int>(absoluteCartesianCoordTmp[1]);
		if (static_cast<int>(absoluteCartesianCoordTmp[2]) < NzMin)
			NzMin = static_cast<int>(absoluteCartesianCoordTmp[2]);
	}
	for (unsigned int j = 0; j < cubes_obs.size(); ++j)
	{
		const double* absoluteCartesianCoordTmp(cubes_obs[j].absoluteCartesianCoord);
		if (static_cast<int>(absoluteCartesianCoordTmp[0]) > NxMax)
			NxMax = static_cast<int>(absoluteCartesianCoordTmp[0]);
		if (static_cast<int>(absoluteCartesianCoordTmp[1]) > NyMax)
			NyMax = static_cast<int>(absoluteCartesianCoordTmp[1]);
		if (static_cast<int>(absoluteCartesianCoordTmp[2]) > NzMax)
			NzMax = static_cast<int>(absoluteCartesianCoordTmp[2]);
		if (static_cast<int>(absoluteCartesianCoordTmp[0]) < NxMin)
			NxMin = static_cast<int>(absoluteCartesianCoordTmp[0]);
		if (static_cast<int>(absoluteCartesianCoordTmp[1]) < NyMin)
			NyMin = static_cast<int>(absoluteCartesianCoordTmp[1]);
		if (static_cast<int>(absoluteCartesianCoordTmp[2]) < NzMin)
			NzMin = static_cast<int>(absoluteCartesianCoordTmp[2]);
	}
	NCubesX = NxMax - NxMin + 1;
	NCubesY = NyMax - NyMin + 1;
	NCubesZ = NzMax - NzMin + 1;
	cout << "Level " << this->level << " : NCubesX, NCubesY, NCubesZ = "
		<< NCubesX << ", " << NCubesY << ", " << NCubesZ << endl;
	offsetAlphaIndexX = 0;
	offsetAlphaIndexY = 0;
	offsetAlphaIndexZ = 0;
	cout << "Level " << this->level
		<< " : offsetAlphaIndexX, offsetAlphaIndexY, offsetAlphaIndexZ = "
		<< offsetAlphaIndexX << ", " << offsetAlphaIndexY << ", "
		<< offsetAlphaIndexZ << endl;
}

void Level::alphaTranslationsComputation(const double alphaTranslation_smoothing_factor,
                                         const double alphaTranslation_thresholdRelValueMax,
                                         const double alphaTranslation_RelativeCountAboveThreshold)
{
	//	cout << alphaTranslation_smoothing_factor << " "
	//		<< alphaTranslation_thresholdRelValueMax << " "
	//		<< alphaTranslation_RelativeCountAboveThreshold << endl;

	const int translationOrder = getN(), NThetas = getNThetas(), NPhis = getNPhis();
	const int translationOrder_prime = static_cast<int>(ceil(translationOrder * alphaTranslation_smoothing_factor));

	int Nx = (this->getCeiling()) ? NCubesX : min(NCubesX, 4);
	int Ny = (this->getCeiling()) ? NCubesY : min(NCubesY, 4);
	int Nz = (this->getCeiling()) ? NCubesZ : min(NCubesZ, 4);

	ArrayXXd thetasPhis_all_directions(NThetas * NPhis, 2), thetasPhis;
	ArrayXd weightsThetasPhis_all_directions(NThetas * NPhis), weightsThetasPhis;
	for (int i = 0; i < NThetas; ++i)
	{
		for (int j = 0; j < NPhis; ++j)
		{
			thetasPhis_all_directions(i + j * NThetas, 0) = thetas(i);
			thetasPhis_all_directions(i + j * NThetas, 1) = phis(j);
			weightsThetasPhis_all_directions(i + j * NThetas) = weightsThetas(i) * weightsPhis(j);
		}
	}

	int N_directions;

	N_directions = NThetas * NPhis;
	thetasPhis.resize(N_directions, 2);
	thetasPhis = thetasPhis_all_directions;
	weightsThetasPhis.resize(N_directions);
	weightsThetasPhis = weightsThetasPhis_all_directions;

	thetasPhis_all_directions.resize(0, 0);
	weightsThetasPhis_all_directions.resize(0);

	this->alphaTranslations.resize(boost::extents[Nx][Ny][Nz]);
	this->alphaTranslationsIndexesNonZeros.resize(boost::extents[Nx][Ny][Nz]);
	ArrayXcd alpha(N_directions);
	ArrayXi isAlphaNonZero(N_directions);
	cout << "\n    Process " << ", Level " << this->level << " alpha translations computation" << endl;
	cout << "    alpha.shape() = " << Nx << ", " << Ny << ", " << Nz << ", " << N_directions << endl;
	flush(cout);

	double r_mn[3];
	for (int x = 0; x < Nx; ++x)
	{
		cout << "\r    " << (x + 1) * 100 / Nx << " % computed";
		flush(cout);
		for (int y = 0; y < Ny; ++y)
		{
			for (int z = 0; z < Nz; ++z)
			{
				r_mn[0] = static_cast<double>(x - this->offsetAlphaIndexX);
				r_mn[1] = static_cast<double>(y - this->offsetAlphaIndexY);
				r_mn[2] = static_cast<double>(z - this->offsetAlphaIndexZ);
				if ((abs(r_mn[0]) > 1.0) || (abs(r_mn[1]) > 1.0) || (abs(r_mn[2]) > 1.0))
				{ // if cartesian distance is sufficient
					r_mn[0] *= this->cubeSideLength;
					r_mn[1] *= this->cubeSideLength;
					r_mn[2] *= this->cubeSideLength;
					IT_theta_IT_phi_alpha_C2(alpha, r_mn, getK(), translationOrder, translationOrder_prime, thetasPhis);
					alpha *= weightsThetasPhis;
					//					cout << alpha(3) << endl;
					// we then seek the max of alpha_all_directions
					double max_abs_alpha_local = alpha.abs().maxCoeff();
					double max_abs_alpha = max_abs_alpha_local;
					//					cout << alpha.abs()(3) << endl;
					//					cout << max_abs_alpha << alpha.abs()(0) << endl;
					int countNonZeroLocal = 0, countNonZero;
					for (unsigned int kkk = 0; kkk < alpha.size(); ++kkk)
					{
						if (abs(alpha(kkk)) >= alphaTranslation_thresholdRelValueMax * max_abs_alpha)
						{
							countNonZeroLocal++;
							isAlphaNonZero(kkk) = 1;
						}
						else
						{
							alpha(kkk) = 0.0;
							isAlphaNonZero(kkk) = 0;
						}
					}
					/* To be erased */
					countNonZero = countNonZeroLocal;
					if (countNonZero * 1.0 / (NThetas * NPhis) < alphaTranslation_RelativeCountAboveThreshold)
					{
						this->alphaTranslations[x][y][z].resize(countNonZeroLocal);
						this->alphaTranslationsIndexesNonZeros[x][y][z].resize(countNonZeroLocal);
						int index = 0;
						for (unsigned int kkk = 0; kkk < alpha.size(); ++kkk)
						{
							if (isAlphaNonZero(kkk) == 1)
							{
								this->alphaTranslations[x][y][z](index) = alpha(kkk);
								this->alphaTranslationsIndexesNonZeros[x][y][z](index) = kkk;
								index++;
							}
						}
						if (index != countNonZeroLocal)
						{
							cout << "error in computing alphaTranslations. index = "
								<< index << ", countNonZero = " << countNonZero
								<< ". Exiting..." << endl;
							exit(1);
						}
					}
					else
					{
						this->alphaTranslations[x][y][z].resize(alpha.size());
						this->alphaTranslations[x][y][z] = alpha;
						this->alphaTranslationsIndexesNonZeros[x][y][z].resize(0);
					}
				}
			}
		}
	}
	cout << endl;
	// now we compute the alphaTranslationsIndexes. This is necessary due to the
	// fact that we use the symmetries in alpha computation, which results
	// in less computation time and less memory consumption for exactly the
	// same precision.
	alphaTranslationsIndexes.resize(boost::extents[2][2][2]);
	// if r_p < 0, then corresponding index in alphaTranslationsIndexes is 0.
	// for example, only r_y < 0: coordinates in alphaTranslationsIndexes are: (1, 0, 1)
	ArrayXi oldAlphaIndex(NThetas * NPhis), newAlphaIndexZ(NThetas * NPhis), newAlphaIndexY(NThetas * NPhis), newAlphaIndexX(NThetas * NPhis);
	for (int i = 0; i < NThetas; ++i)
	{
		for (int j = 0; j < NPhis; ++j)
		{
			const int index = i + j * NThetas;
			oldAlphaIndex(index) = index;
		}
	}
	for (int m = 0; m < 2; ++m)
	{
		for (int n = 0; n < 2; ++n)
		{
			for (int p = 0; p < 2; ++p)
			{
				alphaTranslationIndexConstructionX(newAlphaIndexX, oldAlphaIndex, m, NThetas, NPhis);
				alphaTranslationIndexConstructionY(newAlphaIndexY, newAlphaIndexX, n, NThetas, NPhis);
				alphaTranslationIndexConstructionZ(newAlphaIndexZ, newAlphaIndexY, p, NThetas, NPhis);
				alphaTranslationsIndexes[m][n][p] = newAlphaIndexZ;
			}
		}
	}
}

double Level::getAlphaTranslationsSizeMB(void) const
{
	const int Nx = this->alphaTranslations.shape()[0];
	const int Ny = this->alphaTranslations.shape()[1];
	const int Nz = this->alphaTranslations.shape()[2];
	//	cout << Nx << " " << Ny << " " << Nz << endl;
	int N_alpha_elements = 0;
	for (int x = 0; x < Nx; ++x)
	{
		for (int y = 0; y < Ny; ++y)
		{
			for (int z = 0; z < Nz; ++z)
			{
				N_alpha_elements += this->alphaTranslations[x][y][z].size();
				//				cout << this->alphaTranslations[x][y][z].size() << endl;
			}
		}
	}
	const double alphaTranslationsSizeMB = N_alpha_elements * 2.0 * 4.0 / (1024.0 * 1024.0);
	return alphaTranslationsSizeMB;
}

void Level::alphaTranslationIndexConstructionZ(ArrayXi& newAlphaIndex,
                                               const ArrayXi& oldAlphaIndex,
                                               const int alphaCartesianCoordZ,
                                               const int N_theta,
                                               const int N_phi)
/**
 * Symmetry with respect to the \f$ xOy \f$ plane of symmetry. We the have that
 * \f[
 * \alpha \left( \theta > \pi/2, \phi \right) = \alpha \left( \pi - \theta, \phi \right)
 * \f]
 * which, in indexes terms, can be translated as
 *
 */
{
	if (alphaCartesianCoordZ == 0)
	{
		for (int i = 0; i < N_theta; ++i)
		{
			for (int j = 0; j < N_phi; ++j)
			{
				const int index = i + j * N_theta;
				newAlphaIndex(index) = oldAlphaIndex(N_theta - i - 1 + j * N_theta);
			}
		}
	}
	else newAlphaIndex = oldAlphaIndex;
}

void Level::alphaTranslationIndexConstructionY(ArrayXi& newAlphaIndex,
                                               const ArrayXi& oldAlphaIndex,
                                               const int alphaCartesianCoordY,
                                               const int N_theta,
                                               const int N_phi)
/**
 * Symmetry with respect to the \f$ xOz \f$ plane of symmetry. We the have that
 * \f[
 * \alpha \left( \theta, \phi > \pi \right) = \alpha \left( \theta, 2 \pi - \phi \right)
 * \f]
 * which, in indexes terms, can be translated as
 *
 */
{
	if (alphaCartesianCoordY == 0)
	{
		for (int i = 0; i < N_theta; ++i)
		{
			for (int j = 0; j < N_phi; ++j)
			{
				const int index = i + j * N_theta;
				newAlphaIndex(index) = oldAlphaIndex(i + (N_phi - 1 - j) * N_theta);
			}
		}
	}
	else newAlphaIndex = oldAlphaIndex;
}

void Level::alphaTranslationIndexConstructionX(ArrayXi& newAlphaIndex,
                                               const ArrayXi& oldAlphaIndex,
                                               const int alphaCartesianCoordX,
                                               const int N_theta,
                                               const int N_phi)
/**
 * Symmetry with respect to the \f$ yOz \f$ plane of symmetry. We the have that
 * \f[
 * \alpha \left( \theta, \pi/2 < \phi < 3/2 \pi \right)
 * = \alpha \left( \theta, \pi - \phi \right)
 * \f]
 * which, in indexes terms, can be translated as
 *
 */
{
	if (alphaCartesianCoordX == 0)
	{
		for (int i = 0; i < N_theta; ++i)
		{
			for (int j = 0; j < N_phi; ++j)
			{
				const int index = i + j * N_theta;
				if (j < N_phi / 2)
					newAlphaIndex(index) = oldAlphaIndex(i + (N_phi / 2 - 1 - j) * N_theta);
				else
					newAlphaIndex(index)
						= oldAlphaIndex(i + (N_phi - 1 - (j - N_phi / 2)) * N_theta);
			}
		}
	}
	else newAlphaIndex = oldAlphaIndex;
}

void Level::shiftingArraysComputation(void)
/**
 * This is the function that computes the shifting coefficients
 from each corner of the cube
 * towards its center.
 * Therefore, there are as many shifting coefficients as there are elements
 * at the given level times 8, since the cube has 8 corners.
 *
 * If one wants to go from the center towards a corner,
 * as it happens in down-shifting,
 * the coefficients must be conjugated.
 *
 * The "shiftingArrays" is in fact a 2D array,
 * where the first column corresponds to the (0,0,0) corner
 * and the last column is the (1,1,1) corner. The numbering is as follows:
 * (0,0,0) -> column 0
 * (0,0,1) -> column 1
 * (0,1,0) -> column 2
 * (0,1,1) -> column 3
 * (1,0,0) -> column 4
 * (1,0,1) -> column 5
 * (1,1,0) -> column 6
 * (1,1,1) -> column 7
 */
{
	const int N_theta = thetas.size(), N_phi = phis.size();
	shiftingArrays.resize(8, N_theta * N_phi);
	double sinTheta, cosTheta, shiftingDistance = this->cubeSideLength / 4.0;
	double k_hat[3], shiftingVector[3];
	Array<ArrayXd, Dynamic, 1> shiftingVectors(8);
	for (int p = 0; p < 8; ++p) shiftingVectors(p).resize(3);
	shiftingVectors(0) << -1.0 , -1.0 , -1.0;
	shiftingVectors(1) << -1.0 , -1.0 , 1.0;
	shiftingVectors(2) << -1.0 , 1.0 , -1.0;
	shiftingVectors(3) << -1.0 , 1.0 , 1.0;
	shiftingVectors(4) << 1.0 , -1.0 , -1.0;
	shiftingVectors(5) << 1.0 , -1.0 , 1.0;
	shiftingVectors(6) << 1.0 , 1.0 , -1.0;
	shiftingVectors(7) << 1.0 , 1.0 , 1.0;
	for (int p = 0; p < 8; ++p)
	{
		for (int i = 0; i < 3; i++) shiftingVector[i] = shiftingVectors(p)(i) * shiftingDistance;
		for (int m = 0; m < N_theta; ++m)
		{
			sinTheta = sin(thetas(m));
			cosTheta = cos(thetas(m));
			for (int n = 0; n < N_phi; ++n)
			{
				k_hat[0] = sinTheta * cos(phis(n));
				k_hat[1] = sinTheta * sin(phis(n));
				k_hat[2] = cosTheta;
				const double temp(k_hat[0] * shiftingVector[0] + k_hat[1] * shiftingVector[1] + k_hat[2] * shiftingVector[2]);
				shiftingArrays(p, m + n * N_theta) = static_cast<complex<double>>(exp(-I * k * temp));
			}
		}
	}
}

Array3i Level::getAlphaTranslationsExtents(void) const
{
	Array3i dimensions;
	dimensions << alphaTranslations.shape()[0] ,
		alphaTranslations.shape()[1] ,
		alphaTranslations.shape()[2];
	return dimensions;
}

void Level::sortCubesByParents(void)
{
	sort(cubes.begin(), cubes.end());
	for (int j = 0; j < getLevelSize(); ++j) cubes[j].setIndex(j);
	// we write the new values of the indexes
	sort(cubes_obs.begin(), cubes_obs.end());
	for (unsigned int j = 0; j < cubes_obs.size(); ++j) cubes_obs[j].setIndex(j);
	// we write the new values of the indexes
}

void Level::updateFatherIndexes(const Level& fatherLevel)
{
	vector<int> sonsIndexes;
	for (int i = 0; i < fatherLevel.getLevelSize(); ++i)
	{
		sonsIndexes = fatherLevel.getCube(i).sonsIndexes;
		for (unsigned int j = 0; j < sonsIndexes.size(); ++j)
		{
			if (cubes[sonsIndexes[j]].getFatherNumber()
				!= fatherLevel.getCube(i).getNumber())
			{
				cout << "Level::updateFatherIndexes: "
					<< "father number in sons and father do not match!" << endl;
				exit(1);
			}
			cubes[sonsIndexes[j]].setFatherIndex(i);
		}
	}

	for (unsigned int i = 0; i < fatherLevel.cubes_obs.size(); ++i)
	{
		sonsIndexes = fatherLevel.cubes_obs[i].sonsIndexes;
		for (unsigned int j = 0; j < sonsIndexes.size(); ++j)
		{
			if (cubes_obs[sonsIndexes[j]].getFatherNumber()
				!= fatherLevel.cubes_obs[i].getNumber())
			{
				cout << "Level::updateFatherIndexes: "
					<< "father number in sons and father do not match!" << endl;
				exit(1);
			}
			cubes_obs[sonsIndexes[j]].setFatherIndex(i);
		}
	}
}

void Level::setNumbersToIndexes(void)
{
	const int N_cubes = getLevelSize(), N_cubes_obs = cubes_obs.size();
	int i, n;
	numbersToIndexes.reserve(N_cubes);
	for (i = 0; i < N_cubes; ++i)
	{
		n = getCube(i).getNumber();
		numbersToIndexes.push_back(Dictionary<int, int>(n, i));
	}
	for (i = 0; i < N_cubes_obs; ++i)
	{
		n = cubes_obs[i].getNumber();
		numbersToIndexes_obs.push_back(Dictionary<int, int>(n, i));
	}
}

void Level::sortNumbersToIndexes(void)
{
	sort(numbersToIndexes.begin(), numbersToIndexes.end());
	sort(numbersToIndexes_obs.begin(), numbersToIndexes_obs.end());
}

//void Level::printNumbersToIndexes(void)
//{
//	const int N_cubes = getLevelSize();
//	const int l = getLevel();
//	for (int i = 0; i < N_cubes; ++i)
//	{
//		cout << "Level " << l << " : numbersToIndexes["
//			<< i << "] = " << numbersToIndexes[i].getKey()
//			<< ", " << numbersToIndexes[i].getVal() << endl;
//	}
//}

void Level::searchCubesNeighborsIndexes(void)
{
	setNumbersToIndexes();
	sortNumbersToIndexes();
	const int N_cubes = getLevelSize();
	for (int i = 0; i < N_cubes; ++i)
	{
		const double* absCartCoord(cubes[i].absoluteCartesianCoord);
		vector<int> neighborsIndexes;
		// we find the neighbors
		for (int x = -1; x < 2; ++x)
		{
			for (int y = -1; y < 2; ++y)
			{
				for (int z = -1; z < 2; ++z)
				{
					int index = -1;
					const double CandidateAbsCartCoord[3] = {absCartCoord[0] + x, absCartCoord[1] + y, absCartCoord[2] + z};
					/// no component of (absoluteCartesianCoord(i) + p) -- where i=0,1,2 and p = x,y,z -- can be:
					/// (1) negative or (2) greater than MaxNumberCubes1D.
					int condition = 1;
					for (int j = 0; j < 3; ++j) condition *= ((CandidateAbsCartCoord[j] >= 0) && (CandidateAbsCartCoord[j] < getMaxNumberCubes1D()));
					/* we also do not want to consider the cube itself */
					condition *= !((x == 0) && (y == 0) && (z == 0));
					if (condition > 0)
					{
						int candidate_number
							= static_cast<int>(CandidateAbsCartCoord[0] * pow(getMaxNumberCubes1D(), 2)
								+ CandidateAbsCartCoord[1] * getMaxNumberCubes1D()
								+ CandidateAbsCartCoord[2]);
						index = getIndexOfNumber(candidate_number);
					}
					if (index > -1) neighborsIndexes.push_back(index);
				}
			}
		}
		// we now trim the excess capacity of neighborsIndexes
		cubes[i].neighborsIndexes.resize(neighborsIndexes.size());
		for (unsigned int j = 0; j < neighborsIndexes.size(); ++j)
		{
			cubes[i].neighborsIndexes[j] = neighborsIndexes[j];
		}
	}

	const int N_cubes_obs = cubes_obs.size();
	for (int i = 0; i < N_cubes_obs; ++i)
	{
		const double* absCartCoord(cubes_obs[i].absoluteCartesianCoord);
		vector<int> neighborsIndexes;
		// we find the neighbors
		for (int x = -1; x < 2; ++x)
		{
			for (int y = -1; y < 2; ++y)
			{
				for (int z = -1; z < 2; ++z)
				{
					int index = -1;
					const double CandidateAbsCartCoord[3] = {absCartCoord[0] + x, absCartCoord[1] + y, absCartCoord[2] + z};
					/// no component of (absoluteCartesianCoord(i) + p) -- where i=0,1,2 and p = x,y,z -- can be:
					/// (1) negative or (2) greater than MaxNumberCubes1D.
					int condition = 1;
					for (int j = 0; j < 3; ++j) condition *= ((CandidateAbsCartCoord[j] >= 0) && (CandidateAbsCartCoord[j] < getMaxNumberCubes1D()));
					/* we also do not want to consider the cube itself */
					condition *= !((x == 0) && (y == 0) && (z == 0));
					if (condition > 0)
					{
						int candidate_number = static_cast<int>(CandidateAbsCartCoord[0] * pow(getMaxNumberCubes1D(), 2) + CandidateAbsCartCoord[1] * getMaxNumberCubes1D() + CandidateAbsCartCoord[2]);
						index = getIndexOfNumber(candidate_number);
					}
					if (index > -1) neighborsIndexes.push_back(index);
				}
			}
		}
		// we now trim the excess capacity of neighborsIndexes
		cubes_obs[i].neighborsIndexes.resize(neighborsIndexes.size());
		for (unsigned int j = 0; j < neighborsIndexes.size(); ++j)
		{
			cubes_obs[i].neighborsIndexes[j] = neighborsIndexes[j];
		}
	}
}

int Level::getIndexOfNumber(const int number) const // "numbersToIndexes" must be ordered
// this is done in the calling function...
{
	int ind_inf, ind_sup, ind_mid, index;
	const int N = getNumbersToIndexesSize();
	if ((number < getNumberToIndex(0)) || (number > getNumberToIndex(N - 1))) index = -1;
	else
	{
		ind_inf = 0;
		ind_sup = N - 1;
		while (ind_sup - ind_inf > 1)
		{
			ind_mid = (ind_sup + ind_inf) / 2;
			if (number > getNumberToIndex(ind_mid)) ind_inf = ind_mid;
			else ind_sup = ind_mid;
		}
		if (number == getNumberToIndex(ind_inf)) index = getIndexToIndex(ind_inf);
		else if (number == getNumberToIndex(ind_sup)) index = getIndexToIndex(ind_sup);
		else index = -1;
	}
	return index;
}

void Level::computeOldIndexesOfCubes(ArrayXi& oldIndexesOfCubes,
                                     ArrayXi& oldIndexesOfCubes_obs)
{
	// this function computes the indexes of the cubes in the original mesh
	int NCubes = cubes.size();
	oldIndexesOfCubes.resize(NCubes);
	for (int i = 0; i < NCubes; ++i)
	{
		oldIndexesOfCubes(i) = cubes[i].getOldIndex();
	}
	NCubes = cubes_obs.size();
	oldIndexesOfCubes_obs.resize(NCubes);
	for (int i = 0; i < NCubes; ++i)
	{
		oldIndexesOfCubes_obs(i) = cubes_obs[i].getOldIndex();
	}
}

void Level::computeGaussLocatedArguments
(const ArrayXi& local_cubes_NRWG,
 const ArrayXi& local_RWG_numbers,
 const ArrayXi& local_RWG_Numbers_surf,
 const ArrayXXi& local_RWGNumbers_signedTriangles,
 const ArrayXXd& local_RWGNumbers_trianglesCoord,
 const int N_Gauss)
{
	if (getLeaf())
	{
		const int NCubes = cubes.size();
		int startIndex_in_localArrays = 0;
		for (int i = 0; i < NCubes; ++i)
		{
			const int NRWG = local_cubes_NRWG(i);
			cubes[i].computeGaussLocatedArguments
			(local_RWG_numbers, local_RWG_Numbers_surf,
			 local_RWGNumbers_signedTriangles, local_RWGNumbers_trianglesCoord,
			 startIndex_in_localArrays, NRWG, N_Gauss);
			startIndex_in_localArrays += NRWG;
		}
	}
}

void Level::computeObsLocatedArguments(const ArrayXi& cubes_Nobs,
                                       const ArrayXi& obs_numbers,
                                       const ArrayXXd& obsNumber_obsCoord)
{
	if (getLeaf())
	{
		const int NCubes = cubes_obs.size();
		int startIndex = 0;
		for (int i = 0; i < NCubes; ++i)
		{
			const int Nobs = cubes_Nobs(i);
			cubes_obs[i].obs_numbers.resize(Nobs);
			cubes_obs[i].r_obs.resize(Nobs, 3);
			for (int j = 0; j < Nobs; ++j)
			{
				int obs_number = obs_numbers(startIndex + j);
				cubes_obs[i].obs_numbers(j) = obs_number;
				cubes_obs[i].r_obs.row(j) = obsNumber_obsCoord.row(obs_number);
			}
			startIndex += Nobs;
		}
	}
}

void Level::RWGs_renumbering(void)
{
	if (getLeaf())
	{
		const int N_cubes = cubes.size();
		int startIndex = 0;
		for (int i = 0; i < N_cubes; ++i)
		{
			const int NRWG = cubes[i].RWG_numbers.size();
			for (int j = 0; j < NRWG; j++) cubes[i].RWG_numbers[j] = startIndex + j;
			startIndex += NRWG;
		}
	}
}

void Level::computeSup(ArrayXXcd& Sup,
                       const complex<double>& k,
                       const ArrayXcd& I_PQ,
                       const Cube& cube,
                       const ArrayXd& thetas,
                       const ArrayXd& phis,
                       const ArrayXi& target_surface)
{
	const int NThetas = thetas.size(), NPhis = phis.size();
	const int NGauss = cube.triangle_GaussCoord.cols() / 3;

	// A. Francavilla (29-05-2013)
	double sum_weigths;
	const double *xi, *eta, *weigths;
	IT_points(xi, eta, weigths, sum_weigths, NGauss);

	ArrayXXd kHats(NThetas * NPhis, 3), thetaHats(NThetas * NPhis, 3);
	ArrayXXd phiHats(NThetas * NPhis, 3);
	ArrayXXcd FC3Components(NThetas * NPhis, 3);
	vector<double> sin_thetas, cos_thetas;
	sin_thetas.resize(NThetas);
	cos_thetas.resize(NThetas);

	for (int p = 0; p < NThetas; ++p)
	{
		sin_thetas[p] = sin(thetas(p));
		cos_thetas[p] = cos(thetas(p));
	}

	// initialisation of arrays
	for (int q = 0; q < NPhis; ++q)
	{
		const double cos_phi = cos(phis(q)), sin_phi = sin(phis(q));
		for (int p = 0; p < NThetas; ++p)
		{
			int index = p + q * NThetas;
			const double sin_theta = sin_thetas[p], cos_theta = cos_thetas[p];
			kHats(index, 0) = sin_theta * cos_phi;
			kHats(index, 1) = sin_theta * sin_phi;
			kHats(index, 2) = cos_theta;
			thetaHats(index, 0) = cos_theta * cos_phi;
			thetaHats(index, 1) = cos_theta * sin_phi;
			thetaHats(index, 2) = -sin_theta;
			phiHats(index, 0) = -sin_phi;
			phiHats(index, 1) = cos_phi;
			phiHats(index, 2) = 0.0;
			for (int i = 0; i < 3; i++) FC3Components(index, i) = 0.0;
		}
	}

	// computation of FC3Components array
	const complex<double> I_k(static_cast<complex<double>>(I * k));
	const double* rCenter = cube.rCenter;
	const int T = cube.Triangle_numberOfRWGs.size();
	int startIndex = 0, startIndex_r_opp = 0;

	for (int i = 0; i < T; i++)
	{
		if ((target_surface != cube.TriangleToSurf[i]).all()) continue;
		const int n_rwg = cube.Triangle_numberOfRWGs[i];

		for (int j = 0; j < NGauss; j++)
		{
			const double r[3] = {cube.triangle_GaussCoord(i, j * 3),
				cube.triangle_GaussCoord(i, j * 3 + 1),
				cube.triangle_GaussCoord(i, j * 3 + 2)};
			complex<double> fj[3] = {0.0, 0.0, 0.0};

			//			if(cube.getIndex()==0)
			//			{
			//				for (int i = 0; i < 3; i++) cout << r[i] << " ";
			//				cout << endl;
			//			}

			// loop on the RWGs for triangle i
			for (int rwg = 0; rwg < n_rwg; rwg++)
			{
				const int RWG_index = cube.TriangleToRWGindex[startIndex + rwg];
				// if (all(target_surface != cube.RWG_numbers_surf[RWG_index])) continue;
				const double weight
					= cube.TriangleToRWGweight[startIndex + rwg] * weigths[j];
				const complex<double> i_pq = I_PQ(cube.RWG_numbers[RWG_index]) * weight;
				const int index = startIndex_r_opp + rwg * 3;
				fj[0] += i_pq * (r[0] - cube.TriangleToRWG_ropp[index]);
				fj[1] += i_pq * (r[1] - cube.TriangleToRWG_ropp[index + 1]);
				fj[2] += i_pq * (r[2] - cube.TriangleToRWG_ropp[index + 2]);
			} // end loop RWGs

			const double nHat[3] = {cube.triangle_nHat(i, 0),
				cube.triangle_nHat(i, 1),
				cube.triangle_nHat(i, 2)};

			const double expArg[3] = {r[0] - rCenter[0],
				r[1] - rCenter[1], r[2] - rCenter[2]};

			for (int q = 0; q < NPhis / 2; q++)
			{
				const int index_1 = q * NThetas;
				const int opp_index_1 = NThetas - 1 + (q + NPhis / 2) * NThetas;
				for (int p = 0; p < NThetas; p++)
				{
					const int index = p + index_1, opp_index = opp_index_1 - p;
					const complex<double> a(I_k * (expArg[0] * kHats(index, 0)
						+ expArg[1] * kHats(index, 1)
						+ expArg[2] * kHats(index, 2)));
					double c, s, e;
					e = (a.real() == 0.0) ? 1.0 : exp(a.real());
					//					sincosf(a.imag(), &s, &c);
					s = sin(a.imag());
					c = cos(a.imag());
					const complex<double> EXP(e * c, e * s);
					const complex<double> conjEXP((1.0 / e) * c, -(1.0 / e) * s);
					FC3Components(index, 0) += fj[0] * EXP;
					FC3Components(index, 1) += fj[1] * EXP;
					FC3Components(index, 2) += fj[2] * EXP;
					FC3Components(opp_index, 0) += fj[0] * conjEXP;
					FC3Components(opp_index, 1) += fj[1] * conjEXP;
					FC3Components(opp_index, 2) += fj[2] * conjEXP;
				}
			} // end q loop
		} // end j Gauss loop

		startIndex += n_rwg;
		startIndex_r_opp += n_rwg * 3;
	}

	for (int i = 0; i < Sup.cols(); ++i)
	{
		Sup(0, i) = thetaHats(i, 0) * FC3Components(i, 0)
			+ thetaHats(i, 1) * FC3Components(i, 1)
			+ thetaHats(i, 2) * FC3Components(i, 2);
		Sup(1, i) = phiHats(i, 0) * FC3Components(i, 0)
			+ phiHats(i, 1) * FC3Components(i, 1)
			+ phiHats(i, 2) * FC3Components(i, 2);
	}
}

void Level::computeSup_obs(ArrayXXcd& Sup,
                           const complex<double>& k,
                           const ArrayXXcd& E_obs,
                           const Cube& cube,
                           const ArrayXd& thetas,
                           const ArrayXd& phis,
                           const bool M_current)
{
	// const int NThetas = thetas.size(),
	// NPhis = phis.size(), NGauss = cube.triangle_GaussCoord.extent(1)/3;
	const int NThetas = thetas.size(), NPhis = phis.size();

	ArrayXXd kHats(NThetas * NPhis, 3), thetaHats(NThetas * NPhis, 3), phiHats(NThetas * NPhis, 3);
	ArrayXXcd FC3Components(NThetas * NPhis, 3);
	vector<double> sin_thetas, cos_thetas;

	sin_thetas.resize(NThetas);
	cos_thetas.resize(NThetas);

	for (int p = 0; p < NThetas; ++p)
	{
		sin_thetas[p] = sin(thetas(p));
		cos_thetas[p] = cos(thetas(p));
	}

	// initialisation of arrays
	for (int q = 0; q < NPhis; ++q)
	{
		const double cos_phi = cos(phis(q)), sin_phi = sin(phis(q));
		for (int p = 0; p < NThetas; ++p)
		{
			int index = p + q * NThetas;
			const double sin_theta = sin_thetas[p], cos_theta = cos_thetas[p];
			kHats(index, 0) = sin_theta * cos_phi;
			kHats(index, 1) = sin_theta * sin_phi;
			kHats(index, 2) = cos_theta;
			thetaHats(index, 0) = cos_theta * cos_phi;
			thetaHats(index, 1) = cos_theta * sin_phi;
			thetaHats(index, 2) = -sin_theta;
			phiHats(index, 0) = -sin_phi;
			phiHats(index, 1) = cos_phi;
			phiHats(index, 2) = 0.0;
			for (int i = 0; i < 3; i++) FC3Components(index, i) = 0.0;
		}
	}

	// computation of integration
	// defining local arrays used for faster computations
	const complex<double> I_k(static_cast<complex<double>>(I * k));
	const double* rCenter = cube.rCenter;
	const int Nobs = cube.obs_numbers.size();

	for (int i = 0; i < Nobs; i++)
	{
		const double r[3] = {cube.r_obs(i, 0), cube.r_obs(i, 1), cube.r_obs(i, 2)};
		// computation of the shifting terms
		const double expArg[3] = {r[0] - rCenter[0], r[1] - rCenter[1], r[2] - rCenter[2]};
		const int obsNumber = cube.obs_numbers(i);

		// for phi>pi, kHat = -kHat(pi-theta, phi-pi)
		for (int q = 0; q < NPhis / 2; q++)
		{
			const int index_1 = q * NThetas,
				opp_index_1 = NThetas - 1 + (q + NPhis / 2) * NThetas;
			for (int p = 0; p < NThetas; p++)
			{
				const int index = p + index_1, opp_index = opp_index_1 - p;
				const complex<double>
					a(I_k * (expArg[0] * kHats(index, 0)
						+ expArg[1] * kHats(index, 1)
						+ expArg[2] * kHats(index, 2)));

				double c, s, e;
				e = (a.real() == 0.0) ? 1.0 : exp(a.real());
				s = sin(a.imag());
				c = cos(a.imag());

				const complex<double> EXP(e * c, e * s);
				FC3Components(index, 0) += E_obs(obsNumber, 0) * EXP;
				FC3Components(index, 1) += E_obs(obsNumber, 1) * EXP;
				FC3Components(index, 2) += E_obs(obsNumber, 2) * EXP;
				const complex<double> conjEXP((1.0 / e) * c, -(1.0 / e) * s);
				FC3Components(opp_index, 0) += E_obs(obsNumber, 0) * conjEXP;
				FC3Components(opp_index, 1) += E_obs(obsNumber, 1) * conjEXP;
				FC3Components(opp_index, 2) += E_obs(obsNumber, 2) * conjEXP;
			}
		}
	}

	for (int i = 0; i < Sup.cols(); ++i)
	{
		Sup(0, i) = thetaHats(i, 0) * FC3Components(i, 0)
			+ thetaHats(i, 1) * FC3Components(i, 1)
			+ thetaHats(i, 2) * FC3Components(i, 2);
		Sup(1, i) = phiHats(i, 0) * FC3Components(i, 0)
			+ phiHats(i, 1) * FC3Components(i, 1)
			+ phiHats(i, 2) * FC3Components(i, 2);
	}
}

void Level::sphericalIntegration(ArrayXcd& ZI,
                                 const ArrayXXcd& Sdown,
                                 const Cube& cube,
                                 const ArrayXd& thetas,
                                 const ArrayXd& phis,
                                 const double w,
                                 const complex<double>& mu_r,
                                 const complex<double>& k,
                                 const ArrayXcd& CFIE,
                                 const bool M_current,
                                 const ArrayXi& target_surface)
{
	const complex<double> ZZERO(0.0, 0.0);
	const complex<double> EJ_factor(-I * mu_0 * w * mu_r * CFIE(0));
	const complex<double> EM_factor(-I * k * CFIE(0));

	const int NThetas = thetas.size(), NPhis = phis.size();
	const int NGauss = cube.triangle_GaussCoord.cols() / 3;

	// A. Francavilla (29-05-2013)
	double sum_weigths;
	const double *xi, *eta, *weigths;
	IT_points(xi, eta, weigths, sum_weigths, NGauss);

	ArrayXXd kHats(NThetas * NPhis, 3);
	ArrayXXcd GC3Components(NThetas * NPhis, 3), GC3Exp(NThetas * NPhis, 3);
	vector<double> sin_thetas, cos_thetas;
	sin_thetas.resize(NThetas);
	cos_thetas.resize(NThetas);
	for (int p = 0; p < NThetas; ++p)
	{
		sin_thetas[p] = sin(thetas(p));
		cos_thetas[p] = cos(thetas(p));
	}
	// initialisation of arrays
	//	cout << Sdown.row(0) << endl;
	for (int q = 0; q < NPhis; ++q)
	{
		const double cos_phi = cos(phis(q)), sin_phi = sin(phis(q));
		const double phiHat[3] = {-sin_phi, cos_phi, 0.0};
		for (int p = 0; p < NThetas; ++p)
		{
			int index = p + q * NThetas;
			const double sin_theta = sin_thetas[p], cos_theta = cos_thetas[p];
			kHats(index, 0) = sin_theta * cos_phi;
			kHats(index, 1) = sin_theta * sin_phi;
			kHats(index, 2) = cos_theta;
			const double thetaHat[3] = {cos_theta * cos_phi,
				cos_theta * sin_phi, -sin_theta};
			GC3Components(index, 0)
				= Sdown(0, index) * thetaHat[0] + Sdown(1, index) * phiHat[0];
			GC3Components(index, 1)
				= Sdown(0, index) * thetaHat[1] + Sdown(1, index) * phiHat[1];
			GC3Components(index, 2)
				= Sdown(0, index) * thetaHat[2] + Sdown(1, index) * phiHat[2];
		}
	}
	// computation of integration
	// defining local arrays used for faster computations
	const complex<double> minus_I_k(static_cast<complex<double>>(-I * k));
	const double* rCenter = cube.rCenter;
	const int T = cube.Triangle_numberOfRWGs.size();
	int startIndex = 0, startIndex_r_opp = 0;
	for (int i = 0; i < T; i++)
	{
		if ((target_surface != cube.TriangleToSurf[i]).all()) continue;
		const int n_rwg = cube.Triangle_numberOfRWGs[i];
		const double nHat[3] = {cube.triangle_nHat(i, 0),
			cube.triangle_nHat(i, 1), cube.triangle_nHat(i, 2)};

		for (int j = 0; j < NGauss; j++)
		{
			const double r[3] = {cube.triangle_GaussCoord(i, j * 3),
				cube.triangle_GaussCoord(i, j * 3 + 1),
				cube.triangle_GaussCoord(i, j * 3 + 2)};

			// computation of the shifting terms
			const double expArg[3]
				= {r[0] - rCenter[0], r[1] - rCenter[1], r[2] - rCenter[2]};
			complex<double> EJ[3] = {0.0, 0.0, 0.0};

			// for phi>pi, kHat = -kHat(pi-theta, phi-pi)
			for (int q = 0; q < NPhis / 2; q++)
			{
				const int index_1 = q * NThetas,
					opp_index_1 = NThetas - 1 + (q + NPhis / 2) * NThetas;

				for (int p = 0; p < NThetas; p++)
				{
					const int index = p + index_1, opp_index = opp_index_1 - p;
					const complex<double>
						a(minus_I_k * (expArg[0] * kHats(index, 0)
							+ expArg[1] * kHats(index, 1)
							+ expArg[2] * kHats(index, 2)));
					double c, s, e;
					e = (a.real() == 0.0) ? 1.0 : exp(a.real());
					s = sin(a.imag());
					c = cos(a.imag());
					const complex<double> EXP(e * c, e * s);
					GC3Exp(index, 0) = GC3Components(index, 0) * EXP;
					GC3Exp(index, 1) = GC3Components(index, 1) * EXP;
					GC3Exp(index, 2) = GC3Components(index, 2) * EXP;
					const complex<double> conjEXP((1.0 / e) * c, -(1.0 / e) * s);
					GC3Exp(opp_index, 0) = GC3Components(opp_index, 0) * conjEXP;
					GC3Exp(opp_index, 1) = GC3Components(opp_index, 1) * conjEXP;
					GC3Exp(opp_index, 2) = GC3Components(opp_index, 2) * conjEXP;
					EJ[0] += (GC3Exp(index, 0) + GC3Exp(opp_index, 0));
					EJ[1] += (GC3Exp(index, 1) + GC3Exp(opp_index, 1));
					EJ[2] += (GC3Exp(index, 2) + GC3Exp(opp_index, 2));
				}
			}
			//			cout << EJ[0] << " " << EJ[1] << " " << EJ[2] << endl;

			// now for the M_current
			complex<double> GC3_x_kHat[3] = {0.0, 0.0, 0.0};
			if (M_current)
			{
				// we first have to compute GC3_x_kHat
				// for phi>pi, kHat = -kHat(pi-theta, phi-pi)
				for (int q = 0; q < NPhis; q++)
				{
					const int index_1 = q * NThetas;
					for (int p = 0; p < NThetas; p++)
					{
						const int index = p + index_1;
						GC3_x_kHat[0] += GC3Exp(index, 1) * kHats(index, 2)
							- GC3Exp(index, 2) * kHats(index, 1);
						GC3_x_kHat[1] += GC3Exp(index, 2) * kHats(index, 0)
							- GC3Exp(index, 0) * kHats(index, 2);
						GC3_x_kHat[2] += GC3Exp(index, 0) * kHats(index, 1)
							- GC3Exp(index, 1) * kHats(index, 0);
					}
				}
			} // end if for M_current

			// loop on the RWGs for triangle i
			for (int rwg = 0; rwg < n_rwg; rwg++)
			{
				// common EFIE and MFIE
				const int RWG_index = cube.TriangleToRWGindex[startIndex + rwg];
				// if (all(target_surface != cube.RWG_numbers_surf[RWG_index])) continue;

				const double weight
					= cube.TriangleToRWGweight[startIndex + rwg] * weigths[j];
				const int RWGNumber = cube.RWG_numbers[RWG_index];
				// EFIE
				const int index = startIndex_r_opp + rwg * 3;
				const double fj[3] = {(r[0] - cube.TriangleToRWG_ropp[index]),
					(r[1] - cube.TriangleToRWG_ropp[index + 1]),
					(r[2] - cube.TriangleToRWG_ropp[index + 2])};

				if (M_current)
					ZI(RWGNumber)
						+= EM_factor * (fj[0] * GC3_x_kHat[0] + fj[1] * GC3_x_kHat[1]
							+ fj[2] * GC3_x_kHat[2]) * weight;
				else
					ZI(RWGNumber) += EJ_factor * (EJ[0] * fj[0] + EJ[1] * fj[1]
						+ EJ[2] * fj[2]) * weight;
			}// end loop on the RWGs
		}// end Gauss integration points loop

		startIndex += n_rwg;
		startIndex_r_opp += n_rwg * 3;
	}// end triangles loop
}

void Level::sphericalIntegration_obs(ArrayXXcd& E_obs,
                                     const ArrayXXcd& Sdown,
                                     const Cube& cube,
                                     const ArrayXd& thetas,
                                     const ArrayXd& phis,
                                     const double w,
                                     const complex<double>& mu_r,
                                     const complex<double>& k,
                                     const ArrayXcd& CFIE,
                                     const bool Mcurrent)
{
	const int NThetas = thetas.size(), NPhis = phis.size();
	const complex<double> EJ_factor(-I * mu_0 * w * mu_r);
	const complex<double> EM_factor(-I * k);

	ArrayXXd kHats(NThetas * NPhis, 3);
	ArrayXXcd GC3Components(NThetas * NPhis, 3), GC3Exp(NThetas * NPhis, 3);
	vector<double> sin_thetas, cos_thetas;
	sin_thetas.resize(NThetas);
	cos_thetas.resize(NThetas);
	for (int p = 0; p < NThetas; ++p)
	{
		sin_thetas[p] = sin(thetas(p));
		cos_thetas[p] = cos(thetas(p));
	}
	// initialisation of arrays
	for (int q = 0; q < NPhis; ++q)
	{
		const double cos_phi = cos(phis(q)), sin_phi = sin(phis(q));
		const double phiHat[3] = {-sin_phi, cos_phi, 0.0};
		for (int p = 0; p < NThetas; ++p)
		{
			int index = p + q * NThetas;
			const double sin_theta = sin_thetas[p], cos_theta = cos_thetas[p];
			kHats(index, 0) = sin_theta * cos_phi;
			kHats(index, 1) = sin_theta * sin_phi;
			kHats(index, 2) = cos_theta;
			const double
				thetaHat[3] = {cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta};
			GC3Components(index, 0)
				= Sdown(0, index) * thetaHat[0] + Sdown(1, index) * phiHat[0];
			GC3Components(index, 1)
				= Sdown(0, index) * thetaHat[1] + Sdown(1, index) * phiHat[1];
			GC3Components(index, 2)
				= Sdown(0, index) * thetaHat[2] + Sdown(1, index) * phiHat[2];
		}
	}
	// computation of integration
	// defining local arrays used for faster computations
	const complex<double> minus_I_k(static_cast<complex<double>>(-I * k));
	const double* rCenter = cube.rCenter;
	const int Nobs = cube.obs_numbers.size();
	for (int i = 0; i < Nobs; i++)
	{
		const double r[3] = {cube.r_obs(i, 0), cube.r_obs(i, 1), cube.r_obs(i, 2)};

		// computation of the shifting terms
		const double expArg[3] = {r[0] - rCenter[0], r[1] - rCenter[1], r[2] - rCenter[2]};
		const int obsNumber = cube.obs_numbers(i);

		// for phi>pi, kHat = -kHat(pi-theta, phi-pi)
		for (int q = 0; q < NPhis / 2; q++)
		{
			const int index_1 = q * NThetas;
			const int opp_index_1 = NThetas - 1 + (q + NPhis / 2) * NThetas;

			for (int p = 0; p < NThetas; p++)
			{
				const int index = p + index_1, opp_index = opp_index_1 - p;
				const complex<double>
					a(minus_I_k * (expArg[0] * kHats(index, 0)
						+ expArg[1] * kHats(index, 1)
						+ expArg[2] * kHats(index, 2)));
				double c, s, e;
				e = (a.real() == 0.0) ? 1.0 : exp(a.real());
				s = sin(a.imag());
				c = cos(a.imag());
				const complex<double> EXP(e * c, e * s);
				GC3Exp(index, 0) = GC3Components(index, 0) * EXP;
				GC3Exp(index, 1) = GC3Components(index, 1) * EXP;
				GC3Exp(index, 2) = GC3Components(index, 2) * EXP;
				const complex<double> conjEXP((1.0 / e) * c, -(1.0 / e) * s);
				GC3Exp(opp_index, 0) = GC3Components(opp_index, 0) * conjEXP;
				GC3Exp(opp_index, 1) = GC3Components(opp_index, 1) * conjEXP;
				GC3Exp(opp_index, 2) = GC3Components(opp_index, 2) * conjEXP;

				if (!Mcurrent)
				{
					E_obs(obsNumber, 0)
						+= EJ_factor * (GC3Exp(index, 0) + GC3Exp(opp_index, 0));
					E_obs(obsNumber, 1)
						+= EJ_factor * (GC3Exp(index, 1) + GC3Exp(opp_index, 1));
					E_obs(obsNumber, 2)
						+= EJ_factor * (GC3Exp(index, 2) + GC3Exp(opp_index, 2));
				}
			}
		}

		// now for the M current
		if (Mcurrent)
		{
			// we first have to compute GC3_x_kHat
			// for phi>pi, kHat = -kHat(pi-theta, phi-pi)
			complex<double> GC3_x_kHat[3] = {0.0, 0.0, 0.0};
			for (int q = 0; q < NPhis; q++)
			{
				const int index_1 = q * NThetas;
				for (int p = 0; p < NThetas; p++)
				{
					const int index = p + index_1;
					GC3_x_kHat[0] += GC3Exp(index, 1) * kHats(index, 2)
						- GC3Exp(index, 2) * kHats(index, 1);
					GC3_x_kHat[1] += GC3Exp(index, 2) * kHats(index, 0)
						- GC3Exp(index, 0) * kHats(index, 2);
					GC3_x_kHat[2] += GC3Exp(index, 0) * kHats(index, 1)
						- GC3Exp(index, 1) * kHats(index, 0);
				}
			}

			E_obs(obsNumber, 0) += EM_factor * GC3_x_kHat[0];
			E_obs(obsNumber, 1) += EM_factor * GC3_x_kHat[1];
			E_obs(obsNumber, 2) += EM_factor * GC3_x_kHat[2];
		} // end if for M current
	}
}
