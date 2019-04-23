/**********************************************************************
 *
 * level.h
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

#ifndef LEVEL_H
#define LEVEL_H

#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>  // Include STL algorithms for sorting lists and vectors

using namespace std;

#include "dictionary.hpp"
#include "cube.hpp"
#include "interpolation.hpp"

class Level
{
	bool leaf;
	bool ceiling;
	int level;
	int maxNumberCubes1D;
	int NCubesX;
	int NCubesY;
	int NCubesZ;
	int offsetAlphaIndexX, offsetAlphaIndexY, offsetAlphaIndexZ;
	int N;
	double cubeSideLength;
	std::complex<double> k;
	int numberTimesCopied;
public:
	vector<Cube> cubes, cubes_obs;
	vector<Dictionary<int, int>> numbersToIndexes, numbersToIndexes_obs;
	vector<int> localCubesIndexes;
	vector<std::vector<int>> listOfFcToBeReceived;
	vector<std::vector<int>> listOfFcToBeSent;
	ArrayXd thetas, phis;
	ArrayXd weightsThetas, weightsPhis;
	ArrayXXcd shiftingArrays;

	boost::multi_array<ArrayXcd, 3> alphaTranslations;
	boost::multi_array<ArrayXi, 3> alphaTranslationsIndexesNonZeros;
	boost::multi_array<ArrayXi, 3> alphaTranslationsIndexes;
	LagrangeFastInterpolator2D lfi2D;

	Array<ArrayXXcd, Dynamic, 1> Sdown, Sdown_obs;

	Level();
	Level(const int l,
	      const double leaf_side_length,
	      const double big_cube_lower_coord[3],
	      const ArrayXXd& cubes_centroids);

	Level(const int l,
	      const int N_expansion,
	      const double leaf_side_length,
	      const double big_cube_lower_coord[3],
	      const ArrayXXd& cubes_centroids,
	      const ArrayXXd& cubes_centroids_obs,
	      const complex<double>& wavenumber,
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
	      const int VERBOSE);

	void copyLevel(const Level& levelToCopy);
	Level(const Level&); // copy constructor
	Level& operator=(const Level&); // copy assignment operator

	Level(const Level& sonLevel,
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
	      const int VERBOSE); // fromSon constructor

	~Level();

	bool getLeaf(void) const { return leaf; }
	void setCeiling(bool l) { ceiling = l; }
	bool getCeiling(void) const { return ceiling; }
	int getLevel(void) const { return level; }
	// //! maxNumberCubes1D = static_cast<int> pow(2.0, level)
	int getMaxNumberCubes1D(void) const { return maxNumberCubes1D; }
	int getNumberTimesCopied(void) const { return numberTimesCopied; }
	void incrementNumberTimesCopied(void) { numberTimesCopied += 1; }
	void NCubesXYZComputation();
	int getNCubesX(void) const { return NCubesX; }
	int getNCubesY(void) const { return NCubesY; }
	int getNCubesZ(void) const { return NCubesZ; }
	int getOffsetAlphaIndexX(void) const { return offsetAlphaIndexX; }
	int getOffsetAlphaIndexY(void) const { return offsetAlphaIndexY; }
	int getOffsetAlphaIndexZ(void) const { return offsetAlphaIndexZ; }
	int getN(void) const { return N; }
	double getCubeSideLength(void) const { return cubeSideLength; }
	std::complex<double> getK(void) const { return k; }
	void addNode(Cube cube) { cubes.push_back(cube); }
	Cube getCube(const int i) const { return cubes[i]; }
	double getCubesSizeMB(void) const { return cubes.size() * (thetas.size() * phis.size()) * 2 * 2.0 * 4.0 / (1024.0 * 1024.0); }
	std::vector<Dictionary<int, int>> getNumbersToIndexes(void) const { return numbersToIndexes; }
	std::vector<int> getLocalCubesIndexes(void) const { return localCubesIndexes; }
	void computeOldIndexesOfCubes(ArrayXi& oldIndexesOfCubes,
	                              ArrayXi& oldIndexesOfCubes_obs);
	int getNumbersToIndexesSize(void) const { return numbersToIndexes.size(); }
	int getNumberToIndex(const int i) const { return numbersToIndexes[i].getKey(); }
	int getIndexToIndex(const int i) const { return numbersToIndexes[i].getVal(); }
	int getIndexOfNumber(const int) const;
	int getLevelSize(void) const { return cubes.size(); }

	int getSizeOfAlphaTransParticipantsIndexes(void) const
	{
		int result = 0;
		for (unsigned int i = 0; i < cubes.size(); ++i) result += cubes[i].AlphaTransParticipantsIndexes.size();
		return result;
	}

	double getSizeMBOfAlphaTransParticipantsIndexes(void) const
	{
		return getSizeOfAlphaTransParticipantsIndexes() * 4.0 / (1024.0 * 1024.0);
	}

	ArrayXd getThetas(void) const { return thetas; }
	ArrayXd getPhis(void) const { return phis; }
	int getNThetas(void) const { return thetas.size(); }
	int getNPhis(void) const { return phis.size(); }
	void computeGaussLocatedArguments
	(const ArrayXi& local_cubes_NRWG,
	 const ArrayXi& local_RWG_numbers,
	 const ArrayXi& local_RWG_Numbers_surf,
	 const ArrayXXi& local_RWGNumbers_signedTriangles,
	 const ArrayXXd& local_RWGNumbers_trianglesCoord,
	 const int N_Gauss);
	void computeObsLocatedArguments(const ArrayXi& cubes_Nobs,
	                                const ArrayXi& obs_numbers,
	                                const ArrayXXd& obsNumber_obsCoord);
	void RWGs_renumbering(void);
	void shiftingArraysComputation(void);
	double getShiftingArraysSizeMB(void) const { return shiftingArrays.size() * 2.0 * 4.0 / (1024.0 * 1024.0); }
	Array3i getAlphaTranslationsExtents(void) const;
	// blitz::Array< blitz::Array<std::complex<float>, 1>, 3> getAlphaTranslations(void) const {return alphaTranslations;}
	double getAlphaTranslationsSizeMB(void) const;
	double getAlphaTranslationsIndexesSizeMB(void) const { return alphaTranslationsIndexes.size() * 4.0 / (1024.0 * 1024.0); }
	ArrayXXcd getShiftingArrays(void) const { return shiftingArrays; }

	const ArrayXcd getShiftingArray(const double Dx, const double Dy, const double Dz) const
	{
		return shiftingArrays.row((Dx > 0.0) * 4 + (Dy > 0.0) * 2 + (Dz > 0.0));
	}

	void alphaTranslationsComputation(const double alphaTranslation_smoothing_factor,
	                                  const double alphaTranslation_thresholdRelValueMax,
	                                  const double alphaTranslation_RelativeCountAboveThreshold);
	void alphaTranslationIndexConstructionZ(ArrayXi& newAlphaIndex,
	                                        const ArrayXi& oldAlphaIndex,
	                                        const int alphaCartesianCoordZ,
	                                        const int N_theta,
	                                        const int N_phi);
	void alphaTranslationIndexConstructionY(ArrayXi& newAlphaIndex,
	                                        const ArrayXi& oldAlphaIndex,
	                                        const int alphaCartesianCoordY,
	                                        const int N_theta,
	                                        const int N_phi);
	void alphaTranslationIndexConstructionX(ArrayXi& newAlphaIndex,
	                                        const ArrayXi& oldAlphaIndex,
	                                        const int alphaCartesianCoordX,
	                                        const int N_theta,
	                                        const int N_phi);
	// blitz::Array<float, 1> getWeightsThetas(void) const {return weightsThetas;}
	// blitz::Array<float, 1> getWeightsPhis(void) const {return weightsPhis;}
	// LagrangeFastInterpolator2D getLfi2D(void) const {return lfi2D;}

	void sortCubesByParents(void);
	// void printCubesSonsIndexes(void);
	// void printCubesFathersNumbers(void);
	// void printCubesRCenters(void);
	// void printCubesRWG_numbers(void);
	void updateFatherIndexes(const Level& fatherLevel);
	void setNumbersToIndexes(void);
	void sortNumbersToIndexes(void);
	// void printNumbersToIndexes(void);
	void searchCubesNeighborsIndexes(void);
	// void computeLocalCubesIndexes(void);

	void computeSup(ArrayXXcd& Sup,
	                const complex<double>& k,
	                const ArrayXcd& I_PQ,
	                const Cube& cube,
	                const ArrayXd& thetas,
	                const ArrayXd& phis,
	                const ArrayXi&);
	void computeSup_obs(ArrayXXcd& Sup,
	                    const complex<double>& k,
	                    const ArrayXXcd& E_obs,
	                    const Cube& cube,
	                    const ArrayXd& thetas,
	                    const ArrayXd& phis,
	                    const bool);
	void sphericalIntegration(ArrayXcd& ZI,
	                          const ArrayXXcd& Sdown,
	                          const Cube& cube,
	                          const ArrayXd& thetas,
	                          const ArrayXd& phis,
	                          const double w,
	                          const complex<double>& mu_r,
	                          const complex<double>& k,
	                          const ArrayXcd& CFIE,
	                          const bool Mcurrent,
	                          const ArrayXi&);
	void sphericalIntegration_obs
	(ArrayXXcd& E_obs,
	 const ArrayXXcd& Sdown,
	 const Cube& cube,
	 const ArrayXd& thetas,
	 const ArrayXd& phis,
	 const double w,
	 const complex<double>& mu_r,
	 const complex<double>& k,
	 const ArrayXcd& CFIE,
	 const bool Mcurrent);
};

#endif
