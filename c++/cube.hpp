/**********************************************************************
 *
 * cube.h
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

#ifndef CUBE_H
#define CUBE_H

#include <iostream>
#include <boost/multi_array.hpp>
#include <Eigen/Dense>
#include <vector>
#include <algorithm>  // Include STL algorithms for sorting lists and vectors

using namespace std;
using namespace Eigen;

#include "EMConstants.hpp"
#include "GK_triangle.hpp"
#include "dictionary.hpp"

class Cube
{
	int number;
	int index;
	int oldIndex;
	int fatherNumber;
	int fatherIndex;
public:
	int procNumber;
	int fatherProcNumber;
	std::vector<int> sonsIndexes;
	std::vector<int> sonsProcNumbers;
	std::vector<int> neighborsIndexes;
	std::vector<int> AlphaTransParticipantsIndexes, AlphaTransParticipantsIndexes_obs;
	std::vector<int> nonLocalAlphaTransParticipantsIndexes;
	double rCenter[3];
	double absoluteCartesianCoord[3];

	std::vector<int> RWG_numbers;
	std::vector<int> RWG_numbers_surf;
	std::vector<int> Triangle_numberOfRWGs;
	std::vector<int> TriangleToRWGindex;
	std::vector<int> TriangleToSurf;
	std::vector<double> TriangleToRWGweight;
	std::vector<double> TriangleToRWG_ropp;
	ArrayXXd triangle_GaussCoord;
	ArrayXXd triangle_nHat;
	ArrayXi obs_numbers;
	ArrayXXd r_obs;
	// constructors
	Cube(void)
	{
	};
	Cube(const int level,
	     const double sideLength,
	     const double big_cube_lower_coord[3],
	     const double r_c[3]);
	//! allows a cube to be overwritten by the provided cube
	/*!
	  \param cubeToCopy the cube to copy
	*/
	void copyCube(const Cube& cubeToCopy);
	Cube(const Cube&); // copy constructor
	Cube& operator=(const Cube&); // copy assignment operator
	//! the constructor from a children cube
	/*!
	  \param sonCube the children cube
	  \param level the number of the level
	  \param big_cube_lower_coord the lower coordinate of the father of all cubes
	  \param sideLength the length of the cube side
	  \return the constructed cube
	*/
	Cube(const Cube& sonCube,
	     const int level,
	     const double big_cube_lower_coord[3],
	     const double sideLength);
	//! the destructor
	~Cube();

	// specific functions
	void addSon(const Cube&);
	//! returns the sonsIndexes vector
	void setIndex(const int i) { index = i; }
	int getIndex(void) const { return index; }
	void setOldIndex(const int i) { oldIndex = i; }
	int getOldIndex(void) const { return oldIndex; }
	int getNumber(void) const { return number; }
	int getProcNumber(void) const { return procNumber; }
	int getFatherNumber(void) const { return fatherNumber; }
	void setFatherNumber(const int n) { fatherNumber = n; }
	int getFatherIndex(void) const { return fatherIndex; }
	void setFatherIndex(const int i) { fatherIndex = i; }

	//! \brief computes the points locations and values for the arguments for the complex exponentials
	//! that will be used in computing the radiation function of the leaf cube
	/*!
	  \param target_mesh the mesh of the target, C++ object
	  \param N_Gauss int, the number of points per triangles
	  \return void
	*/
	void computeGaussLocatedArguments
	(const ArrayXi& local_RWG_numbers,
	 const ArrayXi& local_RWG_Numbers_surf,
	 const ArrayXXi& local_RWGNumbers_signedTriangles,
	 const ArrayXXd& local_RWGNumbers_trianglesCoord,
	 const int startIndex_in_localArrays,
	 const int NRWG,
	 const int N_Gauss);
	// overloaded operators
	bool operator==(const Cube&) const;
	bool operator<(const Cube&) const;
};

#endif
