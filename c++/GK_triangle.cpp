/**********************************************************************
 *
 * GK_triangle.cpp
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

#include <stdio.h>
#include <stdlib.h>
#include <sstream>
using namespace std;
#include "GK_triangle.hpp"

IT_points_weights::IT_points_weights(int Nquads_tmp)
{
	Nquads = Nquads_tmp;
	FILE *fp;
	stringstream filename;
	filename << "./c++/QuadratureTables/" << Nquads << ".txt";
	fp = fopen(filename.str().c_str(), "r");

	xi = new double[Nquads]; eta = new double[Nquads]; weights = new double[Nquads];

	if (fp == NULL) {
		printf("Error in GK_triangle.cpp: file cannot be opened. \n");
		exit(1);
	}

	double gamma;
	for (int i = 0; i<Nquads; i++) {
		fscanf(fp, "%lf %lf %lf %lf", &weights[i], &xi[i], &eta[i], &gamma);
	}

	fclose(fp);

}
