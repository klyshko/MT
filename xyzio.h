/*
 * xyzio.h
 *
 *  Created on: Apr 7, 2011
 *      Author: zhmurov
 */

#pragma once

#define XYZ_FIELD_SIZE	16

typedef struct {
	char name;
	double x;
	double y;
	double z;
} XYZAtom;

typedef struct {
	int atomCount;
	XYZAtom* atoms;
} XYZ;

void readXYZ(const char* filename, XYZ* xyzData);
void readCoordinatesFromXYZ(const char* filename, double* x, double* y, double* z, int count);
void writeXYZ(const char* filename, XYZ* xyzData);
