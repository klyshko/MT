/*
 * dcdio.h
 *
 *  Created on: Feb 27, 2009
 *      Author: zhmurov
 */

#pragma once

#include <stdio.h>

typedef struct {
	int nfile, npriv, nsavc, nstep, N;
	float delta;
	char remark1[80];
	char remark2[80];
} DCDHeader;

typedef struct {
	float* X;
	float* Y;
	float* Z;
} DCDFrame;

typedef struct {
	double a, b, c;
	double alpha, beta, gamma;
} DCDUnitCell;

typedef struct {
	FILE* file;
	int hasUC;
	DCDUnitCell uc;
	DCDHeader header;
	DCDFrame frame;
} DCD;

void createDCD(DCD* dcd, int atomCount, int frameCount, int firstFrame, float timestep, int dcdFreq, int hazUC, int UCX, int UCY, int UCZ);

void dcdOpenWrite(DCD* dcd, char *dcd_filename);
void dcdOpenAppend(DCD* dcd, char *dcd_filename);
void dcdOpenRead(DCD* dcd, char *dcd_filename);

void dcdWriteHeader(DCD dcd);
void dcdWriteFrame(DCD dcd);
void dcdReadHeader(DCD* dcd);
int dcdReadFrame(DCD* dcd);

void dcdClose(DCD dcd);
