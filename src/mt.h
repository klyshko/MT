/*
 * mt.h
 *
 *  Created on: Apr 16, 2012
 *      Author: zhmurov
 */
#pragma once

#ifndef MT_H_
#define MT_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "pdbio.h"
#include "dcdio.h"
#include "xyzio.h"

#define R_MT 8.12f //nm  microtube radius
#define r_mon 2.0f //nm	 monomer radius
#define lj_cutoff (6.0)
#define LARGENUMBER 999999

#define PF_NUMBER 13
#define ANG_THRES 0.15f
#define R_THRES (r_mon * 0.8)

#define Turn 13 //Monomers per turn(2PI) 

#define real float
#define real4 float4

#define KB 0.0019872041f //kcal/(mol*K)

#define ANGLE_CUTOFF 1.0f

#define xp1_def (0.5f*R_MT*(cosf(2.0f*M_PI/13.0f) - 1.0f))
#define yp1_def (0.5f*R_MT*sinf(2.0f*M_PI/13.0f))
#define zp1_def (-3.0f*r_mon/13.0f)

#define xp2_def (0.5f*R_MT*(cosf(2.0f*M_PI/13.0f) - 1.0f))
#define yp2_def (-0.5f*R_MT*sinf(2.0f*M_PI/13.0f))
#define zp2_def (3.0f*r_mon/13.0f)

//numerical test
#define delta_x		0.001
#define delta_y		0.001
#define delta_z		0.001
#define delta_fi	0.001
#define delta_psi	0.001
#define delta_theta	0.001

//assembly_mode
#define PAIR_CUTOFF 2.0

typedef struct {
	real x;
	real y;
	real z;
	real fi;
	real theta;
	real psi;
} Coord;

typedef struct{
		int* harmonicCount;
		int* harmonic;
		int maxHarmonicPerMonomer;
		int* longitudinalCount;
		int* longitudinal;
		int maxLongitudinalPerMonomer;
		int* lateralCount;
		int* lateral;
		int maxLateralPerMonomer;
		int* LJCount;
        int* LJ;
		int maxLJPerMonomer;
		bool* fixed;
        int* mon_type;
} Topology;

typedef struct{
        double U_harm;
        double U_long;
        double U_lat;
        double U_psi;
        double U_fi;
        double U_teta;
        double U_lj;
} Energies;

#endif /* MT_H_ */
