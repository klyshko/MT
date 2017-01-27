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
#include "Cuda.h"
 #include <cuda.h>

#define R_MT 8.12f //nm  microtube radius
#define r_mon 2.0f //nm	 monomer radius
#define lj_cutoff (6.0)
#define LARGENUMBER 999999

#define PF_NUMBER 13
#define ANG_THRES 1.0f
#define R_THRES (r_mon * 8)

#define Turn 13 //Monomers per turn(2PI) 

#define BLOCK_SIZE 32
#define ZERO 	999999

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
#define PAIR_CUTOFF 2.5f

typedef struct {
	real x;
	real y;
	real z;
	real fi;
	real theta;
	real psi;
	real w;
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
		bool* extra;
        int* mon_type;
        int* gtp;
        int* on_tubule_cur;
        int* on_tubule_prev;
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

struct Tea {
	float4 *rforce; // Precomputed random forces
	float4 *mforce; // Copy of molecular forces for each bead
	float4 *coords; // Copy of coordinates for each bead
	float *d_beta_ij; // $\beta_{ij}$ from eq. (14) in Geyer&Winter, 2009;
	float *h_beta_ij;
	float *d_epsilon; // Array of epsilon values for individual beads, evaluated on device and used to fill d_beta_ij
	float *h_epsilon;
	float4 *d_ci; // Does not actually store $C_i$, only $\sum (Dij/Dii)^2$
    float *d_tensor; // Memory for tensor when running in `exact` mode
	int epsilon_freq; // How often to update epsilon, beta_ij and c_i
	int Ntot; // Number of aminos per trajectory, nothing special
	float a; // Bead hydrodynamic radius, in A
	int capricious; // If != 0, then the simulation will stop if HI tensor has abnormal values. If zero, the simulation will continue anyway (and it is probably perfectly fine).
	int unlisted; // If ==0, then beads will interact hydrodynamically with their friends in covalent, native and pairs lists
    int exact; // If > 0, use Cholesky-based treatment
	float epsmax; // If epsilon exceeds this value, then abort simulation. Default: never [in capricious mode, epsilon > 1.0 will trigger stop anyway]
};

struct float6 {
	float3 xx; // xx, xy, xz
	float3 yz; // yy, yz, zz
};

#define _XX xx.x
#define _XY xx.y
#define _XZ xx.z
#define _YX xx.y
#define _YY yz.x
#define _YZ yz.y
#define _ZX xx.z
#define _ZY yz.y
#define _ZZ yz.z


#endif /* MT_H_ */
