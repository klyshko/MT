
#pragma once

#include "ran2.h"
#include <stdio.h>
#define PRNGNAME "HybridTaus";

#define jflone 0x3f800000
#define jflmsk 0x007fffff
#define EPS 1.0e-8f


void generateSeeds(uint4* seeds, int seed, int Np);
int initRand(int seed, int Np);
__device__ float uintToFloat(unsigned uint);
__device__ unsigned TausStep(unsigned &z, int S1, int S2, int S3, unsigned M);
__device__ unsigned LCGStep(unsigned &z, unsigned A_ht, unsigned C_ht);
__device__ unsigned HybridTaus(uint4 &seed);
__device__ float4 rforce(int d_i);