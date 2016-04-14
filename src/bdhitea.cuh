/*
 * bdhitea.cuh
 *
 *  Created on: Mar 4, 2013
 *	  Author: alekseenko
 */
#ifndef BDHITEA_CUH_
#define BDHITEA_CUH_
#include "Cuda.h"
#include "parameters.h"
#include "configreader.h"
#include "preparator.h"
#include "globals.h" 
#include "ht.cuh"

// TODO: Rename this module, since it not inly implements TEA-HI, but also exact (Cholesky-based) HI

void initTeaIntegrator();
void integrateTea();
void deleteTeaIntegrator();

void createTeaUpdater();
void initTeaUpdater();
void updateTea(long long int step);


// Macroses to facilitate access to float6 members


__global__ void integrateTea_prepare(Coord* d_f, Coord* d_r);
__device__ inline float6 integrateTea_RPY(const float4& dr);
__device__ inline float4 integrateTea_epsilon_local(const float4& coord1, const int idx2, Coord* d_r);
__global__ void integrateTea_epsilon_unlisted(Coord* d_r);
__device__ inline float4 integrateTea_force(const float4& coord1, const int idx2, const float3& ci, const int idx1);
__global__ void integrateTea_kernel_unlisted(Coord* d_f, Coord* d_r);


#endif /* BDHITEA_CUH_ */

// SOP-GPU with TEA --- Standard Operating Procedures for Generator Protection Unit with Technical and Economic Analysis

