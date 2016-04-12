/*
 * ht.cu
 *
 *  Created on: Feb 17, 2010
 *      Author: zhmurov
 */
#include "ht.cuh"

int hasGauss = 0;
float gauss;

struct HybTau{
	uint4* h_seeds;
	uint4* d_seeds;
	uint4 mseed;
};

 HybTau ht;
__device__ __constant__ HybTau c_ht;

int initRand(int seed, int Np){
	printf("Initializing Hybrid Taus PRNG...");
	ht.h_seeds = (uint4*)calloc(Np, sizeof(uint4));
	generateSeeds(ht.h_seeds, seed, Np);
	cudaMalloc((void**)&ht.d_seeds, Np * sizeof(uint4));
	cudaMemcpy(ht.d_seeds, ht.h_seeds, Np * sizeof(uint4), cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(c_ht, &ht, sizeof(HybTau), 0, cudaMemcpyHostToDevice);
	printf("done\n");
	return 1;
}

void generateSeeds(uint4* seeds, int seed, int Np){
	for(int i = 0; i < Np; i++){
		do{
			seeds[i].x = (unsigned)(ran2(&seed)*UINT_MAX);
		} while(seeds[i].x < 128);
		do{
			seeds[i].y = (unsigned)(ran2(&seed)*UINT_MAX);
		} while(seeds[i].y < 128);
		do{
			seeds[i].z = (unsigned)(ran2(&seed)*UINT_MAX);
		} while(seeds[i].z < 128);
		do{
			seeds[i].w = (unsigned)(ran2(&seed)*UINT_MAX);
		} while(seeds[i].w < 128);
	}
	ht.mseed = seeds[0];
}

// Random number generator.
// Taked from GPU Gems 3, Chapter 37

__device__  float uintToFloat(unsigned uint){
	unsigned itemp = jflone | (jflmsk & uint);
	float result = (*(float *)&itemp) - 1.0f;
	if(result == 0){
		return EPS;
	} else {
		return result;
	}
}

__device__  unsigned TausStep(unsigned &z, int S1, int S2, int S3, unsigned M){
	unsigned b = (((z << S1)^z) >> S2);
	return z = (((z & M) << S3) ^b);
}

__device__  unsigned LCGStep(unsigned &z, unsigned A_ht, unsigned C_ht){
	return z = (A_ht * z + C_ht);
}

__device__  unsigned HybridTaus(uint4 &seed){
	return 	TausStep(seed.x, 13, 19, 12, 4294967294UL) ^
		TausStep(seed.y,  2, 25,  4, 4294967288UL) ^
		TausStep(seed.z,  3, 11, 17, 4294967280UL) ^
		LCGStep(seed.w, 1664525, 1013904223UL);
	/*return 2.3283064365387e-10 * (
			TausStep(seed.x, 13, 19, 12, 4294967294UL) ^
			TausStep(seed.y,  2, 25,  4, 4294967288UL) ^
			TausStep(seed.z,  3, 11, 17, 4294967280UL) ^
			LCGStep(seed.w, 1664525, 1013904223UL));*/
}


__device__ float4 rforce(int d_i){
	uint4 seed = c_ht.d_seeds[d_i];
	float4 result;
	float r = sqrtf(-2.0f * logf(uintToFloat(HybridTaus(seed))));
	float theta = 2.0f*M_PI*uintToFloat(HybridTaus(seed));
	result.x = r*__sinf(theta);
	result.y = r*__cosf(theta);
	r = sqrtf(-2.0f * logf(uintToFloat(HybridTaus(seed))));
	theta = 2.0f*M_PI*uintToFloat(HybridTaus(seed));
	result.z = r*__sinf(theta);
	result.w = r*__cosf(theta);
	c_ht.d_seeds[d_i] = seed;
	return result;
}
