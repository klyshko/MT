#pragma once

#include "mt.h"
#include "parameters.h"
#include "Cuda.h"

using namespace std;

extern Coord* r;
extern Coord* f;

extern Parameters par;
extern Topology top;
extern Energies* energies;

extern PDB pdb;
extern PDB coordspdb;
extern PDB pdb_ang;
extern DCD dcd;
extern XYZ xyz;

extern Coord* d_r;
extern Coord* d_f;

extern __device__ __constant__ Parameters c_par;
extern __device__ __constant__ Topology c_top;

extern Topology topGPU;
extern Energies* d_energies;

extern Tea tea;
extern __device__ __constant__ Tea c_tea;