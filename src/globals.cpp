#include "globals.h"

Coord* r;
Coord* f;

Parameters par;
Topology top;
Energies* energies;

PDB pdb, pdb_ang;
DCD dcd;
XYZ xyz;

Coord* d_r;
Coord* d_f;

__device__ __constant__ Parameters c_par;
__device__ __constant__ Topology c_top;

Topology topGPU;
Energies* d_energies;

Tea tea;
__device__ __constant__ Tea c_tea;

