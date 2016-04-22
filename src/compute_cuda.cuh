#include <cuda.h>
#include "wrapper.h"
#include "Cuda.h"
#include "mt.h"
#include "globals.h"
#include "parameters.h"
#include "updater.h"
#include "bdhitea.cuh"

__global__ void compute_kernel(const Coord* d_r, Coord* d_f);
__global__ void pairs_kernel(const Coord* d_r);
__global__ void energy_kernel(const Coord* d_r, Energies* d_energies);
__global__ void LJ_kernel(const Coord* r);
__global__ void integrate_kernel(Coord* d_r, Coord* d_f);


void compute(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies);
void deleteIntegration(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies);
void initIntegration(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies);

