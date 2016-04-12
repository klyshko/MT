/*
 * main.cpp
 *
 *  Created on: 11.04.2012
 *      Author: zhmurov

    Modified on: 02.01.2016
         Author: klyshko
 */
#ifdef USE_MPI
    #include <mpi.h>
#endif 
#include "globals.h"
#include "mt.h"
#include "wrapper.h"
#include "preparator.h"
#include "updater.h"

extern void compute(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies);

#ifdef USE_MPI
int mpi_dev_cur, mpi_rank, mpi_size;
#endif

int main(int argc, char *argv[]){

#ifdef USE_MPI
    MPI::Init(argc, argv);
    mpi_rank = MPI::COMM_WORLD.Get_rank();
    mpi_size = MPI::COMM_WORLD.Get_size();
#endif

    parseParametersFile(argv[1], argc, argv);

#ifdef USE_MPI
    addParameterT<int>(PARAMETER_MPI_RANK, mpi_rank);
    if (getYesNoParameter(PARAMETER_MPI_DEV_AUTO, 0)) {
        int mpi_dpn = getIntegerParameter(PARAMETER_MPI_DEVPERNODE);
        mpi_dev_cur = mpi_rank % mpi_dpn;
    } else if (mpi_size > 1) {
        mpi_dev_cur = getIntegerParameter(PARAMETER_MPI_DEVICE(mpi_rank));
    } else {
        mpi_dev_cur = getIntegerParameter(PARAMETER_MPI_DEVICE(mpi_rank),getIntegerParameter(PARAMETER_DEVICE));
    }
    addParameterT<int>(PARAMETER_MPI_DEV_CUR, mpi_dev_cur);

    if (getYesNoParameter(PARAMETER_MPI_FIRSTRUN, 1)) // We adjust the 'firstrun' parameter so the trajectory numbers are different for each thread
        par.firstrun = mpi_rank * par.Ntr + par.firstrun;

#endif

    initParameters(argc, argv); // Init parameters of model from conf files

#if defined(READ_FROM_DCD)
    if(argc < 4)
    {
        printf("Specify config file and dcds for both coordinates and angles filenames\n");
        exit(-1);
    }
    char* dcd_filename_xyz = argv[2];
    char* dcd_filename_ang = argv[3];
    ReadFromDCD(par, top, dcd_filename_xyz, dcd_filename_ang);
    saveCoordPDB("result_xyz.pdb", "result_ang.pdb");
    exit(0);
#endif

    initTimer();

    if (par.is_assembly){
        AssemblyInit();
    }

    if (par.out_energy){
        energies = (Energies*)malloc(par.Ntot * par.Ntr * sizeof(Energies));
    }
    
    compute(r, f, par, top, energies);                  // <---------- /// ---------------     COMPUTING of DYNAMICS

    saveCoordPDB("result_xyz.pdb", "result_ang.pdb");
    
#ifdef USE_MPI
        MPI::COMM_WORLD.Barrier();
        MPI::Finalize();
#endif


    free(r);
    free(f);
    free(top.lateral);
    free(top.lateralCount);
    free(top.longitudinal);
    free(top.longitudinalCount);
    free(top.harmonic);
    free(top.harmonicCount);
    free(top.LJ);
    free(top.LJCount);
    free(top.fixed);
    free(top.extra);
    if (par.out_energy){
        free(energies);
    }
    return 0;
}










