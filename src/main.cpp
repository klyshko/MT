/*
 * main.cpp
 *
 *  Created on: 11.04.2012
 *      Author: zhmurov

    Modified on: 02.01.2016
         Author: klyshko
 */

#include "mt.h"
#include "parameters.h"
#include "configreader.h"
#include "timer.h"
#include <ctime>
#include "wrapper.h"
#include <stdio.h>

#ifdef USE_MPI
    #include <mpi.h>
#endif 

void read_PDB(const char* filename_xyz, const char* filename_ang);
void saveCoordPDB(const char* pdbfilename_xyz, const char* pdbfilename_ang);
void ReadFromDCD(Parameters par, Topology top, char* dcdfilename_xyz, char* dcdfilename_ang);
void saveCoordDCD();
void update(long long int step, int* mt_len);
int change_conc(int* delta, int* mt_len);
void UpdateLJPairs();
void UpdatePairs();
void AssemblyInit();
void writeRestart(long long int step);
void readRestart();
void OutputAllEnergies(long long int step);
void mt_length(long long int step, int* mt_len);
void OutputSumForce();
void OutputForces();
void initParameters(int argc, char* argv[]);

extern void compute(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies);

Parameters par;
Coord* r;
Coord* v;
Coord* f;

Topology top;
Energies* energies;

PDB pdb, pdb_ang;
DCD dcd;
XYZ xyz;

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
#endif

#ifdef USE_MPI
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

#if defined(ASSEMBLY)
    AssemblyInit();
#endif

#ifdef OUTPUT_EN
    energies = (Energies*)malloc(par.Ntot * par.Ntr * sizeof(Energies));
#endif

    compute(r, f, par, top, energies);
    saveCoordPDB("result_xyz.pdb", "result_ang.pdb");
    
#ifdef USE_MPI
        MPI::COMM_WORLD.Barrier();
        MPI::Finalize();
#endif


    free(r);
    free(f);
    free(v);
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
#ifdef OUTPUT_EN
    free(energies);
#endif
    return 0;
}

void update(long long int step, int* mt_len){
    printf("Saving coordinates at step %lld\n", step);
    printTime(step);
    printEstimatedTimeleft((real)step/(real)par.steps);
    saveCoordDCD();

#ifdef OUTPUT_FORCE
    OutputSumForce();
   OutputForces();
#endif

#ifdef OUTPUT_EN
    OutputAllEnergies(step);
#endif

}

#ifdef CONCENTRATION
int change_conc(int* delta, int* mt_len){
    int flag = 0;

    for(int tr = 0; tr < par.Ntr; tr++){

        int num_of_extra = 0;
        for (int i = 0; i < par.Ntot; i++){
            if (top.extra[i + tr * par.Ntot]){
                num_of_extra++;
            } 
        }
        float Vol = float(3.14 * par.rep_r * par.rep_r * par.zs[tr]);
        int NFreeDimers = (par.Ntot - mt_len[tr] - num_of_extra) / 2;

       
        while (1.0e7 * Nfree / 6.0 < par.conc * Vol){
            for(int i = 0; i < par.Ntot; i+=2){
                if (top.extra[i + tr * par.Ntot] && top.mon_type[i] == 0){
                    top.extra[i + tr * par.Ntot] = false;
                    top.extra[i + tr * par.Ntot+1] = false;
                    float x,y;
                    for(int j = 0; j > -1 ; j++) {
                        srand(j*time(NULL)-i);
                        x = par.rep_r - 2*(rand() % int(par.rep_r));
                        y = par.rep_r - 2*(rand() % int(par.rep_r));
                        if (x*x + y*y <= par.rep_r * par.rep_r){   //// check if its in cylinder
                            printf("New x,y coordinates for extra particle: %f  %f index: %d\n", x, y, i + tr * par.Ntot);
                            num_of_extra -= 2;
                            break;
                        }
                    }
                    float z = par.zs[tr] + par.rep_leftborder + 3*2*r_mon;
                    r[i + tr * par.Ntot].x = x; //
                    r[i + tr * par.Ntot].y = y;
                    r[i + tr * par.Ntot].z = z; ///fix
                    r[i + tr * par.Ntot +1].x = x;
                    r[i + tr * par.Ntot +1].y = y;
                    r[i + tr * par.Ntot +1].z = z + 2*r_mon;
                    flag++;
                    break;
                }

            }

            NFreeDimers += 1;
        }
            
           // par.zs[tr] += 2.0 * r_mon * delt / 13.0;
               
        //Vol = float(3.14 * par.rep_r * par.rep_r * par.zs[tr]);
        //Nfree = par.Ntot - mt_len[tr] - num_of_extra;
        printf("Concentration for tajectory[%d]: %f [muMole / L],\t %f [Dimers / nm^3],\t %d [Dimers / Volume],\t  Volume: %f [nm^3]\n", tr, 1.0e7 * NFreeDimers / (6.0 * Vol), NFreeDimers / Vol, NFreeDimers, Vol);
    }
    return flag;
}
#endif


void initParameters(int argc, char* argv[]){
#ifdef CUDA
#ifdef USE_MPI
    if (getYesNoParameter(PARAMETER_MPI_DEV_AUTO, 0)) {
        int mpi_dpn = getIntegerParameter(PARAMETER_MPI_DEVPERNODE);
        mpi_dev_cur = mpi_rank % mpi_dpn;
    } else if (mpi_size > 1) {
        mpi_dev_cur = getIntegerParameter(PARAMETER_MPI_DEVICE(mpi_rank));
    } else {
        mpi_dev_cur = getIntegerParameter(PARAMETER_MPI_DEVICE(mpi_rank),getIntegerParameter(PARAMETER_DEVICE));
    }
    par.device = mpi_dev_cur;
#else
    par.device = getIntegerParameter(PARAMETER_DEVICE);
#endif
#endif
    par.dt = getFloatParameter(PARAMETER_TIMESTEP);
    par.rseed = getIntegerParameter(PARAMETER_RANDOM_SEED);
    par.steps = getLongIntegerParameter(PARAMETER_NUMSTEPS,-1);
    par.stride = getLongIntegerParameter(PARAMETER_STRIDE,-1);

//#ifdef LJ_on
    par.ljpairscutoff = getFloatParameter(PARAMETER_LJPAIRSCUTOFF);
    par.ljpairsupdatefreq = getIntegerParameter(PARAMETER_LJPAIRSUPDATEFREQ);
//#endif

    getMaskedParameter(par.coordFilename_ang, PARAMETER_COORD_FILE_ANG);
    getMaskedParameter(par.coordFilename_xyz, PARAMETER_COORD_FILE_XYZ);
    getMaskedParameter(par.ffFilename, PARAMETER_FORCEFIELD_FILE);
    getMaskedParameter(par.condFilename, PARAMETER_CONDITIONS_FILE);
    par.fix = getIntegerParameter(PARAMETER_FIX,1);
    par.Ntr = getIntegerParameter(PARAMETER_RUNNUM,1);

    getMaskedParameter(par.restartkey, PARAMETER_RESTARTKEY);

    if(getYesNoParameter(PARAMETER_ISRESTART,0))
        par.is_restart = true;
    else
        par.is_restart = false;
    if(!par.is_restart)
    {
        par.firststep = getLongIntegerParameter(PARAMETER_FIRSTSTEP, 0);
    }
    else
    {
        FILE *keyf = safe_fopen(par.restartkey, "r");
        if (fscanf(keyf,"%lld",&(par.firststep)) != 1)
            DIE("Reading restartkey %s: unable to get firststep", &(par.restartkey));
        fclose(keyf);
    }   
    read_PDB(par.coordFilename_xyz, par.coordFilename_ang);                                 ////// <----- read topology from pdb files

    par.dcdFilename_xyz = (char**)calloc(par.Ntr, sizeof(char*));
    par.dcdFilename_ang = (char**)calloc(par.Ntr, sizeof(char*));
    par.restart_xyzFilename = (char**)calloc(par.Ntr, sizeof(char*));
    par.restart_angFilename = (char**)calloc(par.Ntr, sizeof(char*));

#if not defined(READ_FROM_DCD)
    createDCD(&dcd, par.Ntot, par.steps/par.stride, 1, par.dt, par.stride, 0, 0, 0, 0);
    char trajnum[10];
    for(int traj = 0; traj < par.Ntr; traj++){
            sprintf(trajnum, "%d", traj);
            par.dcdFilename_xyz[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
            par.dcdFilename_ang[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
            getMaskedParameterWithReplacement(par.dcdFilename_ang[traj], PARAMETER_DCD_FILE_ANG, trajnum, "<run>");
            getMaskedParameterWithReplacement(par.dcdFilename_xyz[traj], PARAMETER_DCD_FILE_XYZ, trajnum, "<run>");
            dcdOpenWrite(&dcd, par.dcdFilename_xyz[traj]);
            dcdWriteHeader(dcd);
            dcdClose(dcd);
            dcdOpenWrite(&dcd, par.dcdFilename_ang[traj]);
            dcdWriteHeader(dcd);
            dcdClose(dcd);

            par.restart_xyzFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
            par.restart_angFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
            getMaskedParameterWithReplacement(par.restart_xyzFilename[traj], PARAMETER_RESTART_XYZ_FILE, trajnum, "<run>");
            getMaskedParameterWithReplacement(par.restart_angFilename[traj], PARAMETER_RESTART_ANG_FILE, trajnum, "<run>");

    }
#endif
    if(par.is_restart)
    {
        readRestart();
    }
    parseParametersFile(par.ffFilename, argc, argv);

#if defined(MORSE)
    par.A_long = getFloatParameter(A_LONG);
    par.A_lat  = getFloatParameter(A_LAT);
    par.D_long  = getFloatParameter(D_LONG);
    par.D_lat   = getFloatParameter(D_LAT);
#else
    /*
    par.C = getFloatParameter(PARAMETER_HARMONIC_C);
    par.A_long = getFloatParameter(PARAMETER_LONGITUDINAL_A);
    par.b_long = getFloatParameter(PARAMETER_LONGITUDINAL_B);
    par.c_long = getFloatParameter(PARAMETER_LONGITUDINAL_C);
    par.r0_long = getFloatParameter(PARAMETER_LONGITUDINAL_R_0);
    par.d_long = getFloatParameter(PARAMETER_LONGITUDINAL_D);
    par.A_lat = getFloatParameter(PARAMETER_LATERAL_A);
    par.b_lat = getFloatParameter(PARAMETER_LATERAL_B);
    par.c_lat = getFloatParameter(PARAMETER_LATERAL_C);
    par.r0_lat = getFloatParameter(PARAMETER_LATERAL_R_0);
    par.d_lat = getFloatParameter(PARAMETER_LATERAL_D);
    */
    printf("No long/lat potentials. Continue to compute\n");
    //exit(0);
#endif
    
#if defined(BARR)
    par.a_barr_long = getFloatParameter(A_BARR_LONG);
    par.r_barr_long = getFloatParameter(R_BARR_LONG);
    par.w_barr_long = getFloatParameter(W_BARR_LONG)/( 2*(sqrt(-2*log(0.5))) );
    par.a_barr_lat = getFloatParameter(A_BARR_LAT);
    par.r_barr_lat = getFloatParameter(R_BARR_LAT);
    par.w_barr_lat = getFloatParameter(W_BARR_LAT)/( 2*(sqrt(-2*log(0.5))) );
#endif

 //ANGLES_HARMONIC  
    par.C = getFloatParameter(PARAMETER_HARMONIC_C);

    par.B_fi = getFloatParameter(PARAMETER_BENDING_B_FI);
    par.B_psi = getFloatParameter(PARAMETER_BENDING_B_PSI);
    par.B_theta = getFloatParameter(PARAMETER_BENDING_B_THETA);
    par.fi_0 = getFloatParameter(PARAMETER_BENDING_FI_0);
    par.psi_0 = getFloatParameter(PARAMETER_BENDING_PSI_0);
    par.theta_0 = getFloatParameter(PARAMETER_BENDING_THETA_0);

#ifdef LJ_on
    par.ljscale = getFloatParameter(PARAMETER_LJSCALE);
    par.ljsigma6 = pow(getFloatParameter(PARAMETER_LJSIGMA),6);
#endif

#if defined(REPULSIVE)
    par.rep_leftborder = getFloatParameter(REP_LEFTBORDER);
    par.rep_h = getFloatParameter(REP_H);
    par.rep_r = getFloatParameter(REP_R);
    par.rep_eps = getFloatParameter(REP_EPS); 
    for (int i = 0; i < par.Ntr; i++){
        par.zs[i] = par.rep_h;
    }   
#endif

   

    parseParametersFile(par.condFilename, argc, argv);
    par.Temp = getFloatParameter(PARAMETER_TEMPERATURE);
    par.gammaR = getFloatParameter(PARAMETER_GAMMA_R);
    par.gammaTheta = getFloatParameter(PARAMETER_GAMMA_THETA);
    par.varR = sqrtf(2.0f*KB*par.Temp*par.dt/par.gammaR);
    par.varTheta = sqrtf(2.0f*KB*par.Temp*par.dt/par.gammaTheta);
    par.alpha = getFloatParameter(ALPHA_GEOMETRY);
    par.freeze_temp = getFloatParameter(FREEZE_TEMP);
#ifdef CONCENTRATION
    par.conc = getFloatParameter(CONC);
#endif 
}

void read_PDB(const char* filename_xyz, const char* filename_ang){
    int i, j;
    readPDB(filename_xyz, &pdb);
    printf("Building topology....\n");
    
    par.Ntot = pdb.atomCount;
    
    r = (Coord*)calloc(par.Ntot*par.Ntr, sizeof(Coord));
    f = (Coord*)calloc(par.Ntot*par.Ntr, sizeof(Coord));

    top.mon_type = (int*)calloc(par.Ntot, sizeof(int));
    for(i = 0; i < par.Ntot; i++){
        if(pdb.atoms[i].name[1] == 'A')
            top.mon_type[i] = 0;
        else if(pdb.atoms[i].name[1] == 'B')
            top.mon_type[i] = 1;
    }

     //initialising fixed atoms list
    top.fixed = (bool*)calloc(par.Ntot, sizeof(bool));
    for(i = 0; i < par.Ntot; i++){
            if( pdb.atoms[i].resid <= par.fix)
                    top.fixed[i] = true;
            else
                    top.fixed[i] = false;
    }

     //initialising extra atoms list
    top.extra = (bool*)calloc(par.Ntot * par.Ntr, sizeof(bool));

    for(int traj = 0; traj < par.Ntr; traj++){
        for(i = 0; i < par.Ntot; i++){
            if(pdb.atoms[i].chain == 'X')
                    top.extra[i + traj*par.Ntot] = true;
            else
                    top.extra[i + traj*par.Ntot] = false;
        }
    }

    for(int traj = 0; traj < par.Ntr; traj++){
            for(i = 0; i < par.Ntot; i++){
                r[i+traj*par.Ntot].x = pdb.atoms[i].x;
                r[i+traj*par.Ntot].y = pdb.atoms[i].y;
                r[i+traj*par.Ntot].z = pdb.atoms[i].z;
            }
    }
    readPDB(filename_ang, &pdb_ang);
    for(int traj = 0; traj < par.Ntr; traj++){
            for(i = 0; i < par.Ntot; i++){
                r[i+traj*par.Ntot].fi = pdb_ang.atoms[i].x;
                r[i+traj*par.Ntot].psi = pdb_ang.atoms[i].y;
                r[i+traj*par.Ntot].theta = pdb_ang.atoms[i].z;
            }
    }

#ifdef LJ_on
    top.LJCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
#endif
    top.harmonicCount = (int*)calloc(par.Ntot, sizeof(int));
    top.longitudinalCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
    top.lateralCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
    
    // Covalent bonds map
    for(i = 0; i < par.Ntot; i++) top.harmonicCount[i] = 0;
    for(i = 0; i < par.Ntot; i++){
        for(j = 0; j < par.Ntot; j++){
            if(pdb.atoms[i].resid == pdb.atoms[j].resid &&
                    pdb.atoms[i].chain == pdb.atoms[j].chain
                    && (i!=j) && (abs(i-j) < 2) ){
                top.harmonicCount[i]++;
                if(top.maxHarmonicPerMonomer < top.harmonicCount[i])
                        top.maxHarmonicPerMonomer = top.harmonicCount[i];
            }
        }
    }
    top.harmonic = (int*)calloc(top.maxHarmonicPerMonomer*par.Ntot, sizeof(int));
    for(i = 0; i < par.Ntot; i++) top.harmonicCount[i] = 0;
    for(i = 0; i < top.maxHarmonicPerMonomer*par.Ntot; i++) top.harmonic[i] = -1;
    for(i = 0; i < par.Ntot; i++){
        for(j = 0; j < par.Ntot; j++){
            if(pdb.atoms[i].resid == pdb.atoms[j].resid &&
                    pdb.atoms[i].chain == pdb.atoms[j].chain &&
                    (i!=j)  && (abs(i-j) < 2)){
                if(pdb.atoms[i].id > pdb.atoms[j].id){
                    top.harmonic[top.maxHarmonicPerMonomer*i + top.harmonicCount[i]] = j; 
                }
                else{
                    top.harmonic[top.maxHarmonicPerMonomer*i + top.harmonicCount[i]] = -j; 
                }
                top.harmonicCount[i]++;
            }
        }
    }

#ifndef ASSEMBLY    
    // Longitudal exponential
    for(i = 0; i < par.Ntot * par.Ntr; i++) top.longitudinalCount[i] = 0;

    for (int traj = 0; traj < par.Ntr; traj++){
        for(i = 0; i < par.Ntot; i++){
            if (top.extra[i + traj * par.Ntot]) { //////////////////////////////////////////////////
                continue;
            }
            for(j = 0; j < par.Ntot; j++){
                if( (pdb.atoms[i].resid == (pdb.atoms[j].resid + 1)) &&
                    (pdb.atoms[i].chain == pdb.atoms[j].chain) &&
                    (pdb.atoms[i].id == (pdb.atoms[j].id + 1)) && (i != j) ){
                    top.longitudinalCount[i + par.Ntot * traj]++;
                    if(top.maxLongitudinalPerMonomer < top.longitudinalCount[i + par.Ntot * traj])
                            top.maxLongitudinalPerMonomer = top.longitudinalCount[i + par.Ntot * traj];
                }
            }
        }
    }
    top.longitudinal = (int*)calloc(top.maxLongitudinalPerMonomer * par.Ntot * par.Ntr, sizeof(int));
    
    for(i = 0; i < par.Ntot * par.Ntr; i++) top.longitudinalCount[i] = 0;

    for(int traj = 0; traj < par.Ntr; traj++){
        for(i = 0; i < par.Ntot; i++){
            if (top.extra[i + traj * par.Ntot]) {
                continue;                              /////////////////////////////////////////////////
            }
            for(j = 0; j < par.Ntot; j++){
                if(abs(pdb.atoms[i].resid - pdb.atoms[j].resid) == 1 &&
                    pdb.atoms[i].chain == pdb.atoms[j].chain &&
                    abs(pdb.atoms[i].id -  pdb.atoms[j].id) == 1){
                        if(pdb.atoms[i].id > pdb.atoms[j].id){
                            top.longitudinal[top.maxLongitudinalPerMonomer * par.Ntot * traj + i * top.maxLongitudinalPerMonomer + top.longitudinalCount[i + par.Ntot * traj]] = j;
                            //top.longitudinal[top.maxLongitudinalPerMonomer*i + top.longitudinalCount[i]] = j; // << BUG: indexation
                        }
                        else{
                            top.longitudinal[top.maxLongitudinalPerMonomer * par.Ntot * traj + i * top.maxLongitudinalPerMonomer + top.longitudinalCount[i + par.Ntot * traj]] = -j;
                            //top.longitudinal[top.maxLongitudinalPerMonomer*i + top.longitudinalCount[i]] = -j; /* BUG: same here */
                        }
                        top.longitudinalCount[i + par.Ntot * traj]++;
                }
            }
        }

    }
    
    // Lateral exponential
    for(i = 0; i < par.Ntot * par.Ntr; i++) top.lateralCount[i] = 0;
    
    real dr, xi, xj, yi, yj, zi, zj;
    real cos_fii, cos_fij, sin_fii, sin_fij,
          cos_psii, cos_psij, sin_psii, sin_psij,
          cos_thetai, cos_thetaj, sin_thetai, sin_thetaj;
    real xp1 = xp1_def;
    real yp1 = yp1_def;
    real zp1 = zp1_def;
    real xp2 = xp2_def;
    real yp2 = yp2_def;
    real zp2 = zp2_def;
    for(int traj = 0; traj < par.Ntr; traj++){

        for(i = 0; i < par.Ntot; i++){
            if (top.extra[i + traj * par.Ntot]) {
                continue;
            }
            xi = r[i].x;    
            yi = r[i].y;    
            zi = r[i].z;
            sin_fii = sinf(r[i].fi);
            cos_fii = cosf(r[i].fi);
            sin_psii = sinf(r[i].psi);
            cos_psii = cosf(r[i].psi);      
            sin_thetai = sinf(r[i].theta);
            cos_thetai = cosf(r[i].theta);  
            for(j = 0; j < par.Ntot; j++){
                xj = r[j].x;    
                yj = r[j].y;    
                zj = r[j].z;
                sin_fij = sinf(r[j].fi);
                cos_fij = cosf(r[j].fi);
                sin_psij = sinf(r[j].psi);
                cos_psij = cosf(r[j].psi);      
                sin_thetaj = sinf(r[j].theta);
                cos_thetaj = cosf(r[j].theta);  
                for(int ind = 0; ind < 2; ind++)
                {
                    if(ind == 0)
                    {
                        xp1 = xp2_def;
                        yp1 = yp2_def;
                        zp1 = zp2_def;
                        xp2 = xp1_def;
                        yp2 = yp1_def;
                        zp2 = zp1_def;
                    }
                    else
                    {
                        xp1 = xp1_def;
                        yp1 = yp1_def;
                        zp1 = zp1_def;
                        xp2 = xp2_def;
                        yp2 = yp2_def;
                        zp2 = zp2_def;
                    }

                        
                    dr = sqrtf(pow(zi - zj + zp2 * cos_fii * cos_thetai -
                        zp1 * cos_fij * cos_thetaj + yp2 * cos_thetai * sin_fii -
                        yp1 * cos_thetaj * sin_fij - xp2 * sin_thetai + xp1 * sin_thetaj,2) +
                        pow(xi - xj - yp2 * cos_fii * sin_psii + zp2 * sin_fii * sin_psii +
                        yp1 * cos_fij * sin_psij - zp1 * sin_fij * sin_psij +
                        cos_psii * (xp2 * cos_thetai + zp2 * cos_fii * sin_thetai +
                        yp2 * sin_fii * sin_thetai) - cos_psij * (xp1 * cos_thetaj + zp1 * cos_fij * sin_thetaj +
                        yp1 * sin_fij * sin_thetaj),2) + pow(yi - yj - zp2 * cos_psii * sin_fii + zp1 * cos_psij * sin_fij +
                        xp2 * cos_thetai * sin_psii - xp1 * cos_thetaj * sin_psij + yp2 * sin_fii * sin_psii * sin_thetai +
                        cos_fii * (yp2 * cos_psii + zp2 * sin_psii * sin_thetai) -
                        yp1 * sin_fij * sin_psij * sin_thetaj - cos_fij * (yp1 * cos_psij +
                        zp1 * sin_psij * sin_thetaj),2));
                    
                    if(dr < r_mon){
                        top.lateralCount[i + traj * par.Ntot]++;
                        if(top.maxLateralPerMonomer < top.lateralCount[i + traj * par.Ntot])
                                top.maxLateralPerMonomer = top.lateralCount[i + traj * par.Ntot];
                    }
                }
            }   
        }
    }
    
    top.lateral = (int*)calloc(top.maxLateralPerMonomer * par.Ntot * par.Ntr, sizeof(int));
    for(i = 0; i < par.Ntot * par.Ntr; i++) top.lateralCount[i] = 0;

    for(int traj = 0; traj < par.Ntr; traj++){
        for(i = 0; i < par.Ntot; i++){
             if (top.extra[i + traj * par.Ntot]) {
                continue;                                  /////////////////////////////////////////////////
            }
            xi = r[i].x;    
            yi = r[i].y;    
            zi = r[i].z;
            sin_fii = sinf(r[i].fi);
            cos_fii = cosf(r[i].fi);
            sin_psii = sinf(r[i].psi);
            cos_psii = cosf(r[i].psi);      
            sin_thetai = sinf(r[i].theta);
            cos_thetai = cosf(r[i].theta);  

            for(j = 0; j < par.Ntot; j++){
                xj = r[j].x;    
                yj = r[j].y;    
                zj = r[j].z;
                sin_fij = sinf(r[j].fi);
                cos_fij = cosf(r[j].fi);
                sin_psij = sinf(r[j].psi);
                cos_psij = cosf(r[j].psi);      
                sin_thetaj = sinf(r[j].theta);
                cos_thetaj = cosf(r[j].theta);  
                for(int ind = 0; ind < 2; ind++)
                {   
                    if(ind == 0)
                    {
                        xp1 = xp2_def;
                        yp1 = yp2_def;
                        zp1 = zp2_def;
                        xp2 = xp1_def;
                        yp2 = yp1_def;
                        zp2 = zp1_def;
                    }
                    else
                    {
                        xp1 = xp1_def;
                        yp1 = yp1_def;
                        zp1 = zp1_def;
                        xp2 = xp2_def;
                        yp2 = yp2_def;
                        zp2 = zp2_def;
                    }

                    dr = sqrtf(pow(zi - zj + zp2 * cos_fii * cos_thetai -
                    zp1 * cos_fij * cos_thetaj + yp2 * cos_thetai * sin_fii -
                    yp1 * cos_thetaj * sin_fij - xp2 * sin_thetai + xp1 * sin_thetaj,2) +
                    pow(xi - xj - yp2 * cos_fii * sin_psii + zp2 * sin_fii * sin_psii +
                    yp1 * cos_fij * sin_psij - zp1 * sin_fij * sin_psij +
                    cos_psii * (xp2 * cos_thetai + zp2 * cos_fii * sin_thetai +
                    yp2 * sin_fii * sin_thetai) - cos_psij * (xp1 * cos_thetaj + zp1 * cos_fij * sin_thetaj +
                    yp1 * sin_fij * sin_thetaj),2) + pow(yi - yj - zp2 * cos_psii * sin_fii + zp1 * cos_psij * sin_fij +
                    xp2 * cos_thetai * sin_psii - xp1 * cos_thetaj * sin_psij + yp2 * sin_fii * sin_psii * sin_thetai +
                    cos_fii * (yp2 * cos_psii + zp2 * sin_psii * sin_thetai) -
                    yp1 * sin_fij * sin_psij * sin_thetaj - cos_fij * (yp1 * cos_psij +
                    zp1 * sin_psij * sin_thetaj),2));
                    if(dr < r_mon){
                        if(ind == 1){
                            top.lateral[top.maxLateralPerMonomer * par.Ntot * traj + i * top.maxLateralPerMonomer + top.lateralCount[i + par.Ntot * traj]] = j;
                        }else{
                            top.lateral[top.maxLateralPerMonomer * par.Ntot * traj + i * top.maxLateralPerMonomer + top.lateralCount[i + par.Ntot * traj]] = -j;
                        }
                        top.lateralCount[i + par.Ntot * traj]++;
                    }
                }
            }
        }
    }

#endif
    
    printf("done building topology without LJ.\n");
}

void AssemblyInit()
{
    free(top.lateral);
    free(top.longitudinal);
    top.maxLateralPerMonomer = 16;
    top.maxLongitudinalPerMonomer = 8;

    top.lateral = (int*)calloc(top.maxLateralPerMonomer * par.Ntot * par.Ntr, sizeof(int));
    top.longitudinal = (int*)calloc(top.maxLongitudinalPerMonomer * par.Ntot * par.Ntr, sizeof(int));
    for(int i = 0; i < par.Ntot * par.Ntr; i++)
    {
        top.longitudinalCount[i] = 0;
        top.lateralCount[i] = 0;

    }
}

void mt_length(long long int step, int* mt_len){

    for (int traj = 0; traj < par.Ntr; traj++){
        int sum = 0;

        int protofilamentsLength[PF_NUMBER] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
        for (int i = traj * par.Ntot; i < (traj + 1) * par.Ntot; i++) {

            float rad = sqrt(r[i].x * r[i].x + r[i].y * r[i].y);
            if ((rad < R_MT + R_THRES) && (rad > R_MT - R_THRES)) {
               
            // detect particles inside microtubule 
                if (cosf(r[i].theta) > cosf(ANG_THRES)){
                    sum++;
                }

            /*  
                for (int j = 0; j < PF_NUMBER; j++){

                    float angle = j * 2 * M_PI / PF_NUMBER;
                    float psi = 0.0;

                    if (r[i].psi >= -ANG_THRES) {
                        psi = r[i].psi;
                    }
                    else {
                        psi = 2 * M_PI + r[i].psi;
                    }

                    if (psi < angle + ANG_THRES && psi > angle - ANG_THRES ) {

                        if ( (r[i].theta < ANG_THRES && r[i].theta > -ANG_THRES) 
                           ||    (( 2 * M_PI - fabs(r[i].theta)) < ANG_THRES &&  (2 * M_PI - fabs(r[i].theta)) > -ANG_THRES) ) {
                            if (r[i].fi < ANG_THRES && r[i].fi > - ANG_THRES){
                                protofilamentsLength[j]++;
                                //count++;
                                //printf("%d\t%f\t%f\t%f\t\t%f\t%f\t%f\t\t%d%s\t%d\t%d\n", j, r[i].x, r[i].y,r[i].z,r[i].fi,r[i].psi,r[i].theta, pdb.atoms[i].resid, pdb.atoms[i].name, count, count2);
                                break;
                            }
                            
                        }
                    }

                }
                */
            }
        }
        /*
        int sum = 0;
        for (int i = 0; i < PF_NUMBER; i++){
            //if (protofilamentsLength[i] % 2) protofilamentsLength[i]--;
            //printf("PF[%d] = %d\n", i, protofilamentsLength[i]);
            sum += protofilamentsLength[i];

        }*/
        mt_len[traj] = sum;
        printf("tubule[%d]: %d\n", traj, sum);
        /*char fileName[64];
        sprintf(fileName, "mtlength%d.dat", traj);
        */
    }
    FILE* mtLenFile = fopen("mt_len.dat", "a");
    fprintf(mtLenFile, "%lld\t" , step);
    for (int traj = 0; traj < par.Ntr; traj++){
        fprintf(mtLenFile, "%f\t" , 4.0 * (float)mt_len[traj] / PF_NUMBER);
    }
    fprintf(mtLenFile, "\n");
    fclose(mtLenFile);
        
}



void OutputAllEnergies(long long int step){

    Energies* fullEnergy = (Energies*)malloc(par.Ntr * sizeof(Energies));

    for (int tr = 0; tr < par.Ntr; tr++){
        fullEnergy[tr].U_harm = 0;
        fullEnergy[tr].U_long = 0;
        fullEnergy[tr].U_lat = 0;
        fullEnergy[tr].U_psi = 0;
        fullEnergy[tr].U_fi = 0;
        fullEnergy[tr].U_teta = 0;
        fullEnergy[tr].U_lj = 0;

        for (int i = 0; i < par.Ntot; i++){
            fullEnergy[tr].U_harm += energies[i + par.Ntot * tr].U_harm;
            fullEnergy[tr].U_long += energies[i + par.Ntot * tr].U_long;
            fullEnergy[tr].U_lat += energies[i + par.Ntot * tr].U_lat;
            fullEnergy[tr].U_psi += energies[i + par.Ntot * tr].U_psi;
            fullEnergy[tr].U_fi += energies[i + par.Ntot * tr].U_fi;
            fullEnergy[tr].U_teta += energies[i + par.Ntot * tr].U_teta;
            fullEnergy[tr].U_lj += energies[i + par.Ntot * tr].U_lj;
        }

        if (!isfinite(fullEnergy[tr].U_harm) || !isfinite(fullEnergy[tr].U_long) || !isfinite(fullEnergy[tr].U_lat) || !isfinite(fullEnergy[tr].U_lj)) {
            printf("Some energy in %d trajectory is NaN. NOT Exit program\n", tr);
            //exit(0);
        }

        //char fileName[64];
        //sprintf(fileName, "energies%d.dat", tr);
        //FILE* energyFile = fopen(fileName, "a");
        printf("Energies[%d]:\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", tr, 
            fullEnergy[tr].U_harm, fullEnergy[tr].U_long, fullEnergy[tr].U_lat, fullEnergy[tr].U_psi, fullEnergy[tr].U_fi, fullEnergy[tr].U_teta, fullEnergy[tr].U_lj);

        //fprintf(energyFile, "%lld\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\n",
        // step, fullEnergy[tr].U_harm, fullEnergy[tr].U_long, fullEnergy[tr].U_lat, fullEnergy[tr].U_psi, fullEnergy[tr].U_fi, fullEnergy[tr].U_teta, fullEnergy[tr].U_lj);
        //fclose(energyFile);
    }
    
}

void OutputSumForce(){
    for(int tr = 0; tr < par.Ntr; tr++){
        Coord f_sum;
        f_sum.x = 0; f_sum.y = 0; f_sum.z = 0;
        for(int i = tr * par.Ntot; i < (tr + 1) * par.Ntot; i++){
            f_sum.x += f[i].x;
            f_sum.y += f[i].y;
            f_sum.z += f[i].z;
        }
        printf("SF for %d traj %.16f\n", tr, sqrt(f_sum.x * f_sum.x + f_sum.y * f_sum.y + f_sum.z * f_sum.z));
        //printf("SF for %d traj: %f %f %f \n", tr, f_sum.x, f_sum.y, f_sum.z);// * f_sum.x + f_sum.y * f_sum.y + f_sum.z * f_sum.z));
    } 
}

void OutputForces(){

    for (int i = 0; i < par.Ntot * par.Ntr; i++){
        //printf("Force[%d].x = %f, y = %f, z = %f\n", i, f[i].x, f[i].y, f[i].z);
        //if (f[i].theta)
        printf("Force[%d].theta = %f, fi = %f, psi = %f\n", i, f[i].theta, f[i].fi, f[i].psi);
        printf("Angle[%d].theta = %f, fi = %f, psi = %f\n", i, r[i].theta, r[i].fi, r[i].psi);
    }
    
}



void saveCoordPDB(const char* pdbfilename_xyz, const char* pdbfilename_ang){
    int i;
    for(i = 0; i < par.Ntot; i++){
        pdb.atoms[i].x = r[i].x;
        pdb.atoms[i].y = r[i].y;
        pdb.atoms[i].z = r[i].z;
    }
    writePDB(pdbfilename_xyz, &pdb);
    for(i = 0; i < par.Ntot; i++){
        pdb.atoms[i].x = r[i].fi;
        pdb.atoms[i].y = r[i].psi;
        pdb.atoms[i].z = r[i].theta;
    }
    writePDB(pdbfilename_ang, &pdb);
}

void saveCoordDCD(){
    int i;
    for(int j=0; j < par.Ntr; j++){
        for(i = 0; i < par.Ntot; i++){
            dcd.frame.X[i] = r[i+par.Ntot*j].x;
            dcd.frame.Y[i] = r[i+par.Ntot*j].y;
            dcd.frame.Z[i] = r[i+par.Ntot*j].z;
        }
        dcdOpenAppend(&dcd, par.dcdFilename_xyz[j]);
        dcdWriteFrame(dcd);
        dcdClose(dcd);
        for(i = 0; i < par.Ntot; i++){
            dcd.frame.X[i] = r[i+par.Ntot*j].fi;
            dcd.frame.Y[i] = r[i+par.Ntot*j].psi;
            dcd.frame.Z[i] = r[i+par.Ntot*j].theta;
        }
        dcdOpenAppend(&dcd, par.dcdFilename_ang[j]);
        dcdWriteFrame(dcd);
        dcdClose(dcd);
    }
}



///////////////// Very rarely used functions 


void ReadFromDCD(Parameters par, Topology top, char* dcdfilename_xyz, char* dcdfilename_ang)
{

    DCD dcd_xyz, dcd_ang;
    dcd_xyz.frame.X = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_xyz.frame.Y = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_xyz.frame.Z = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_ang.frame.X = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_ang.frame.Y = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_ang.frame.Z = (float*)calloc(pdb.atomCount, sizeof(float));
    
    dcdOpenRead(&dcd_xyz, dcdfilename_xyz);
    dcdOpenRead(&dcd_ang, dcdfilename_ang);
    dcdReadHeader(&dcd_xyz);
    dcdReadHeader(&dcd_ang);
    while( ( dcdReadFrame(&dcd_xyz) == 0 ) && ( dcdReadFrame(&dcd_ang) == 0 ) )
    {
        for(int i = 0; i < pdb.atomCount; i++)
        {
            r[i].x = dcd_xyz.frame.X[i];
            r[i].y = dcd_xyz.frame.Y[i];
            r[i].z = dcd_xyz.frame.Z[i];
            r[i].fi    = dcd_ang.frame.X[i];
            r[i].psi   = dcd_ang.frame.Y[i];
            r[i].theta = dcd_ang.frame.Z[i];
        }
    }
}

void readRestart()
{
    for(int traj = 0; traj < par.Ntr; traj++)
    {
        readXYZ(par.restart_xyzFilename[traj], &xyz);
        for(int i = 0; i < par.Ntot; i++)
        {
            r[i + par.Ntot*traj].x = xyz.atoms[i].x;
            r[i + par.Ntot*traj].y = xyz.atoms[i].y;
            r[i + par.Ntot*traj].z = xyz.atoms[i].z;
        }
        readXYZ(par.restart_angFilename[traj], &xyz);
        for(int i = 0; i < par.Ntot; i++)
        {
            r[i + par.Ntot*traj].fi = xyz.atoms[i].x;
            r[i + par.Ntot*traj].psi = xyz.atoms[i].y;
            r[i + par.Ntot*traj].theta = xyz.atoms[i].z;
        }
    }
    free(xyz.atoms);
}

void writeRestart(long long int step)
{
    xyz.atomCount = par.Ntot;
    xyz.atoms = (XYZAtom*)calloc(xyz.atomCount, sizeof(XYZAtom));
    for(int traj = 0; traj < par.Ntr; traj++)
    {
        for(int i = 0; i < par.Ntot; i++)
        {
            xyz.atoms[i].x = r[i + par.Ntot*traj].x;
            xyz.atoms[i].y = r[i + par.Ntot*traj].y;
            xyz.atoms[i].z = r[i + par.Ntot*traj].z;
            xyz.atoms[i].name = pdb.atoms[i].name[0];
        }
        writeXYZ(par.restart_xyzFilename[traj], &xyz);
        for(int i = 0; i < par.Ntot; i++)
        {
            xyz.atoms[i].x = r[i + par.Ntot*traj].fi;
            xyz.atoms[i].y = r[i + par.Ntot*traj].psi;
            xyz.atoms[i].z = r[i + par.Ntot*traj].theta;
            xyz.atoms[i].name = pdb.atoms[i].name[0];
        }
        writeXYZ(par.restart_angFilename[traj], &xyz);
    }
    free(xyz.atoms);
    FILE *keyf = safe_fopen(par.restartkey, "w");
    if (fprintf(keyf,"%lld", step) == 0)
        DIE("Reading restartkey %s: unable to write firststep", &(par.restartkey));
    fclose(keyf);
}