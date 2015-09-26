/*
 * main.cpp
 *
 *  Created on: 11.04.2012
 *      Author: zhmurov
 */

#include "mt.h"
#include "parameters.h"
#include "configreader.h"
#include "timer.h"
#include "wrapper.h"
#include "num_test.h"
#include "output.h"
#ifdef USE_MPI
    #include <mpi.h>
#endif 
//void createMT();
void read_PDB(const char* filename_xyz, const char* filename_ang);
void saveCoordPDB(const char* pdbfilename_xyz, const char* pdbfilename_ang);
void ReadFromDCD(Parameters par, Topology top, char* dcdfilename_xyz, char* dcdfilename_ang);
void saveCoordDCD();
void update(long long int step);
void UpdateLJPairs();
void UpdatePairs();
void AssemblyInit();
void writeRestart(long long int step);
void readRestart();
void OutputEnergies(int index, int traj);
void OutputAllEnergies();

extern void init(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies);
extern void compute(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies);
#ifdef TEST_OPT
extern void num_test(Coord* r, Coord* f, Parameters &par, Topology &top);
#endif

#ifdef TEST_ENERGY
extern double U(Coord* r, Parameters &par, Topology &top);
#endif

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

void initParameters(int argc, char* argv[]);

int main(int argc, char *argv[]){
    printf("Start parsing...\n");
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

    //init(r, f, par, top, energies);

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

//#ifndef  TEST_OPT
    initTimer();

#if defined(ASSEMBLY)
    AssemblyInit();
    //UpdatePairs();
#endif

#ifdef OUTPUT_EN
    energies = (Energies*)malloc(par.Ntot * par.Ntr * sizeof(Energies));
#endif

    compute(r, f, par, top, energies);
//#else
  //  num_test(r, f, par, top);
//#endif

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
#ifdef OUTPUT_EN
    free(energies);
#endif
    return 0;
}

void update(long long int step){
    printf("Saving coordinates at step %lld\n", step);
    printTime(step);
    printEstimatedTimeleft((real)step/(real)par.steps);
    saveCoordDCD();
 /*   
#ifdef TEST_ENERGY
    printf("%f\n", U(r, par, top));
#endif
*/
#ifdef OUTPUT_FORCE
    Coord f_sum;
    f_sum.x = 0; f_sum.y = 0; f_sum.z = 0;
    for(int i = 0; i < par.Ntot; i++)
    {
        f_sum.x += f[i].x;
        f_sum.y += f[i].y;
        f_sum.z += f[i].z;
    }
    printf("SF %.16f\n", sqrt(f_sum.x * f_sum.x + f_sum.y * f_sum.y + f_sum.z * f_sum.z));
#endif

#ifdef OUTPUT_EN
    OutputAllEnergies();
#endif

}

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
#ifndef TEST_OPT
    par.steps = getLongIntegerParameter(PARAMETER_NUMSTEPS,-1);
#else
    par.steps = 1;
#endif

#ifdef T3LS
    par.steps = 1;
#endif
    par.stride = getLongIntegerParameter(PARAMETER_STRIDE,-1);
#ifdef LJ_on
    par.ljpairscutoff = getFloatParameter(PARAMETER_LJPAIRSCUTOFF);
    par.ljpairsupdatefreq = getIntegerParameter(PARAMETER_LJPAIRSUPDATEFREQ);
#endif
    getMaskedParameter(par.coordFilename_ang, PARAMETER_COORD_FILE_ANG);
    getMaskedParameter(par.coordFilename_xyz, PARAMETER_COORD_FILE_XYZ);
    getMaskedParameter(par.ffFilename, PARAMETER_FORCEFIELD_FILE);
    getMaskedParameter(par.condFilename, PARAMETER_CONDITIONS_FILE);
    par.fix = getIntegerParameter(PARAMETER_FIX,1);
#ifndef TEST_OPT
    par.Ntr = getIntegerParameter(PARAMETER_RUNNUM,1);
#else
    par.Ntr = 1;
#endif
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
    read_PDB(par.coordFilename_xyz, par.coordFilename_ang);

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
    par.C = getFloatParameter(PARAMETER_HARMONIC_C);
    par.A_long = getFloatParameter(A_LONG);
    par.A_lat  = getFloatParameter(A_LAT);
    par.D_long  = getFloatParameter(D_LONG);
    par.D_lat   = getFloatParameter(D_LAT);
#else
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

#if defined(BOUNDARIES)
    par.xmax_bound = getFloatParameter(XMAX_BOUND);
    par.xmin_bound = getFloatParameter(XMIN_BOUND);
    par.ymax_bound = getFloatParameter(YMAX_BOUND);
    par.ymin_bound = getFloatParameter(YMIN_BOUND);
    par.zmax_bound = getFloatParameter(ZMAX_BOUND);
    par.zmin_bound = getFloatParameter(ZMIN_BOUND);
    par.ks_bound = getFloatParameter(KS_BOUND);
    
#endif

    parseParametersFile(par.condFilename, argc, argv);
    par.Temp = getFloatParameter(PARAMETER_TEMPERATURE);
    par.gammaR = getFloatParameter(PARAMETER_GAMMA_R);
    par.gammaTheta = getFloatParameter(PARAMETER_GAMMA_THETA);
    par.varR = sqrtf(2.0f*KB*par.Temp*par.dt/par.gammaR);
    par.varTheta = sqrtf(2.0f*KB*par.Temp*par.dt/par.gammaTheta);
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



void read_PDB(const char* filename_xyz, const char* filename_ang){
    int i, j;
    readPDB(filename_xyz, &pdb);
    printf("Building topology....\n");
    par.Ntot = pdb.atomCount;
    r = (Coord*)calloc(par.Ntot*par.Ntr, sizeof(Coord));
    f = (Coord*)calloc(par.Ntot*par.Ntr, sizeof(Coord));

    top.mon_type = (int*)calloc(par.Ntot, sizeof(int));
    for(i = 0; i < par.Ntot; i++){
        if(pdb.atoms[i].name == "CA")
            top.mon_type[i] = 0;
        else if(pdb.atoms[i].name == "CB")
            top.mon_type[i] = 1;

    }

#ifdef LJ_on
    top.LJCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
#endif
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
    top.harmonicCount = (int*)calloc(par.Ntot, sizeof(int));
    top.longitudinalCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
    top.lateralCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
    
    // Covalent bonds map
    for(i = 0; i < par.Ntot; i++) top.harmonicCount[i] = 0;
    for(i = 0; i < par.Ntot; i++){
        for(j = 0; j < par.Ntot; j++){
            if(pdb.atoms[i].resid == pdb.atoms[j].resid &&
                    pdb.atoms[i].chain == pdb.atoms[j].chain
                    && (i!=j) ){
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
                    (i!=j)){
                if(pdb.atoms[i].id > pdb.atoms[j].id){
                    top.harmonic[top.maxHarmonicPerMonomer*i + top.harmonicCount[i]] = j; // << BUG: indexation
                }
                else{
                    top.harmonic[top.maxHarmonicPerMonomer*i + top.harmonicCount[i]] = -j; /* BUG: same here */
                }
                top.harmonicCount[i]++;
            }
        }
    }

    // Longitudal exponential
    for(i = 0; i < par.Ntot * par.Ntr; i++) top.longitudinalCount[i] = 0;

    for (int traj = 0; traj < par.Ntr; traj++){
        for(i = 0; i < par.Ntot; i++){
            for(j = 0; j < par.Ntot; j++){
                if( (pdb.atoms[i].resid == (pdb.atoms[j].resid + 1)) &&
                    (pdb.atoms[i].chain == pdb.atoms[j].chain) &&
                    (pdb.atoms[i].id == (pdb.atoms[j].id + 1))){
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
            for(j = 0; j < par.Ntot; j++){
                if(abs(pdb.atoms[i].resid - pdb.atoms[j].resid) == 1 &&
                    pdb.atoms[i].chain == pdb.atoms[j].chain &&
                    abs(pdb.atoms[i].id -  pdb.atoms[j].id) == 1){
                        if(pdb.atoms[i].id > pdb.atoms[j].id){
                            top.longitudinal[top.maxLongitudinalPerMonomer*i + top.longitudinalCount[i]] = j; // << BUG: indexation
                        }
                        else{
                            top.longitudinal[top.maxLongitudinalPerMonomer*i + top.longitudinalCount[i]] = -j; /* BUG: same here */
                        }
                    top.longitudinalCount[i]++;
                }
            }
        }
        
    }
        for(i = 0; i < par.Ntot; i++){
            for(j = 0; j < par.Ntot; j++){
                if(abs(pdb.atoms[i].resid - pdb.atoms[j].resid) == 1 &&
                    pdb.atoms[i].chain == pdb.atoms[j].chain &&
                    abs(pdb.atoms[i].id -  pdb.atoms[j].id) == 1){
                        if(pdb.atoms[i].id > pdb.atoms[j].id){
                            top.longitudinal[top.maxLongitudinalPerMonomer*i + top.longitudinalCount[i]] = j; // << BUG: indexation
                        }
                        else{
                            top.longitudinal[top.maxLongitudinalPerMonomer*i + top.longitudinalCount[i]] = -j; /* BUG: same here */
                        }
                    top.longitudinalCount[i]++;
                }
            }
        }
    
    // Lateral exponential
    for(i = 0; i < par.Ntot; i++) top.lateralCount[i] = 0;
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
    for(i = 0; i < par.Ntot; i++){
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
            for(int ind=0; ind < 2; ind++)
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
        /*      if(dr < 10){
                    char nitro = 'N';
                    printf("%-*c%*f%*f%*f\n",16,nitro,16,
                    xi - yp2 * cos_fii * sin_psii + zp2 * sin_fii * sin_psii + cos_psii * (xp2 * cos_thetai + zp2 * cos_fii * sin_thetai + yp2 * sin_fii * sin_thetai),
                    16,
                    yi - zp2 * cos_psii * sin_fii + xp2 * cos_thetai * sin_psii + yp2 * sin_fii * sin_psii * sin_thetai + cos_fii * (yp2 * cos_psii + zp2 * sin_psii * sin_thetai),
                    16,
                    zi + zp2 * cos_fii * cos_thetai + yp2 * cos_thetai * sin_fii - xp2 * sin_thetai
                    );
        //          printf("%f %f %f\n", xp2,yp2,zp2);

                }*/
                if(dr < 1*r_mon){
                    top.lateralCount[i]++;
                    if(top.maxLateralPerMonomer < top.lateralCount[i])
                            top.maxLateralPerMonomer = top.lateralCount[i];
                }
            }
        }   
    }
    top.lateral = (int*)calloc(top.maxLateralPerMonomer*par.Ntot, sizeof(int));
    for(i = 0; i < par.Ntot; i++) top.lateralCount[i] = 0;
    for(i = 0; i < par.Ntot; i++){
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
            for(int ind=0; ind < 2; ind++)
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
                        top.lateral[top.maxLateralPerMonomer*i + top.lateralCount[i]] = j;
                    }else{
                        top.lateral[top.maxLateralPerMonomer*i + top.lateralCount[i]] = -j;
                    }
                    top.lateralCount[i]++;
                }
            }
        }
    }
    //todo initialising fixed atoms list
    top.fixed = (bool*)calloc(par.Ntot, sizeof(bool));
    for(i=0; i < par.Ntot; i++){
            if( pdb.atoms[i].resid <= par.fix )
                    top.fixed[i] = true;
            else
                    top.fixed[i] = false;
    }
#ifdef LJ_on
    //UpdateLJPairs();
#endif
    printf("done building topology without LJ.\n");
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

void OutputAllEnergies(){

    Energies fullEnergy;
    fullEnergy.U_harm = 0;
    fullEnergy.U_long = 0;
    fullEnergy.U_lat = 0;
    fullEnergy.U_psi = 0;
    fullEnergy.U_fi = 0;
    fullEnergy.U_teta = 0;
    fullEnergy.U_lj = 0;

    for (int i = 0; i < par.Ntot * par.Ntr; i++){
        fullEnergy.U_harm += energies[i].U_harm;
        fullEnergy.U_long += energies[i].U_long;
        fullEnergy.U_lat += energies[i].U_lat;
        fullEnergy.U_psi += energies[i].U_psi;
        fullEnergy.U_fi += energies[i].U_fi;
        fullEnergy.U_teta += energies[i].U_teta;
        fullEnergy.U_lj += energies[i].U_lj;
    }
    printf("ENERGIES harm = %f long = %f lat = %f psi = %f fi = %f theta = %f lj = %f\n", fullEnergy.U_harm, fullEnergy.U_long, fullEnergy.U_lat, fullEnergy.U_psi, fullEnergy.U_fi, fullEnergy.U_teta, fullEnergy.U_lj);
}

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
       // U(r, par, top);
    }
}

void AssemblyInit()
{
    free(top.lateral);
    free(top.longitudinal);

    top.lateral = (int*)calloc(2*par.Ntot*par.Ntr, sizeof(int));
    top.longitudinal = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
    for(int i = 0; i < par.Ntot*par.Ntr; i++)
    {
        top.longitudinalCount[i] = 0;
        top.lateralCount[i] = 2;

    }
    top.maxLongitudinalPerMonomer = 1;
    top.maxLateralPerMonomer = 2;
}