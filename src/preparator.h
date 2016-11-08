#include "parameters.h"
#include "configreader.h"
#include "globals.h"
#include "wrapper.h"

void initParameters(int argc, char* argv[]);
void read_PDB(const char* filename_xyz, const char* filename_ang);
void saveCoordPDB(const char* pdbfilename_xyz, const char* pdbfilename_ang);
void ReadFromDCD(Parameters par, Topology top, char* dcdfilename_xyz, char* dcdfilename_ang);
void saveCoordDCD();
void writeRestart(long long int step);
void readRestart();
void AssemblyInit();

void appendCoordPDB();