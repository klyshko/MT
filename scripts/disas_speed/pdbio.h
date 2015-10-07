/*
 * pdbio.h
 *
 *  Created on: Jul 5, 2009
 *      Author: zhmurov
 */

#pragma once

#include <stdio.h>

/*
 * Structures
 */

typedef struct {

  int id;
  char   name[5], chain, resName[4], altLoc;
  int    resid;
  double x, y, z;

  double occupancy;
  double beta;

} PDBAtom;

typedef struct {

	int resid1;
	char chain1;

	int resid2;
	char chain2;

} PDBSSBond;

typedef struct {
	int* connectCount;
	int* connectMap;
} PDBConnect;

typedef struct {
	int atomCount;
	PDBAtom* atoms;
	int ssCount;
	PDBSSBond* ssbonds;
	PDBConnect connections;
} PDB;


/*
 * Public methods
 */
void readPDB(const char* filename, PDB* pdbData);
void readCoordinatesFromPDB(const char* filename, double* x, double* y, double* z, int count);
void writePDB(const char* filename, PDB* pdbData);
void writePDB(const char* filename, PDB* pdbData, int printConnections);
void printAtom(PDBAtom atomData);
void printAtomToFile(FILE* file, PDBAtom atomData);
