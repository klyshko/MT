#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "pdbio.h"

/*
 * pdbio.c
 *
 *  Created on: Nov 9, 2008
 *      Author: zhmurov
 */

#define BUF_SIZE 80
//#define DEBUGPDBIO

/*
 * Private methods
 */
void parseAtomLine(PDB* pdbData, char* line, int currentAtom);
void parseSSBondLine(PDB* pdbData, char* line, int currentSSBond);


/*
 * Parses data from PDB (Protein Data Bank) file format into PDB object.
 * Only ATOM and SSBOND entries are considered, some fields are ignored
 * (marked befor corresponding method).
 *
 * For detailed information about PDB file format, please visit
 * http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Parameters:
 * 		filename: name of a pdb file to parse
 * 		pdbData: pointer of an object to save data into
 */

void readPDB(const char* filename, PDB* pdbData){
	printf("Reading %s.\n", filename);
	int ss_count = 0, atoms_count = 0;
	char buffer[BUF_SIZE];
	FILE* file = fopen(filename, "r");
	if ( file != NULL ){
		while(fgets(buffer, BUF_SIZE, file) != NULL){
			if(strncmp(buffer,"SSBOND",6) == 0){
				ss_count++;
			}
			if(strncmp(buffer, "ATOM", 4) == 0){
				atoms_count++;
			}
		}
		printf("Found %d atoms.\n", atoms_count);

		pdbData->atomCount = atoms_count;
		pdbData->ssCount = ss_count;
		pdbData->atoms = (PDBAtom*)malloc(atoms_count*sizeof(PDBAtom));
		pdbData->ssbonds = (PDBSSBond*)malloc(ss_count*sizeof(PDBSSBond));

		int current_atom = 0;
		int current_ss = 0;

		rewind(file);

		while(fgets(buffer, BUF_SIZE, file) != NULL){
			char* pch = strtok(buffer, " ");
			if(strcmp(pch, "SSBOND") == 0){
				parseSSBondLine(pdbData, buffer, current_ss);
				current_ss++;
			}
			if(strcmp(pch, "ATOM") == 0){
				parseAtomLine(pdbData, buffer, current_atom);
				current_atom ++;
			}

		}
	printf("Done reading '%s'.\n", filename);
	fclose(file);
	} else {
		perror(filename);
		exit(0);
	}
}

/*
 * Reads ONLY coordinates from pdb file.
 *
 * Parameters:
 * 		filename: name of a pdb file to parse
 * 		x, y, z: pointers to arrays to save data into
 * 		count: number of atoms for validation
 */
void readCoordinatesFromPDB(const char* filename, double* x, double* y, double* z, int count){
	int atomsCount = 0;
	char buffer[BUF_SIZE];
	FILE* file = fopen(filename, "r");
	if(file != NULL){
		while(fgets(buffer, BUF_SIZE, file) != NULL){
			if(strncmp(buffer, "ATOM", 4) == 0){
				atomsCount++;
			}
		}
		if(atomsCount != count){
			printf("ERROR: Number of atoms in pdb file is wrong: %s\n", filename);
			exit(0);
		}

		int currentAtom = 0;

		rewind(file);

		while(fgets(buffer, BUF_SIZE, file) != NULL){
			char* pch = strtok(buffer, " ");
			if(strcmp(pch, "ATOM") == 0){

				char xstring[9], ystring[9], zstring[9];

				strncpy(xstring, &buffer[30], 8);
				xstring[8] = '\0';

				strncpy(ystring, &buffer[38], 8);
				ystring[8] = '\0';

				strncpy(zstring, &buffer[46], 8);
				zstring[8] = '\0';

				x[currentAtom] = atof(xstring);
				y[currentAtom] = atof(ystring);
				z[currentAtom] = atof(zstring);

				currentAtom ++;
			}

		}
	fclose(file);
	} else {
		perror(filename);
		exit(0);
	}
}

/*
 * Parses single line of 'ATOM' entry from pdb file.
 * ATOM entry format in PDB:
 *
 * COLUMNS      DATA TYPE        FIELD      DEFINITION
 * ------------------------------------------------------
 *  1 -  6      Record name      "ATOM    "
 *  7 - 11      Integer          id		    Atom serial number.
 * 13 - 16      Atom             name       Atom name.
 * 17           Character        altLoc     Alternate location indicator.
 * 18 - 20      Residue name     resName    Residue name.
 * 22           Character        chainID    Chain identifier.
 * 23 - 26      Integer          resSeq     Residue sequence number.
 * 27           AChar            iCode      Code for insertion of residues. (Ignored)
 * 31 - 38      Real(8.3)        x          Orthogonal coordinates for X in Angstroms
 * 39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in Angstroms
 * 47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in Angstroms
 * 55 - 60      Real(6.2)        occupancy  Occupancy.
 * 61 - 66      Real(6.2)        tempFactor Temperature factor.
 * 77 - 78      LString(2)       element    Element symbol, right-justified. (Ignored)
 * 79 - 80      LString(2)       charge     Charge on the atom. (Ignored)
 *
 * For more details visit http://www.wwpdb.org/documentation/format23/sect9.html
 *
 *
 * Method parameters:
 * 		pdbData:	pointer to pdbData object to read into
 * 		line: line from the file to parse
 * 		currentSSBond: indicates location in array of ssbonds
 * 			where new bond should be saved
 */

void parseAtomLine(PDB* pdbData, char* line, int currentAtom){
	char id[6];
	char atomName[5];
	char resName[4], chain, altLoc;
	char resid[5];
	char x[9], y[9], z[9];
	char occupancy[7];
	char beta[7];
	strncpy(id, &line[6], 5);
	id[5] = '\0';
	strncpy(atomName, &line[12], 4);
	atomName[4] = '\0';
	altLoc = line[16];
	strncpy(resName, &line[17], 3);
	resName[3] = '\0';
	strncpy(&chain, &line[21], 1);
	strncpy(resid, &line[22], 4);
	resid[4] = '\0';
	strncpy(x, &line[30], 8);
	x[8] = '\0';
	strncpy(y, &line[38], 8);
	y[8] = '\0';
	strncpy(z, &line[46], 8);
	z[8] = '\0';
	strncpy(occupancy, &line[54], 6);
	occupancy[6] = '\0';
	strncpy(beta, &line[60], 6);
	beta[6] = '\0';
	strcpy(pdbData->atoms[currentAtom].name, strtok(atomName, " "));
	pdbData->atoms[currentAtom].altLoc = altLoc;
	pdbData->atoms[currentAtom].name[4] = 0;
	pdbData->atoms[currentAtom].chain = chain;
	pdbData->atoms[currentAtom].resid = atoi(resid);
	strcpy(pdbData->atoms[currentAtom].resName, resName);
	pdbData->atoms[currentAtom].resName[3] = 0;
	pdbData->atoms[currentAtom].id = atoi(id);
	pdbData->atoms[currentAtom].x = atof(x);
	pdbData->atoms[currentAtom].y = atof(y);
	pdbData->atoms[currentAtom].z = atof(z);
	pdbData->atoms[currentAtom].occupancy = atof(occupancy);
	pdbData->atoms[currentAtom].beta = atof(beta);
#ifdef DEBUGPDBIO
	printAtom(pdbData->atoms[currentAtom]);
#endif
}

/*
 * Parses single line of 'SSBOND' entry from pdb file.
 * SSBOND entry format in PDB:
 *
 * COLUMNS        DATA TYPE       FIELD         DEFINITION
 * -------------------------------------------------------------------
 *  1 -  6        Record name     "SSBOND"
 *  8 - 10        Integer         serNum       Serial number.
 * 12 - 14        LString(3)      "CYS"        Residue name.
 * 16             Character       chainID1     Chain identifier.
 * 18 - 21        Integer         seqNum1      Residue sequence number.
 * 22             AChar           icode1       Insertion code. (Ignored)
 * 26 - 28        LString(3)      "CYS"        Residue name.
 * 30             Character       chainID2     Chain identifier.
 * 32 - 35        Integer         seqNum2      Residue sequence number.
 * 36             AChar           icode2       Insertion code. (Ignored)
 * 60 - 65        SymOP           sym1         Symmetry oper for 1st resid (Ignored)
 * 67 - 72        SymOP           sym2         Symmetry oper for 2nd resid (Ignored)
 *
 * For more details visit http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Method parameters:
 * 		pdbData:	pointer to pdbData object to read into
 * 		line: line from the file to parse
 * 		currentSSBond: indicates location in array of ssbonds
 * 			where new bond should be saved
 */
void parseSSBondLine(PDB* pdbData, char* line, int currentSSBond){
	char chain;
	char resid[4];
	strncpy(&chain, &line[15], 1);
	strncpy(resid, &line[17], 4);
	pdbData->ssbonds[currentSSBond].chain1 = chain;
	pdbData->ssbonds[currentSSBond].resid1 = atoi(resid);
	strncpy(&chain, &line[29], 1);
	strncpy(resid, &line[31], 4);
	pdbData->ssbonds[currentSSBond].chain2 = chain;
	pdbData->ssbonds[currentSSBond].resid2 = atoi(resid);
#ifdef DEBUG
	printf("SS #%d: %c%d - %c%d\n", current_ss, pdbData->ssbonds[current_ss].chain1, pdbData->ssbonds[current_ss].resid1,
			pdbData->ssbonds[current_ss].chain2, pdbData->ssbonds[current_ss].resid2);
#endif
}

/*
 * Saves data from PDB object into a PDB file format
 *
 * For detailed information about PDB file format, please visit
 * http://www.wwpdb.org/documentation/format23/sect6.html
 *
 * Parameters:
 * 		filename: name of a file to save into (will be overwritten/created)
 * 		pdbData: pointer of an object to get data from
 */
void writePDB(const char* filename, PDB* pdbData){
	writePDB(filename, pdbData, 0);
}


void writePDB(const char* filename, PDB* pdbData, int printConnections){
	printf("Saving PDB '%s'...\n", filename);
	FILE* file = fopen(filename, "w");
	int i, j;
	for(i = 0; i < pdbData->ssCount; i++){
		fprintf(file, "SSBOND %3d CYS %c %4d    CYS %c %4d\n",
								i + 1,
								pdbData->ssbonds[i].chain1,
								pdbData->ssbonds[i].resid1,
								pdbData->ssbonds[i].chain2,
								pdbData->ssbonds[i].resid2);
	}
	if(pdbData->atomCount < 100000){
		for(i = 0; i < pdbData->atomCount; i++){
			fprintf(file, "ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
									i + 1,
									pdbData->atoms[i].name,
									pdbData->atoms[i].altLoc,
									pdbData->atoms[i].resName,
									pdbData->atoms[i].chain,
									pdbData->atoms[i].resid,
									pdbData->atoms[i].x,
									pdbData->atoms[i].y,
									pdbData->atoms[i].z,
									pdbData->atoms[i].occupancy,
									pdbData->atoms[i].beta);
		}
	} else {
		for(i = 0; i < pdbData->atomCount; i++){
			fprintf(file, "ATOM  %5x %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
									i + 1,
									pdbData->atoms[i].name,
									pdbData->atoms[i].altLoc,
									pdbData->atoms[i].resName,
									pdbData->atoms[i].chain,
									pdbData->atoms[i].resid,
									pdbData->atoms[i].x,
									pdbData->atoms[i].y,
									pdbData->atoms[i].z,
									pdbData->atoms[i].occupancy,
									pdbData->atoms[i].beta);
		}
	}

	if(printConnections){
		for(i = 0; i < pdbData->atomCount; i++){
			fprintf(file, "CONECT%5d", i+1);
			for(j = 0; j < pdbData->connections.connectCount[i]; j++){
				fprintf(file, "%5d", pdbData->connections.connectMap[j*pdbData->atomCount + i]+1);
			}
			fprintf(file, "\n");
		}
	}
	fprintf(file, "END");
	fclose(file);
	printf("Done saving PDB.\n");
}

/*
 * Prints PDBAtom object in a PDB format.
 * Parameters:
 * 		atomData: PDBAtom object to print
 */
void printAtom(PDBAtom atomData){
	printf("ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
			atomData.id,
			atomData.name,
			atomData.altLoc,
			atomData.resName,
			atomData.chain,
			atomData.resid,
			atomData.x,
			atomData.y,
			atomData.z,
			atomData.occupancy,
			atomData.beta);
}

/*
 * Prints PDBAtom object in a PDB format into a file.
 * Parameters:
 * 		atomData: PDBAtom object to print
 */
void printAtomToFile(FILE* file, PDBAtom atomData){
	fprintf(file, "ATOM  %5d %-4s%c%3s %c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %c\n",
			atomData.id,
			atomData.name,
			atomData.altLoc,
			atomData.resName,
			' ',
			atomData.resid,
			atomData.x,
			atomData.y,
			atomData.z,
			atomData.occupancy,
			atomData.beta,
			atomData.chain);
}
