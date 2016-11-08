/*
 *  dcdio.c
 *
 * Code to work with dcd files.
 *
 *  Created on: Oct 29, 2008
 *      Author: zhmurov
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "dcdio.h"

void createDCD(DCD* dcd, int atomCount, int frameCount, int firstFrame, float timestep, int dcdFreq, int hazUC, int UCX, int UCY, int UCZ){
	dcd->header.N = atomCount;
	dcd->header.delta = timestep;
	dcd->header.nfile = frameCount;
	dcd->header.npriv = firstFrame;
	dcd->header.nsavc = dcdFreq;
	dcd->header.nstep = frameCount*dcdFreq;
	sprintf(dcd->header.remark1, "REMARKS CREATED BY dcdio.c");
	time_t rawtime;
	struct tm * timeinfo;
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	sprintf(dcd->header.remark2, "REMARKS DATE: %s", asctime(timeinfo));
	dcd->hasUC = hazUC;
	dcd->uc.a = UCX;
	dcd->uc.b = UCY;
	dcd->uc.c = UCZ;
	dcd->uc.alpha = 0;
	dcd->uc.beta = 0;
	dcd->uc.gamma = 0;
	dcd->frame.X = (float*)calloc(atomCount, sizeof(float));
	dcd->frame.Y = (float*)calloc(atomCount, sizeof(float));
	dcd->frame.Z = (float*)calloc(atomCount, sizeof(float));
}

/*
 * Extend the string to be 80 bytes length (for remarks in header)
 * Taken from dcdlib.C (NAMD source)
 */
void pad(char *s, int len);

/*
 * Open DCD file to write data into.
 * The FILE* pointer is returned and must be used for further calls
 */
void dcdOpenWrite(DCD* dcd, char *dcd_filename){
	dcd->file = fopen(dcd_filename, "w");
}

void dcdOpenAppend(DCD* dcd, char *dcd_filename){
	dcd->file = fopen(dcd_filename, "a");
}

/*
 * Write header to dcd file.
 * Inputs are:
 *  dcd_file - FILE* object from dcd_open_write call
 *  dcd_filename - dcd filename to write into remarks
 *  N - Number of atoms/beads in a system
 *	NFILE - Number of frames (will not be updated in current implementation)
 *	NPRIV - Starting timestep of DCD file - NOT ZERO
 *	NSAVC - Timesteps between DCD saves
 *	NSTEP - Number of timesteps
 *	DELTA - length of a timestep
 *
 *  DCD header format is the following:
 *
 *  Position	Length			data
 *  (bytes/4)	(bytes)
 *  ------------------------------------
 *  1			4				'84'
 *  2			4				'CORD'
 *  3			4			  	NFILE
 *  4			4				NPRIV
 *  5			4				NSAVC
 *  6			4				NPRIV-NSAVC or NSTEP
 *  7-11		5*4				Zeros
 *  12			4				DELTA
 *  13			4				With/without unit cell
 *  14-21		8				Unit cell description (or zeros)
 *  22			4				'24'
 *  23			4				'84'
 *  24			4				'164'
 *  25			4				'2'
 *  26-45		80				Remarks #1 ("REMARK CREATED BY NAMD FILENAME='...'")
 *  46-65		80				Remarks #2 ("REMARK" + date + username)
 *  66			4				'164'
 *  67			4				'4'
 *  68			4				N
 *  69			4				'4'
 *
 */
void dcdWriteHeader( DCD dcd){
	FILE* dcd_file = dcd.file;
	int iout;
	float fout;
	char cout[5];
	iout = 84;
	fwrite(&iout, 4, 1, dcd_file);
	sprintf(cout, "CORD");
	fwrite(&cout, 4, 1, dcd_file);
	iout = dcd.header.nfile;
	fwrite(&iout, 4, 1, dcd_file);
	iout = dcd.header.npriv;
	fwrite(&iout, 4, 1, dcd_file);
	iout = dcd.header.nsavc;
	fwrite(&iout, 4, 1, dcd_file);
	iout = dcd.header.npriv-dcd.header.nsavc;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 0;
	for(int i = 0; i < 5; i++){
		fwrite(&iout, 4, 1, dcd_file);
	}
	fout = dcd.header.delta;
	fwrite(&fout, 4, 1, dcd_file);
	fwrite(&dcd.hasUC, 4, 1, dcd_file);
	for(int i = 0; i < 8; i++){
		fwrite(&iout, 4, 1, dcd_file);
	}
	iout = 24;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 84;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 164;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 2;
	fwrite(&iout, 4, 1, dcd_file);
	char title[81];
	strncpy(title,dcd.header.remark1,sizeof(title));
//	sprintf(title, "%s", dcd.header.remark1);
	pad(title, 80);
	fwrite(title, 80, 1, dcd_file);
	strncpy(title,dcd.header.remark2,sizeof(title));
//	sprintf(title, "%s", dcd.header.remark2);
	pad(title, 80);
	fwrite(&title, 80, 1, dcd_file);
	iout = 164;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 4;
	fwrite(&iout, 4, 1, dcd_file);
	iout = dcd.header.N;
	fwrite(&iout, 4, 1, dcd_file);
	iout = 4;
	fwrite(&iout, 4, 1, dcd_file);
}

/*
 * Writes frame into dcd
 * Input
 *  dcd_file - file to write in
 *  N - number of atoms
 *  X, Y, Z - pointers to arrays with coordinates
 *
 *  Writing format is following
 *
 *  Length			Data
 *  (bytes)
 *  ---------------------------(If unit cell is defined - not implemented in this code)
 *  4				'48'
 *  12*4			Unit cell data
 *  4				'48'
 *  ------------------------(If unit cell is defined - not implemented in this code)
 *  4				N*4
 *  N*4				X
 *  4				N*4
 *  4				N*4
 *  N*4				Y
 *  4				N*4
 *  4				N*4
 *  N*4				Z
 *  4				N*4
 *
 */
void dcdWriteFrame(DCD dcd){
	FILE* dcd_file = dcd.file;
	int iout;
	if(dcd.hasUC){
		iout = 48;
		fwrite(&iout, 4, 1, dcd_file);
		fwrite(&dcd.uc.a, 8, 1, dcd_file);
		fwrite(&dcd.uc.alpha, 8, 1, dcd_file);
		fwrite(&dcd.uc.b, 8, 1, dcd_file);
		fwrite(&dcd.uc.beta, 8, 1, dcd_file);
		fwrite(&dcd.uc.gamma, 8, 1, dcd_file);
		fwrite(&dcd.uc.c, 8, 1, dcd_file);
		fwrite(&iout, 4, 1, dcd_file);
	}
	iout = dcd.header.N*4;
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(dcd.frame.X, 4*dcd.header.N, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(dcd.frame.Y, 4*dcd.header.N, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
	fwrite(dcd.frame.Z, 4*dcd.header.N, 1, dcd_file);
	fwrite(&iout, 4, 1, dcd_file);
}

/*
 * Close DCD after writing
 */
void dcdClose(DCD dcd){
	fclose(dcd.file);
}


void dcdOpenRead(DCD* dcd, char *dcd_filename){
	dcd->file = fopen(dcd_filename, "r");
}

void dcdReadHeader(DCD* dcd){
	FILE* dcd_file = dcd->file;
	int iin;
	char cin[5];
	fread(&iin, 4, 1, dcd_file);
	if(iin != 84){
		printf("Error! Wrong DCD file: DCD supposed to have '84' in first 4 bytes, but it hasn't.\n");
		exit(0);
	}
	fread(&cin, 4, 1, dcd_file);
	char cord_char[5];
	sprintf(cord_char, "CORD");
	if(strncmp(cin, cord_char, 4)){
		printf("Error! Wrong DCD file: no 'CORD' sequence at the beginning of the file. Found: %s.\n", cin);
		exit(0);
	}

	fread(&dcd->header.nfile, 4, 1, dcd_file);
	fread(&dcd->header.npriv, 4, 1, dcd_file);
	fread(&dcd->header.nsavc, 4, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	for(int i = 0; i < 5; i++){
		fread(&iin, 4, 1, dcd_file);
	}
	fread(&dcd->header.delta, 4, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	dcd->hasUC = iin;
	for(int i = 0; i < 8; i++){
		fread(&iin, 4, 1, dcd_file);
		//printf("%d: %d\n", i+1, iin);
	}
	//24
	fread(&iin, 4, 1, dcd_file);
	//84;
	fread(&iin, 4, 1, dcd_file);
	//164;
	fread(&iin, 4, 1, dcd_file);
	//2;
	fread(&iin, 4, 1, dcd_file);
	fread(&dcd->header.remark1, 80, 1, dcd_file);
	//printf("Title1: %s\n", title);
	fread(&dcd->header.remark2, 80, 1, dcd_file);
	//printf("Title2: %s\n", title);
	//164
	fread(&iin, 4, 1, dcd_file);
	//4
	fread(&iin, 4, 1, dcd_file);
	//N
	fread(&dcd->header.N, 4, 1, dcd_file);
	//4
	fread(&iin, 4, 1, dcd_file);
}

int dcdReadFrame(DCD* dcd){
	FILE* dcd_file = dcd->file;
	int iin;
	double din;
	fread(&iin, 4, 1, dcd_file);
	if(dcd->hasUC){
		fread(&dcd->uc.a, 4, 1, dcd_file);
		fread(&dcd->uc.alpha, 4, 1, dcd_file);
		fread(&dcd->uc.b, 4, 1, dcd_file);
		fread(&dcd->uc.beta, 4, 1, dcd_file);
		fread(&dcd->uc.gamma, 4, 1, dcd_file);
		fread(&dcd->uc.c, 4, 1, dcd_file);
		/*for(int i = 0; i < 6; i++){
			fread(&din, 8, 1, dcd_file);
			printf("%d: %f\n", i+1, din);
		}*/
		fread(&iin, 4, 1, dcd_file);
		fread(&iin, 4, 1, dcd_file);
	}
	//printf("%d\n", iin);
	fread(dcd->frame.X, 4*dcd->header.N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(dcd->frame.Y, 4*dcd->header.N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	fread(dcd->frame.Z, 4*dcd->header.N, 1, dcd_file);
	fread(&iin, 4, 1, dcd_file);
	//printf("%d\n", iin);
	if(feof(dcd_file) == 0){
		return 0;
	} else {
		return -1;
	}
}

void pad(char *s, int len){
	int curlen;
	int i;
	curlen = strlen(s);
	if (curlen > len){
		s[len-1] = '\0';
		return;
	}
	for (i = curlen; i < len; i++){
		s[i] = '\0';
	}
	s[i] = '\0';
}
