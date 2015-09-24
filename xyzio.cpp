/*
 * xyzio.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: zhmurov
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "wrapper.h"
#include "xyzio.h"

#define BUF_SIZE 256

void readXYZ(const char* filename, XYZ* xyzData){
    printf("Reading %s.\n", filename);
    char buffer[BUF_SIZE];
    FILE* file = safe_fopen(filename, "r");
    if(file != NULL){
        safe_fgets(buffer, BUF_SIZE, file);
        xyzData->atomCount = atoi(buffer);
        safe_fgets(buffer, BUF_SIZE, file);
        xyzData->atoms = (XYZAtom*)calloc(xyzData->atomCount, sizeof(XYZAtom));
        int i;
        char* pch;
        for(i = 0; i < xyzData->atomCount; i++){
            safe_fgets(buffer, BUF_SIZE, file);
            pch = strtok(buffer, " \t\r\n");
            xyzData->atoms[i].name = pch[0];
            pch = strtok(NULL, " \t\r\n");
            xyzData->atoms[i].x = atof(pch);
            pch = strtok(NULL, " \t\r\n");
            xyzData->atoms[i].y = atof(pch);
            pch = strtok(NULL, " \t\r\n");
            xyzData->atoms[i].z = atof(pch);
        }
        printf("Done reading '%s'.\n", filename);
        fclose(file);
    } else {
        DIE("Error opening '%s': %s", filename, strerror(errno));
    }
}

void readCoordinatesFromXYZ(const char* filename, double* x, double* y, double* z, int count){
    char buffer[BUF_SIZE];
    FILE* file = safe_fopen(filename, "r");
    if(file != NULL){
        safe_fgets(buffer, BUF_SIZE, file);
        int atomsCount = atoi(buffer);
        if(atomsCount != count){
            DIE("ERROR: Number of atoms in xyz file is wrong: %s", filename);
        }
        safe_fgets(buffer, BUF_SIZE, file);
        int i;
        char* pch;
        for(i = 0; i < atomsCount; i++){
            safe_fgets(buffer, BUF_SIZE, file);
            pch = strtok(buffer, " \t\r\n");
            pch = strtok(NULL, " \t\r\n");
            x[i] = atof(pch);
            pch = strtok(NULL, " \t\r\n");
            y[i] = atof(pch);
            pch = strtok(NULL, " \t\r\n");
            z[i] = atof(pch);
        }
    fclose(file);
    } else {
        DIE("Error opening '%s': %s", filename, strerror(errno));
    }
}

void writeXYZ(const char* filename, XYZ* xyzData){
    printf("Writing %s.\n", filename);
    FILE* file = safe_fopen(filename, "w");
    if(file != NULL){
        int i;
        fprintf(file, "%d\n", xyzData->atomCount);
        fprintf(file, "Created by 'xyzio.cpp'\n");
        for(i = 0; i < xyzData->atomCount; i++){
            fprintf(file, "%-*c%*f%*f%*f\n",
                    XYZ_FIELD_SIZE, xyzData->atoms[i].name,
                    XYZ_FIELD_SIZE, xyzData->atoms[i].x,
                    XYZ_FIELD_SIZE, xyzData->atoms[i].y,
                    XYZ_FIELD_SIZE, xyzData->atoms[i].z);
        }
        printf("Done writing '%s'.\n", filename);
        fclose(file);
    } else {
        DIE("Error opening '%s': %s", filename, strerror(errno));
    }
}
