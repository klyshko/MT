#include "pdbio.h"
#include "dcdio.h"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <vector>

PDB pdbdata;
DCD dcdxyz, dcdang, dcdout;

void allocdcd(DCD* dcd, int n){
		dcd->frame.X = (float*)calloc(n, sizeof(float));
		dcd->frame.Y = (float*)calloc(n, sizeof(float));
		dcd->frame.Z = (float*)calloc(n, sizeof(float));
}

int main(int argc, char** argv){
		if(argc < 3){
				printf("Usage: %s infile.pdb infile.xyz.dcd infile.ang.dcd outfile.dcd\n", argv[0]);
				exit(1);
		}
		readPDB(argv[1], &pdbdata);
		dcdOpenRead(&dcdxyz, argv[2]);
		dcdOpenRead(&dcdang, argv[3]);
		dcdOpenWrite(&dcdout, argv[4]);
		dcdReadHeader(&dcdxyz);
		dcdReadHeader(&dcdang);
        createDCD(&dcdout, pdbdata.atomCount, dcdxyz.header.nfile, 0, 1.0, dcdxyz.header.nsavc, 0, 0, 0, 0);
        dcdWriteHeader(dcdout);
        allocdcd(&dcdxyz, pdbdata.atomCount);
        allocdcd(&dcdang, pdbdata.atomCount);
        allocdcd(&dcdout, pdbdata.atomCount);
		long int frames = 0;
		do{
			dcdReadFrame(&dcdxyz); dcdReadFrame(&dcdang);
            frames ++;
			for(int i=0; i < pdbdata.atomCount; i++){
                float x = dcdxyz.frame.X[i];
                float y = dcdxyz.frame.Y[i];
                float z = dcdxyz.frame.Z[i];
                float fi    = dcdang.frame.X[i];
                float psi   = dcdang.frame.Y[i];
                float theta = dcdang.frame.Z[i];
                dcdout.frame.X[i] = sqrtf(x*x + y*y);
                dcdout.frame.Y[i] = z;
                dcdout.frame.Z[i] = theta;
            }
            dcdWriteFrame(dcdout);
		}while(!feof(dcdxyz.file) && !feof(dcdang.file));
        printf("%ld frames processed!\nHave a nice day!\n", frames);
        // What, you really expected me to deallocate memory? 
		return 0;
}
