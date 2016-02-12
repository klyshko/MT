#include "pdbio.h"
#include "dcdio.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>

PDB pdbdata;
DCD dcd_xyz, dcd_ang, dcd_xyz_old, dcd_ang_old;
int stride;


using namespace std;

int max(int a, int b)
{
	return (a>b)?a:b;
}
int min(int a, int b)
{
	return (a<b)?a:b;
}

int main(int argc, char** argv){
	if(argc < 5){
		printf("Usage: %s infile.pdb infile_xyz.dcd infile_ang.dcd stride\n", argv[0]);
		exit(1);
	}
    stride = atoi(argv[4]);
	readPDB(argv[1], &pdbdata);
	dcdOpenRead(&dcd_xyz, argv[2]);
	dcdOpenRead(&dcd_ang, argv[3]);
	printf("%s\n%s\n", argv[2], argv[3]);
	dcdReadHeader(&dcd_xyz);
	dcdReadHeader(&dcd_ang);
	long int frames = 0;
    double t_xyz = 0.0, t_rot = 0.0, t_x = 0.0, t_y = 0.0, t_z = 0.0, tr_x = 0.0, tr_y = 0.0, tr_z = 0.0, *tmp;
    double da, db, dg, dx, dy, dz;
	dcd_xyz.frame.X = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_xyz.frame.Y = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_xyz.frame.Z = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_ang.frame.X = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_ang.frame.Y = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_ang.frame.Z = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_xyz_old.frame.X = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_xyz_old.frame.Y = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_xyz_old.frame.Z = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_ang_old.frame.X = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_ang_old.frame.Y = (float*)calloc(pdbdata.atomCount, sizeof(float));
	dcd_ang_old.frame.Z = (float*)calloc(pdbdata.atomCount, sizeof(float));
	do{
        dcdReadFrame(&dcd_xyz);
        dcdReadFrame(&dcd_ang);
        frames++;
        t_xyz = 0.0;
        t_rot = 0.0;
        t_x = 0.0; t_y = 0.0; t_z = 0.0; tr_x = 0.0; tr_y = 0.0; tr_z = 0.0;
        if(frames == 1)
        {
            swap(dcd_xyz.frame.X, dcd_xyz_old.frame.X);
            swap(dcd_xyz.frame.Y, dcd_xyz_old.frame.Y);
            swap(dcd_xyz.frame.Z, dcd_xyz_old.frame.Z);
            swap(dcd_ang.frame.X, dcd_ang_old.frame.X);
            swap(dcd_ang.frame.Y, dcd_ang_old.frame.Y);
            swap(dcd_ang.frame.Z, dcd_ang_old.frame.Z);
            continue;
        }
        for(int i = 0; i < pdbdata.atomCount; i ++)
        {
            dx = dcd_xyz.frame.X[i] - dcd_xyz_old.frame.X[i];
            dy = dcd_xyz.frame.Y[i] - dcd_xyz_old.frame.Y[i];
            dz = dcd_xyz.frame.Z[i] - dcd_xyz_old.frame.Z[i];
            da = dcd_ang.frame.X[i] - dcd_ang_old.frame.X[i];
            db = dcd_ang.frame.Y[i] - dcd_ang_old.frame.Y[i];
            dg = dcd_ang.frame.Z[i] - dcd_ang_old.frame.Z[i];
            t_xyz += ( 
                pow(dx ,2) + 
                pow(dy, 2) +
                pow(dz, 2)
                ); 
            t_rot += (
                pow(da, 2) + 
                pow(db, 2) +
                pow(dg, 2) -
                2*(da)*(dg)*
                pow(cos(dcd_ang.frame.Y[i]),2)
                );
            t_x += pow(dx ,2);
            t_y += pow(dy ,2);
            t_z += pow(dz ,2);
            tr_x += pow(da,2);//pow(-da*sin(dcd_ang.frame.Y[i])*cos(dcd_ang.frame.Z[i]) + db*sin(dcd_ang.frame.Z[i]) ,2);
            tr_y += pow(db,2);//pow( da*sin(dcd_ang.frame.Y[i])*sin(dcd_ang.frame.Z[i]) + db*cos(dcd_ang.frame.Z[i]) ,2);
            tr_z += pow(dg,2);//pow(da*cos(dcd_ang.frame.Y[i]) + dg, 2);
        }

        float gammaR = 1.06e+06;
        float gammaT = 5e+06;
        int N = 1560;
        float K = 0.002;
//        t_xyz = t_xyz*1e-18/(pow((double)stride*2e-10,2)) /2 / 910;
        t_xyz *= (gammaR / (6*stride * 200 * N * K) );
        t_rot *= (gammaT / (6*stride * 200 * N * K) );
        t_x   *= (gammaR / (2*stride * 200 * N * K) );
        t_y   *= (gammaR / (2*stride * 200 * N * K) );
        t_z   *= (gammaR / (2*stride * 200 * N * K) );
        tr_x  *= (gammaT / (2*stride * 200 * N * K) );
        tr_y  *= (gammaT / (2*stride * 200 * N * K) );
        tr_z  *= (gammaT / (2*stride * 200 * N * K) );
        printf("%d %f %f %15f %f %f %15f %f %f\n", frames, t_xyz, t_rot, t_x, t_y, t_z, tr_x, tr_y, tr_z);
        t_rot = 0.0;
        t_xyz = 0.0;
        swap(dcd_xyz.frame.X, dcd_xyz_old.frame.X);
        swap(dcd_xyz.frame.Y, dcd_xyz_old.frame.Y);
        swap(dcd_xyz.frame.Z, dcd_xyz_old.frame.Z);
        swap(dcd_ang.frame.X, dcd_ang_old.frame.X);
        swap(dcd_ang.frame.Y, dcd_ang_old.frame.Y);
        swap(dcd_ang.frame.Z, dcd_ang_old.frame.Z);
	}while(!(feof(dcd_xyz.file)) && !(feof(dcd_ang.file)));
	return 0;
}
