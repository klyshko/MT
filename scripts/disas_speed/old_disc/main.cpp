#include "pdbio.h"
#include "dcdio.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <algorithm>
#define thres 10.12
#define pf_number 13

PDB pdbdata;
DCD dcd;

int main(int argc, char** argv){
		if(argc != 3){
				printf("Usage: %s infile.pdb infile.dcd\n", argv[0]);
				exit(1);
		}
		readPDB(argv[1], &pdbdata);
		dcdOpenRead(&dcd, argv[2]);
		dcdReadHeader(&dcd);
		long int frames =0;
		dcd.frame.X = (float*)calloc(pdbdata.atomCount, sizeof(float));
		dcd.frame.Y = (float*)calloc(pdbdata.atomCount, sizeof(float));
		dcd.frame.Z = (float*)calloc(pdbdata.atomCount, sizeof(float));

		int chain_lenghts[pf_number];
		int max_number[pf_number];
		do{ dcdReadFrame(&dcd); frames ++;

				std::fill_n(chain_lenghts, pf_number, 0);
				std::fill_n(max_number, pf_number, 0);
				for(int i=0; i < pdbdata.atomCount; i++) chain_lenghts[pdbdata.atoms[i].chain - 'A']++;
//				int current_pf = 0;

				for(int i=0; i < pdbdata.atomCount; i++){
//					printf("%f %d\n", dcd.frame.X[i], frames);
					if(dcd.frame.X[i] < thres){
							if( (max_number[pdbdata.atoms[i].chain - 'A'] < 2*pdbdata.atoms[i].resid) && 
									(fabs(dcd.frame.Z[i]) < 0.05) ){
									max_number[pdbdata.atoms[i].chain - 'A'] = 2*pdbdata.atoms[i].resid;
//									printf("%d %d\n",max_number[pdbdata.atoms[i].chain - 'A'], frames);
							}
							//		printf("%d %c\n",max_number[pdbdata.atoms[i].chain - 'A'], pdbdata.atoms[i].chain);
/*							max_number[current_pf] = 2*pdbdata.atoms[i].resid;
							current_pf ++;
							int offset = 0;
							for(int j = 0; j < current_pf; j++){
								offset += chain_lenghts[j];
							}
							i = offset-1;*/
					}
		}				
		int min_number = pdbdata.atomCount;																				
		for(int j = 0; j < pf_number; j++){
//			printf("%d %d %d\n",min_number,max_number[j], frames);
			if(min_number > max_number[j]) min_number = max_number[j];
		}
		if(min_number < 10){
				 break;
		}
		}while(!feof(dcd.file));
/*		int min = pdbdata.atomCount, max = 0;
		for(int i=0; i < pdbdata.atomCount; i++){
				if(max < pdbdata.atoms[i].resid)
						max = pdbdata.atoms[i].resid;
				if(dcd.frame.X[i] > thres){
						if(min > pdbdata.atoms[i].resid)
								min = pdbdata.atoms[i].resid;
				}
		}
		printf("%d %d %d\n", min, max, (max -min) * 2);*/
		int mean = 0;
		for(int i=0; i< pf_number; i++){chain_lenghts[i]=chain_lenghts[i]-max_number[i]; mean += chain_lenghts[i]; }
		printf("%f %d\n", (float)mean/pf_number, frames);
		return 0;
}
