#include "pdbio.h"
#include "dcdio.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#define thres 5.0
#define hor_thres 2.0
#define theta_thres 0.2
#define pf_number 13

PDB pdbdata;
DCD dcd;

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
		if(argc < 3){
				printf("Usage: %s infile.pdb infile.dcd [hist|timeline]\n", argv[0]);
				exit(1);
		}
		readPDB(argv[1], &pdbdata);
		dcdOpenRead(&dcd, argv[2]);
		printf("%s\n", argv[2]);
		dcdReadHeader(&dcd);
		long int frames = 0;
		long int frames_curled = 0;
		bool start_curl_record = false;
		dcd.frame.X = (float*)calloc(pdbdata.atomCount, sizeof(float));
		dcd.frame.Y = (float*)calloc(pdbdata.atomCount, sizeof(float));
		dcd.frame.Z = (float*)calloc(pdbdata.atomCount, sizeof(float));

		int chain_lenghts[pf_number];
		int mt_end_number[pf_number];
		float mt_end_y[pf_number];
		float curl_min_y[pf_number];
		int pf_end_number[pf_number];
		int curled_start[pf_number];
		vector<int> curled_lenghts, hist;

		FILE* lfile;

		if( (argc==4) && (!strcmp(argv[3], "timeline"))){
			char filename[128];
			sprintf(filename, "%s.dat", argv[2]);
			lfile = fopen(filename,"w");
		}
		do{
			dcdReadFrame(&dcd); frames ++;
			std::fill_n(chain_lenghts, pf_number, 0);
			std::fill_n(mt_end_number, pf_number, 0);
			std::fill_n(mt_end_y, pf_number, 0);
			std::fill_n(curl_min_y, pf_number, pdbdata.atomCount*2.0f);
			std::fill_n(pf_end_number, pf_number, pdbdata.atomCount);
			std::fill_n(curled_start, pf_number, pdbdata.atomCount);
			for(int i=0; i < pdbdata.atomCount; i++) chain_lenghts[pdbdata.atoms[i].chain - 'A']++;

			for(int i=0; i < pdbdata.atomCount; i++)
			{
				for(int j=0; j < pdbdata.atomCount; j++)
				{
					int resi = pdbdata.atoms[i].resid;
					int resj = pdbdata.atoms[j].resid;
					int idi = pdbdata.atoms[i].id;
					int idj = pdbdata.atoms[j].id;
					char chaini = pdbdata.atoms[i].chain;
					char chainj = pdbdata.atoms[j].chain;
					char* namei = pdbdata.atoms[i].name;
					char* namej = pdbdata.atoms[j].name;
					if( (abs(resi - resj) == 1) && (chaini==chainj) && (sqrt(pow(dcd.frame.X[i] - dcd.frame.X[j],2) + pow(dcd.frame.Y[i] - dcd.frame.Y[j],2)) > thres) 
							&& (pf_end_number[chaini - 'A'] > min(resi, resj)) && (namei[1]!=namej[1]) && (abs(idi-idj)==1) )
					{
							pf_end_number[chaini - 'A'] = min(resi, resj);
					}
				}
			}
			for(int i = 0; i < pf_number; i++)
			{
				if(pf_end_number[i] > chain_lenghts[i]/2)
				{
						pf_end_number[i] = chain_lenghts[i]/2;
				}
			}
		
			//for(int i = 0; i<pf_number; i++)printf("MAX_X %f\n", mt_max_x[i]);
			for(int i = 0; i < pdbdata.atomCount; i++)
			{
				if( (pdbdata.atoms[i].resid <= pf_end_number[pdbdata.atoms[i].chain - 'A']) && (dcd.frame.Y[i] > mt_end_y[pdbdata.atoms[i].chain - 'A']) && (dcd.frame.Z[i] < hor_thres) ) 
				{
					mt_end_y[pdbdata.atoms[i].chain - 'A'] = dcd.frame.Y[i];
					mt_end_number[pdbdata.atoms[i].chain - 'A'] = pdbdata.atoms[i].resid;
				}
				if( (dcd.frame.Z[i] > theta_thres) && (pdbdata.atoms[i].resid < pf_end_number[pdbdata.atoms[i].chain - 'A']) && (pdbdata.atoms[i].resid < curled_start[pdbdata.atoms[i].chain - 'A']) )
						curled_start[pdbdata.atoms[i].chain - 'A'] = pdbdata.atoms[i].resid;
			}
			int min_number = pdbdata.atomCount;	
			int max_pf_number = 0;
			for(int j = 0; j < pf_number; j++)
			{
				if(max_pf_number < chain_lenghts[j]/2)
				{
					max_pf_number = chain_lenghts[j]/2;
				}
			}
			for(int j = 0; j < pf_number; j++)
			{
				if(curled_start[j] > pf_end_number[j])
				{
					curled_start[j] = pf_end_number[j];
				}
				if(min_number > mt_end_number[j])
				{
					 min_number = mt_end_number[j];
				}
			}
			int max_curl_length=0;
			for(int i=0; i < pf_number; i++)
			{
				if( max_curl_length < (pf_end_number[i]-curled_start[i])) max_curl_length = (pf_end_number[i]-curled_start[i]);
			}
			if( ((max_pf_number - min_number) >= 10) || (max_curl_length < 2)) start_curl_record = true;
			for(int i=0; i < pf_number; i++)
			{
				if(start_curl_record)
				{
					 curled_lenghts.push_back((pf_end_number[i]-curled_start[i]));
				}
			}
			if(start_curl_record) frames_curled++;
			if(min_number < 12)
			{
					 break;
			}
			if( (argc==4) && (!strcmp(argv[3], "timeline")) )
			{
				float lt = 0;
				for(int i=0; i < pf_number; i++) lt += (float)mt_end_number[i]/13.0;
				fprintf(lfile, "%ld %f\n", frames, 2*lt);
			}
		} while(!feof(dcd.file));

		//fclose(lfile);
		int mean = 0;
		int mean_curled = 0;
		for(int i=0; i< pf_number; i++)
		{
				chain_lenghts[i]=chain_lenghts[i]/2-mt_end_number[i]; mean += chain_lenghts[i]; 
				curled_start[i]=pf_end_number[i]-curled_start[i]; mean_curled += curled_start[i];
		}
		if(argc == 3)
		{
			printf("%f %ld %f", (float)2*mean/pf_number, frames, (float)2*mean_curled/pf_number);
		}
		else if( (argc==4) && (!strcmp(argv[3], "hist")) )
		{
			if(!start_curl_record) return 0;
			int max = *max_element(curled_lenghts.begin(), curled_lenghts.end());
			hist.resize(max+1,0);
			for(vector<int>::iterator it = curled_lenghts.begin(); it != curled_lenghts.end(); ++it){
				hist[*it]+=1;
			}
			printf("%ld ", frames_curled);
			for(vector<int>::iterator it = hist.begin(); it != hist.end(); ++it)
				printf("%f ", (float)*it/frames_curled);
			printf("\n");
		}
		return 0;
}
