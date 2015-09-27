/*
 * compute_cuda.cu
 *
 *  Created on: 04.06.2012
 *      Author: zhmurov
 */

#include <cuda.h>
#include "Cuda.h"
#include "mt.h"
#include "parameters.h"
#include "HybridTaus.cu"
#include "output.h"
#include "num_test.h"

#define BLOCK_SIZE 32
#define MAX_F 10.0

extern void update(long long int step);
extern void UpdateLJPairs();
extern void UpdatePairs();

__device__ __constant__ Parameters c_par;
__device__ __constant__ Topology c_top;

#ifdef T3LS
__device__ Coord d_F[1000][1000];
#endif


void init(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies){

	
}


#if defined(MORSE)
__device__ real dmorse(real D, real a, real x)
{
    return (2*a*D*(1-exp(-a*x))*exp(-a*x));
}

__device__ real morse_en(real D, real a, real x)
{
    return D * (1 - exp(-a*x)) * (1 - exp(-a*x)) - D;
}
#endif

#if defined(BARR)
__device__ real dbarr(real a, real r, real w, real x)
{
    return ( - a*exp(-(x-r)*(x-r)/(2*w*w)) * (x-r)/(w*w));
}

__device__ real barr(real a, real r, real w, real x)
{
    return a*exp(-(x-r)*(x-r)/(2*w*w));
}
#endif


__global__ void compute_kernel(const Coord* d_r, Coord* d_f){
	const int p = blockIdx.x*blockDim.x + threadIdx.x;
	const int ind = p%c_par.Ntot;
	const int traj = p/c_par.Ntot;
	real cos_fii, cos_fij, sin_fii, sin_fij, 
		  cos_psii, cos_psij, sin_psii, sin_psij,
		  cos_thetai, cos_thetaj, sin_thetai, sin_thetaj;
	real xi, xj, yi, yj, zi, zj;
	real dUdr, dr, gradx, grady, gradz, gradfi, gradpsi, gradtheta, expon, expon2, df, dft;
	int i,j;
	Coord ri, rj, fi = (Coord){0.0,0.0,0.0,0.0,0.0,0.0}, fj = (Coord){0.0,0.0,0.0,0.0,0.0,0.0};

	real xp1 = xp1_def;
	real yp1 = yp1_def;
	real zp1 = zp1_def;
	real xp2 = xp2_def;
	real yp2 = yp2_def;
	real zp2 = zp2_def;
	real R_MON = r_mon;

#ifdef T3LS
    for(int count=0; count<c_par.Ntot; count++)
    {
        d_F[ind][count].x     = 0.0f;
        d_F[ind][count].y     = 0.0f;
        d_F[ind][count].x     = 0.0f;
        d_F[ind][count].psi   = 0.0f;
        d_F[ind][count].fi    = 0.0f;
        d_F[ind][count].theta = 0.0f;
    }
#endif

	if(ind < c_par.Ntot && traj < c_par.Ntr){
		ri = d_r[p];
		cos_fii = cosf(ri.fi); 
		sin_fii = sinf(ri.fi);
		cos_psii = cosf(ri.psi);
		sin_psii = sinf(ri.psi);
		cos_thetai = cosf(ri.theta);
		sin_thetai = sinf(ri.theta);
		xi = ri.x;
		yi = ri.y;
		zi = ri.z;
		// harmonic
		for(int k = 0; k < c_top.harmonicCount[ind]; k++){
			j = c_top.harmonic[c_top.maxHarmonicPerMonomer*ind+k];
			if(j<0){
				R_MON = r_mon;
				j *= -1;
			}
			else{
				R_MON = -r_mon;
			}
			rj = d_r[j + traj*c_par.Ntot];
			cos_fij = cosf(rj.fi);
			sin_fij = sinf(rj.fi);
			cos_psij = cosf(rj.psi);
			sin_psij = sinf(rj.psi);
			cos_thetaj = cosf(rj.theta);
			sin_thetaj = sinf(rj.theta);
			xj = rj.x;
			yj = rj.y;
			zj = rj.z;
			dr = sqrtf(pow(-zi + zj -
		 		R_MON * cos_fii * cos_thetai - R_MON * cos_fij * cos_thetaj,2) +
	   			pow(-xi + xj -
		 		R_MON * (sin_fii * sin_psii + cos_fii * cos_psii * sin_thetai) -
				R_MON * (sin_fij * sin_psij + cos_fij * cos_psij * sin_thetaj),2) +
				pow(-yi + yj -
				R_MON * (-cos_psii * sin_fii + cos_fii * sin_psii * sin_thetai) -
				R_MON * (-cos_psij * sin_fij + cos_fij * sin_psij * sin_thetaj),2));

			gradx = -((-xi + xj - R_MON* (sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					R_MON* (sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj)));


			grady = -((-yi + yj -
					R_MON* (-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
					R_MON* (-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj)));


			gradz = -((-zi + zj - R_MON* cos_fii*cos_thetai - R_MON* cos_fij*cos_thetaj));


			gradtheta = ( R_MON*(-zi + zj - R_MON* cos_fii*cos_thetai - R_MON* cos_fij*cos_thetaj)*cos_fii*sin_thetai -
						R_MON* (-xi + xj -
						R_MON* (sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
						R_MON* (sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj))*cos_fii*cos_psii*cos_thetai -
						R_MON* (-yi + yj - R_MON* (-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
						R_MON* (-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj))*cos_fii*cos_thetai*sin_psii);


			gradpsi = (-R_MON*(-xi + xj - R_MON*(sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
						R_MON*(sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj))*(cos_psii*sin_fii - cos_fii*sin_psii*sin_thetai) -
						R_MON*(-yi + yj - R_MON*(-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
						R_MON*(-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj))*(sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai));


			gradfi = ( R_MON* (-zi + zj -
					R_MON* cos_fii*cos_thetai - R_MON* cos_fij*cos_thetaj)*cos_thetai*sin_fii -
					R_MON* (-xi + xj -
					R_MON* (sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					R_MON* (sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj))*(cos_fii*sin_psii - cos_psii*sin_fii*sin_thetai) -
					R_MON* (-yi + yj - R_MON* (-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
					R_MON* (-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj))*(-cos_fii*cos_psii - sin_fii*sin_psii*sin_thetai));

			fi.x     += -c_par.C * gradx;
			fi.y     += -c_par.C * grady;
			fi.z     += -c_par.C * gradz;
			fi.fi    += -c_par.C * gradfi;
			fi.psi   += -c_par.C * gradpsi;
			fi.theta += -c_par.C * gradtheta;

		
            if(dr < ANGLE_CUTOFF )
            {

                if(R_MON > 0){
                    fi.psi   += c_par.B_psi		*	(rj.psi 	-	ri.psi 		- c_par.psi_0	);
                    fi.fi	 += c_par.B_fi		*	(rj.fi 		-	ri.fi 		- c_par.fi_0	);
                    fi.theta += c_par.B_theta	*	(rj.theta 	- 	ri.theta 	- c_par.theta_0	);
                }
                else{
                    fi.psi   -= c_par.B_psi		*	(ri.psi 	- 	rj.psi 		- c_par.psi_0	);
                    fi.fi	 -= c_par.B_fi		*	(ri.fi 		- 	rj.fi 		- c_par.fi_0	);
                    fi.theta -= c_par.B_theta	*	(ri.theta 	- 	rj.theta 	- c_par.theta_0	);
                }
            }
            
		
		}
		
		for(int k = 0; k < c_top.longitudinalCount[ind]; k++){
			j = c_top.longitudinal[c_top.maxLongitudinalPerMonomer*ind+k];
			if(j < 0){
				R_MON = r_mon;
				j *= -1;
			}
			else{
				R_MON = -r_mon;
			}
			rj = d_r[j + traj*c_par.Ntot];
			cos_fij = cosf(rj.fi);
			sin_fij = sinf(rj.fi);
			cos_psij = cosf(rj.psi);
			sin_psij = sinf(rj.psi);
			cos_thetaj = cosf(rj.theta);
			sin_thetaj = sinf(rj.theta);
			xj = rj.x;
			yj = rj.y;
			zj = rj.z;
			dr = sqrtf(pow(-zi + zj -
		 		R_MON * cos_fii * cos_thetai - R_MON * cos_fij * cos_thetaj,2) +
	   			pow(-xi + xj -
		 		R_MON * (sin_fii * sin_psii + cos_fii * cos_psii * sin_thetai) -
				R_MON * (sin_fij * sin_psij + cos_fij * cos_psij * sin_thetaj),2) +
				pow(-yi + yj -
				R_MON * (-cos_psii * sin_fii + cos_fii * sin_psii * sin_thetai) -
				R_MON * (-cos_psij * sin_fij + cos_fij * sin_psij * sin_thetaj),2));
			
			gradx = -((-xi + xj - R_MON* (sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					R_MON* (sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj)));

			grady = -((-yi + yj -
					R_MON* (-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
					R_MON* (-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj)));


			gradz = -((-zi + zj - R_MON* cos_fii*cos_thetai - R_MON* cos_fij*cos_thetaj));


			gradtheta = ( R_MON* (-zi + zj -
						R_MON* cos_fii*cos_thetai - R_MON* cos_fij*cos_thetaj)*cos_fii*sin_thetai -
						R_MON* (-xi + xj -
						R_MON* (sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
						R_MON* (sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj))*cos_fii*cos_psii*cos_thetai -
						R_MON* (-yi + yj - R_MON* (-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
						R_MON* (-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj))*cos_fii*cos_thetai*sin_psii);


			gradpsi = (- R_MON* (-xi + xj -
						R_MON* (sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
						R_MON* (sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj))*(cos_psii*sin_fii - cos_fii*sin_psii*sin_thetai) -
						R_MON* (-yi + yj -
						R_MON* (-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
						R_MON* (-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj))*(sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai));


			gradfi = ( R_MON* (-zi + zj -
					R_MON* cos_fii*cos_thetai - R_MON* cos_fij*cos_thetaj)*cos_thetai*sin_fii -
					R_MON* (-xi + xj -
					R_MON* (sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					R_MON* (sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj))*(cos_fii*sin_psii - cos_psii*sin_fii*sin_thetai) -
					R_MON* (-yi + yj - R_MON* (-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
					R_MON* (-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj))*(-cos_fii*cos_psii - sin_fii*sin_psii*sin_thetai));


#if defined(MORSE)
            dUdr = dmorse(c_par.D_long, c_par.A_long, dr)/dr;
#endif

#if defined(BARR)
            dUdr += dbarr(c_par.a_barr_long, c_par.r_barr_long, c_par.w_barr_long, dr)/dr;
#endif
        
			fi.x     += -dUdr*gradx;
			fi.y     += -dUdr*grady;
			fi.z     += -dUdr*gradz;
			fi.fi    += -dUdr*gradfi;
			fi.psi   += -dUdr*gradpsi;
			fi.theta += -dUdr*gradtheta;

            if(dr < ANGLE_CUTOFF)
            {
                if(R_MON > 0){
                    fi.psi   += c_par.B_psi		*	(rj.psi 	- ri.psi 	- c_par.psi_0	);
                    fi.fi	 += c_par.B_fi		*	(rj.fi 		- ri.fi 	- c_par.fi_0	);
                    fi.theta += c_par.B_theta	*	(rj.theta 	- ri.theta 	- c_par.theta_0	);

                }
                else{
                    fi.psi   -= c_par.B_psi		*	(ri.psi 	- rj.psi 	- c_par.psi_0	);
                    fi.fi	 -= c_par.B_fi		*	(ri.fi 		- rj.fi 	- c_par.fi_0	);
                    fi.theta -= c_par.B_theta	*	(ri.theta 	- rj.theta 	- c_par.theta_0	);
                }
            }
            
		}

		
		for(int k = 0; k < 2; k++){//c_top.lateralCount[ind]; k++){
			j = c_top.lateral[2*ind + k];//c_top.maxLateralPerMonomer*ind+k];
			if (j != LARGENUMBER) {

				if(j <= 0){
					j*=-1;
					xp1 = xp2_def;
					yp1 = yp2_def;
					zp1 = zp2_def;
					xp2 = xp1_def;
					yp2 = yp1_def;
					zp2 = zp1_def;
				} else {
					xp1 = xp1_def;
					yp1 = yp1_def;
					zp1 = zp1_def;
					xp2 = xp2_def;
					yp2 = yp2_def;
					zp2 = zp2_def;
				}

				rj = d_r[j + traj*c_par.Ntot];
				cos_fij = cosf(rj.fi);
				sin_fij = sinf(rj.fi);
				cos_psij = cosf(rj.psi);
				sin_psij = sinf(rj.psi);
				cos_thetaj = cosf(rj.theta);
				sin_thetaj = sinf(rj.theta);
				xj = rj.x;
				yj = rj.y;
				zj = rj.z;
				
				dr = sqrtf(pow(zi - zj + zp2 * cos_fii * cos_thetai -
			 		zp1 * cos_fij * cos_thetaj + yp2 * cos_thetai * sin_fii -
			 		yp1 * cos_thetaj * sin_fij - xp2 * sin_thetai + xp1 * sin_thetaj,2) +
		   			pow(xi - xj - yp2 * cos_fii * sin_psii + zp2 * sin_fii * sin_psii +
			  		yp1 * cos_fij * sin_psij - zp1 * sin_fij * sin_psij +
				 	cos_psii * (xp2 * cos_thetai + zp2 * cos_fii * sin_thetai +
					yp2 * sin_fii * sin_thetai) - cos_psij * (xp1 * cos_thetaj + zp1 * cos_fij * sin_thetaj +
			   		yp1 * sin_fij * sin_thetaj),2) +
					pow(yi - yj - zp2 * cos_psii * sin_fii + zp1 * cos_psij * sin_fij +
					xp2 * cos_thetai * sin_psii - xp1 * cos_thetaj * sin_psij +
					yp2 * sin_fii * sin_psii * sin_thetai +
					cos_fii * (yp2 * cos_psii + zp2 * sin_psii * sin_thetai) -
					yp1 * sin_fij * sin_psij * sin_thetaj - cos_fij * (yp1 * cos_psij +
					zp1 * sin_psij * sin_thetaj),2));

				gradx = -((-xi + xj - xp2*cos_psii*cos_thetai +
					  xp1*cos_psij*cos_thetaj -
					  zp2*(sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					  yp2*(-cos_fii*sin_psii + cos_psii*sin_fii*sin_thetai) +
					  zp1*(sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj) +
					  yp1*(-cos_fij*sin_psij +
						 cos_psij*sin_fij*sin_thetaj)));


				grady = -((-yi + yj -
					  xp2*cos_thetai*sin_psii + xp1*cos_thetaj*sin_psij -
					  zp2*(-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
					  yp2*(cos_fii*cos_psii + sin_fii*sin_psii*sin_thetai) +
					  zp1*(-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj) +
					  yp1*(cos_fij*cos_psij +
						 sin_fij*sin_psij*sin_thetaj)));


				gradz = -((-zi + zj -
					  zp2*cos_fii*cos_thetai + zp1*cos_fij*cos_thetaj -
					  yp2*cos_thetai*sin_fii + yp1*cos_thetaj*sin_fij +
					  xp2*sin_thetai - xp1*sin_thetaj));


				gradtheta = ((-zi + zj - zp2*cos_fii*cos_thetai + zp1*cos_fij*cos_thetaj - yp2*cos_thetai*sin_fii + yp1*cos_thetaj*sin_fij +
					  xp2*sin_thetai - xp1*sin_thetaj)*(xp2*cos_thetai +
					  zp2*cos_fii*sin_thetai + yp2*sin_fii*sin_thetai) +
					 (-xi + xj - xp2*cos_psii*cos_thetai +
					  xp1*cos_psij*cos_thetaj -
					  zp2*(sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					  yp2*(-cos_fii*sin_psii + cos_psii*sin_fii*sin_thetai) +
					  zp1*(sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj) +
					  yp1*(-cos_fij*sin_psij + cos_psij*sin_fij*sin_thetaj))*(-zp2*cos_fii*cos_psii*cos_thetai - yp2*cos_psii*cos_thetai*sin_fii +
					  xp2*cos_psii*sin_thetai) +
					(-yi + yj - xp2*cos_thetai*sin_psii +
					  xp1*cos_thetaj*sin_psij -
					  zp2*(-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
					  yp2*(cos_fii*cos_psii + sin_fii*sin_psii*sin_thetai) +
					  zp1*(-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj) +
					  yp1*(cos_fij*cos_psij +
						 sin_fij*sin_psij*sin_thetaj))*(-zp2*cos_fii*cos_thetai*sin_psii - yp2*cos_thetai*sin_fii*sin_psii +
					  xp2*sin_psii*sin_thetai));


				gradpsi = ( (-xi + xj -
					  xp2*cos_psii*cos_thetai + xp1*cos_psij*cos_thetaj -
					  zp2*(sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					  yp2*(-cos_fii*sin_psii + cos_psii*sin_fii*sin_thetai) +
					  zp1*(sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj) +
					  yp1*(-cos_fij*sin_psij +
						 cos_psij*sin_fij*sin_thetaj))*(xp2*cos_thetai*sin_psii -
					   zp2*(cos_psii*sin_fii - cos_fii*sin_psii*sin_thetai) -
					  yp2*(-cos_fii*cos_psii -
						 sin_fii*sin_psii*sin_thetai)) +
					(-yi + yj - xp2*cos_thetai*sin_psii +
					  xp1*cos_thetaj*sin_psij -
					  zp2*(-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
					  yp2*(cos_fii*cos_psii + sin_fii*sin_psii*sin_thetai) +
					  zp1*(-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj) +
					  yp1*(cos_fij*cos_psij +
						 sin_fij*sin_psij*sin_thetaj))*(-xp2*cos_psii*cos_thetai -
					  zp2*(sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					  yp2*(-cos_fii*sin_psii +
						 cos_psii*sin_fii*sin_thetai)));


				gradfi = ( (-zi + zj -
					  zp2*cos_fii*cos_thetai + zp1*cos_fij*cos_thetaj -
					  yp2*cos_thetai*sin_fii + yp1*cos_thetaj*sin_fij +
					  xp2*sin_thetai - xp1*sin_thetaj)*(-yp2*cos_fii*cos_thetai +
					  zp2*cos_thetai*sin_fii) +
					(-xi + xj - xp2*cos_psii*cos_thetai +
					  xp1*cos_psij*cos_thetaj -
					  zp2*(sin_fii*sin_psii + cos_fii*cos_psii*sin_thetai) -
					  yp2*(-cos_fii*sin_psii + cos_psii*sin_fii*sin_thetai) +
					  zp1*(sin_fij*sin_psij + cos_fij*cos_psij*sin_thetaj) +
					  yp1*(-cos_fij*sin_psij +
						 cos_psij*sin_fij*sin_thetaj))*(-yp2*(sin_fii*sin_psii +
						 cos_fii*cos_psii*sin_thetai) -
					  zp2*(cos_fii*sin_psii -
						 cos_psii*sin_fii*sin_thetai)) +
					(-yi + yj - xp2*cos_thetai*sin_psii +
					  xp1*cos_thetaj*sin_psij -
					  zp2*(-cos_psii*sin_fii + cos_fii*sin_psii*sin_thetai) -
					  yp2*(cos_fii*cos_psii + sin_fii*sin_psii*sin_thetai) +
					  zp1*(-cos_psij*sin_fij + cos_fij*sin_psij*sin_thetaj) +
					  yp1*(cos_fij*cos_psij +
						 sin_fij*sin_psij*sin_thetaj))*(-yp2*(-cos_psii*sin_fii +
						  cos_fii*sin_psii*sin_thetai) -
					  zp2*(-cos_fii*cos_psii -
						 sin_fii*sin_psii*sin_thetai)));


#if defined(MORSE)
				if (c_top.harmonic[ind] * c_top.harmonic[j] < 0) {
					dUdr = dmorse(c_par.D_lat / 2, c_par.A_lat, dr)/dr;
				}
	            else {
	            	dUdr = dmorse(c_par.D_lat, c_par.A_lat, dr)/dr;
	            }
#endif	            

#if defined(BARR)
	            dUdr += dbarr(c_par.a_barr_lat, c_par.r_barr_lat, c_par.w_barr_lat, dr)/dr;
#endif		

				fi.x     += -dUdr*gradx;
				fi.y     += -dUdr*grady;
				fi.z     += -dUdr*gradz;
				fi.fi    += -dUdr*gradfi;
				fi.psi   += -dUdr*gradpsi;
				fi.theta += -dUdr*gradtheta;
			}

		}
		
		

#ifdef LJ_on
		for(int k=0; k < c_top.LJCount[p]; k++)
		{
			rj = d_r[c_top.LJ[k*c_par.Ntr*c_par.Ntot+p]];
			dr = sqrt(pow(ri.x-rj.x,2)+pow(ri.y-rj.y,2)+pow(ri.z-rj.z,2));


			if( dr < lj_cutoff )
            {
				fi.x += c_par.ljscale*c_par.ljsigma6*(6/pow(dr,8))*(ri.x-rj.x);
				fi.y += c_par.ljscale*c_par.ljsigma6*(6/pow(dr,8))*(ri.y-rj.y);
				fi.z += c_par.ljscale*c_par.ljsigma6*(6/pow(dr,8))*(ri.z-rj.z); 
			}
		}
#endif

#if defined(REPULSIVE)
    if (ri.z < - REP_H / 2){
        fi.z += 0.005 * REP_H;
    } else if (ri.z > REP_H) {
    	fi.z += -0.005 * REP_H;
    }
    real rad2 = ri.x * ri.x + ri.y * ri.y;

    if (rad2 > REP_R * REP_R){

        real coeff = -0.005 * (sqrt(rad2) - REP_R);
        fi.x += ri.x * coeff ;
        fi.y += ri.y * coeff;
    }
#endif

#if defined(BOUNDARIES)

    if (ri.x < xmin_bound) {
    	fi.x += -ks_bound * (x - xmin_bound);
    } else if (ri.x > xmax_bound) {
    	fi.x += -ks_bound * (x - xmax_bound);
    } 

    if (ri.y < ymin_bound) {
    	fi.y += -ks * (y - ymin_bound);
    } else if (ri.y > ymax_bound) {
    	fi.y += -ks_bound * (y - ymax_bound);
    } 

    if (ri.z < zmin_bound) {
    	fi.z += -ks_bound * (z - zmin_bound);
    } else if (ri.z > zmax_bound) {
    	fi.z += -ks_bound * (z - zmax_bound);
    } 
#endif
    	
		d_f[p] = fi;
		fi = (Coord){0.0,0.0,0.0,0.0,0.0,0.0};
		
	}
}

__global__ void pairs_kernel(const Coord* d_r){
	const int p = blockIdx.x*blockDim.x + threadIdx.x;
	const int i = p % c_par.Ntot;
	const int traj = p/c_par.Ntot;
	real cos_fii, cos_fij, sin_fii, sin_fij, 
		  cos_psii, cos_psij, sin_psii, sin_psij,
		  cos_thetai, cos_thetaj, sin_thetai, sin_thetaj;
	real xi, xj, yi, yj, zi, zj;
	real dr, dr2;
    Coord ri, rj;
    real xp1, yp1, zp1, xp2, yp2, zp2;
    real R_MON;

    if(i < c_par.Ntot && traj < c_par.Ntr){

    	c_top.lateral[2*i] = LARGENUMBER;
	    c_top.lateral[2*i + 1] = LARGENUMBER;
		
		if(c_top.harmonic[i] < 0)
	        R_MON = -r_mon;
	    else
	        R_MON = r_mon;
	        
	    ri = d_r[i];
	    cos_fii = cosf(ri.fi); 
	    sin_fii = sinf(ri.fi);
	    cos_psii = cosf(ri.psi);
	    sin_psii = sinf(ri.psi);
	    cos_thetai = cosf(ri.theta);
	    sin_thetai = sinf(ri.theta);
	    xi = ri.x;
	    yi = ri.y;
	    zi = ri.z;

	    real curMinDist = PAIR_CUTOFF;

	    for(int j = 0; j < c_par.Ntot; j++){
	        if(c_top.harmonic[i] * c_top.harmonic[j] <= 0){
	            rj = d_r[j];
	            cos_fij = cosf(rj.fi);
	            sin_fij = sinf(rj.fi);
	            cos_psij = cosf(rj.psi);
	            sin_psij = sinf(rj.psi);
	            cos_thetaj = cosf(rj.theta);
	            sin_thetaj = sinf(rj.theta);
	            xj = rj.x;
	            yj = rj.y;
	            zj = rj.z;
	            dr2 = pow(-zi + zj -
	                R_MON * cos_fii * cos_thetai - R_MON * cos_fij * cos_thetaj,2) +
	                pow(-xi + xj -
	                R_MON * (sin_fii * sin_psii + cos_fii * cos_psii * sin_thetai) -
	                R_MON * (sin_fij * sin_psij + cos_fij * cos_psij * sin_thetaj),2) +
	                pow(-yi + yj -
	                R_MON * (-cos_psii * sin_fii + cos_fii * sin_psii * sin_thetai) -
	                R_MON * (-cos_psij * sin_fij + cos_fij * sin_psij * sin_thetaj),2);
	            dr  = sqrt(dr2);

	            if(dr < curMinDist)
	            {
	                curMinDist = dr;
	                if(c_top.harmonic[i] < 0)
	                    c_top.longitudinal[i] = j;
	                else
	                    c_top.longitudinal[i] = -j;

	                c_top.longitudinalCount[i] = 1;
	            }
	        }
	    }

	    float curMinDistArr[2] = {PAIR_CUTOFF, PAIR_CUTOFF};
	    int latFlag[2] = {0, 0};

	    for(int j = 0; j < c_par.Ntot; j++){
	        if (i != j) {
	        	rj = d_r[j];
	            xj = rj.x;    
	            yj = rj.y;    
	            zj = rj.z;
	            sin_fij = sinf(rj.fi);
	            cos_fij = cosf(rj.fi);
	            sin_psij = sinf(rj.psi);
	            cos_psij = cosf(rj.psi);      
	            sin_thetaj = sinf(rj.theta);
	            cos_thetaj = cosf(rj.theta);  
	            for(int ind = 0; ind < 2; ind++){   
	                if(ind == 0){
	                    xp1 = xp2_def;
	                    yp1 = yp2_def;
	                    zp1 = zp2_def;
	                    xp2 = xp1_def;
	                    yp2 = yp1_def;
	                    zp2 = zp1_def;
	                } else {
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
	                
	                if (ind == 0) {
	                    if (dr < curMinDistArr[0]) {
	                        curMinDistArr[0] = dr;
	                        c_top.lateral[2 * i + 0] = -j;
	                        latFlag[0] = 1; //
	                    }
	                } else {
	                    if ((dr < curMinDistArr[1]) && (c_top.lateral[2 * i + 0] + j != 0)) {
	                        curMinDistArr[1] = dr;
	                        c_top.lateral[2 * i + 1] = j;
	                        latFlag[1] = 1; //
	                    }
	                }
	            }
	        }           
	    }
    }
    
}

__global__ void energy_kernel(const Coord* d_r, Energies* d_energies){
	
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	const int ind = p % c_par.Ntot;
	const int traj = p / c_par.Ntot;
	real cos_fii, cos_fij, sin_fii, sin_fij, 
		  cos_psii, cos_psij, sin_psii, sin_psij,
		  cos_thetai, cos_thetaj, sin_thetai, sin_thetaj;
	real xi, xj, yi, yj, zi, zj;
	int i,j;
	Coord ri, rj;
	real dr, dr2;
	Energies en; //= (Energies){0,0,0,0,0,0,0};
	real U_lat = 0.0, U_long = 0.0, U_harm = 0.0, U_fi = 0.0, U_psi = 0.0, U_teta = 0.0, U_lj = 0.0;

	real xp1 = xp1_def;
	real yp1 = yp1_def;
	real zp1 = zp1_def;
	real xp2 = xp2_def;
	real yp2 = yp2_def;
	real zp2 = zp2_def;
	real R_MON = r_mon;

	if (ind < c_par.Ntot && traj < c_par.Ntr){
		ri = d_r[p];
		cos_fii = cosf(ri.fi); 
		sin_fii = sinf(ri.fi);
		cos_psii = cosf(ri.psi);
		sin_psii = sinf(ri.psi);
		cos_thetai = cosf(ri.theta);
		sin_thetai = sinf(ri.theta);
		xi = ri.x;
		yi = ri.y;
		zi = ri.z;

		for(int k = 0; k < c_top.harmonicCount[ind]; k++){
			j = c_top.harmonic[c_top.maxHarmonicPerMonomer*ind+k];
			if (j < 0){
				R_MON = r_mon;
				j *= -1;
			} else {
				R_MON = -r_mon;
			}
			rj = d_r[j + traj * c_par.Ntot];
			cos_fij = cosf(rj.fi);
			sin_fij = sinf(rj.fi);
			cos_psij = cosf(rj.psi);
			sin_psij = sinf(rj.psi);
			cos_thetaj = cosf(rj.theta);
			sin_thetaj = sinf(rj.theta);
			xj = rj.x;
			yj = rj.y;
			zj = rj.z;
			dr = sqrt(pow(-zi + zj -
				R_MON * cos_fii * cos_thetai - R_MON * cos_fij * cos_thetaj,2) +
				pow(-xi + xj -
				R_MON * (sin_fii * sin_psii + cos_fii * cos_psii * sin_thetai) -
				R_MON * (sin_fij * sin_psij + cos_fij * cos_psij * sin_thetaj),2) +
				pow(-yi + yj -
				R_MON * (-cos_psii * sin_fii + cos_fii * sin_psii * sin_thetai) -
				R_MON * (-cos_psij * sin_fij + cos_fij * sin_psij * sin_thetaj),2));
#ifdef SIGMOID
            U_harm += c_par.A_xx + c_par.b_xx / (c_par.A_xx + exp(- c_par.c_xx * (dr - c_par.d_xx)));
#else
			U_harm += (c_par.C / 2) * pow(dr,2);
#endif
            if(dr < ANGLE_CUTOFF){
                if(R_MON > 0){
                    U_psi  	 += c_par.B_psi		/2	*	pow(rj.psi 		-	ri.psi 		- c_par.psi_0		,2);
                    U_fi	 += c_par.B_fi		/2	*	pow(rj.fi 		-	ri.fi 		- c_par.fi_0		,2);
                    U_teta	 += c_par.B_theta		/2	*	pow(rj.theta 	- 	ri.theta 	- c_par.theta_0	,2);
                }
                else{
                    U_psi  	 += c_par.B_psi	  	/2	*	pow(ri.psi 		- 	rj.psi 		- c_par.psi_0		,2);
                    U_fi	 += c_par.B_fi	  	/2	*	pow(ri.fi 		- 	rj.fi 		- c_par.fi_0		,2);
                    U_teta	 += c_par.B_theta		/2	*	pow(ri.theta 	- 	rj.theta 	- c_par.theta_0	,2);
                }
            }
		}
		for(int k = 0; k < c_top.longitudinalCount[ind]; k++){
			j = c_top.longitudinal[c_top.maxLongitudinalPerMonomer * ind + k];
			if (j < 0) {
				R_MON = r_mon;
				j *= -1;
			}
			else{
				R_MON = -r_mon;
			}
			rj = d_r[j + traj*c_par.Ntot];
			cos_fij = cosf(rj.fi);
			sin_fij = sinf(rj.fi);
			cos_psij = cosf(rj.psi);
			sin_psij = sinf(rj.psi);
			cos_thetaj = cosf(rj.theta);
			sin_thetaj = sinf(rj.theta);
			xj = rj.x;
			yj = rj.y;
			zj = rj.z;
			dr2 = pow(-zi + zj -
				R_MON * cos_fii * cos_thetai - R_MON * cos_fij * cos_thetaj,2) +
				pow(-xi + xj -
				R_MON * (sin_fii * sin_psii + cos_fii * cos_psii * sin_thetai) -
				R_MON * (sin_fij * sin_psij + cos_fij * cos_psij * sin_thetaj),2) +
				pow(-yi + yj -
				R_MON * (-cos_psii * sin_fii + cos_fii * sin_psii * sin_thetai) -
				R_MON * (-cos_psij * sin_fij + cos_fij * sin_psij * sin_thetaj),2);
			dr	= sqrt(dr2);

#if defined(MORSE)
            U_long += morse_en(c_par.D_long, c_par.A_long, dr);
#else
			U_long += (c_par.A_long*(c_par.b_long * dr2 * exp(-dr / c_par.r0_long) - c_par.c_long*exp(-dr2/(c_par.d_long*c_par.r0_long)))); 
#endif

#if defined(BARR)
            U_long += barr(c_par.a_barr_long, c_par.r_barr_long, c_par.w_barr_long, dr);
#endif
			if(dr < ANGLE_CUTOFF){
                if(R_MON > 0){
                    U_psi  	 += c_par.B_psi		/2	*	pow(rj.psi 		-	ri.psi 		- c_par.psi_0		,2);
                    U_fi	 += c_par.B_fi		/2	*	pow(rj.fi 		-	ri.fi 		- c_par.fi_0		,2);
                    U_teta	 += c_par.B_theta		/2	*	pow(rj.theta 	- 	ri.theta 	- c_par.theta_0	,2);
                }
                else{
                    U_psi  	 += c_par.B_psi	  	/2	*	pow(ri.psi 		- 	rj.psi 		- c_par.psi_0		,2);
                    U_fi	 += c_par.B_fi	  	/2	*	pow(ri.fi 		- 	rj.fi 		- c_par.fi_0		,2);
                    U_teta	 += c_par.B_theta		/2	*	pow(ri.theta 	- 	rj.theta 	- c_par.theta_0	,2);
                }
            }
		}

		for(int k = 0; k < c_top.lateralCount[ind]; k++){
			j = c_top.lateral[c_top.maxLateralPerMonomer * ind + k];
			if (j == LARGENUMBER) {
				break;
			}
			if (j <= 0){
				j *= -1;
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

			rj = d_r[j + traj*c_par.Ntot];
			cos_fij = cosf(rj.fi);
			sin_fij = sinf(rj.fi);
			cos_psij = cosf(rj.psi);
			sin_psij = sinf(rj.psi);
			cos_thetaj = cosf(rj.theta);
			sin_thetaj = sinf(rj.theta);
			xj = rj.x;
			yj = rj.y;
			zj = rj.z;
			
			dr = sqrt(pow(zi - zj + zp2 * cos_fii * cos_thetai -
				zp1 * cos_fij * cos_thetaj + yp2 * cos_thetai * sin_fii -
				yp1 * cos_thetaj * sin_fij - xp2 * sin_thetai + xp1 * sin_thetaj,2) +
				pow(xi - xj - yp2 * cos_fii * sin_psii + zp2 * sin_fii * sin_psii +
				yp1 * cos_fij * sin_psij - zp1 * sin_fij * sin_psij +
				cos_psii * (xp2 * cos_thetai + zp2 * cos_fii * sin_thetai +
				yp2 * sin_fii * sin_thetai) - cos_psij * (xp1 * cos_thetaj + zp1 * cos_fij * sin_thetaj +
				yp1 * sin_fij * sin_thetaj),2) +
				pow(yi - yj - zp2 * cos_psii * sin_fii + zp1 * cos_psij * sin_fij +
				xp2 * cos_thetai * sin_psii - xp1 * cos_thetaj * sin_psij +
				yp2 * sin_fii * sin_psii * sin_thetai +
				cos_fii * (yp2 * cos_psii + zp2 * sin_psii * sin_thetai) -
				yp1 * sin_fij * sin_psij * sin_thetaj - cos_fij * (yp1 * cos_psij +
				zp1 * sin_psij * sin_thetaj),2));
			dr2 = dr*dr;
		/*
            float c_lat, d_lat, b_lat, A_lat;
            int type = (top.mon_type[p]==0 && top.mon_type[j]==0) ? 0: //AA
                        ((top.mon_type[p]==0 && top.mon_type[j]==1)? 1 : //AB
                        ((top.mon_type[p]==1 && top.mon_type[j]==0)? 2 : 3)); //BA:BB
            A_lat = par.A_lat[type];
            b_lat = par.b_lat[type];
            c_lat = par.c_lat[type];
            d_lat = par.d_lat[type];
            U_lat += A_lat+b_lat/(A_lat+exp(-c_lat*(dr-d_lat)));
            */
#if defined(MORSE)
            if (c_top.harmonic[ind] * c_top.harmonic[j] < 0) {
				U_lat += morse_en(c_par.D_lat / 2, c_par.A_lat, dr); 	//dUdr = dmorse(c_par.D_lat / 2, c_par.A_lat, dr)/dr;
			}
            else {
            	U_lat += morse_en(c_par.D_lat, c_par.A_lat, dr);
            }
#else
			U_lat += (c_par.A_lat * (c_par.b_lat * dr2 * exp(-dr / c_par.r0_lat) - c_par.c_lat * exp(-dr2 / ( c_par.d_lat * c_par.r0_lat))));
#endif

#if defined(BARR)
            U_lat += barr(c_par.a_barr_lat, c_par.r_barr_lat, c_par.w_barr_lat, dr);
#endif
		}

#ifdef LJ_on
        for(int k = 0; k < c_top.LJCount[p]; k++){
            rj = d_r[c_top.LJ[k * c_par.Ntr * c_par.Ntot + p]];
            dr = sqrt(pow(ri.x - rj.x, 2) + pow(ri.y - rj.y, 2) + pow(ri.z - rj.z, 2));
            if(dr < lj_cutoff){
                U_lj += c_par.ljscale * c_par.ljsigma6 / pow(dr,6);
            }
        }
#endif       
	
        //U = U_harm + U_long + U_lat + U_lj + U_psi + U_fi + U_teta;
		en.U_harm = U_harm / 2;
		en.U_long = U_long / 2;
		en.U_lat = U_lat / 2;
		en.U_lj = U_lj / 2;
		en.U_psi = U_psi / 2;
		en.U_fi = U_fi / 2;
		en.U_teta = U_teta / 2;

		d_energies[p] = en;

	}	

}

__global__ void LJ_kernel(const Coord* r){

    const int p = blockIdx.x*blockDim.x + threadIdx.x;
    const int i = p % c_par.Ntot;
    const int tr = p/c_par.Ntot;
    real dr;
    Coord ri, rj;
    
    if(i < c_par.Ntot && tr < c_par.Ntr){

        ri = r[i];
        c_top.LJCount[i + tr * c_par.Ntot] = 0;

        for(int j = 0; j < c_par.Ntot; j++){
            rj = r[j];
            real dr = sqrt(pow(ri.x - rj.x,2)+
                            pow(ri.y - rj.y,2)+
                            pow(ri.z - rj.z,2));

            if((dr < c_par.ljpairscutoff) && (i != j)){
                c_top.LJCount[i + tr * c_par.Ntot]++;
                c_top.LJ[(c_top.LJCount[i + tr * c_par.Ntot] - 1) * c_par.Ntot * c_par.Ntr + c_par.Ntot * tr + i] = j + c_par.Ntot * tr;
            }
        }   
    }
}


__global__ void integrate_kernel(Coord* d_r, Coord* d_f){
	const int p = blockIdx.x*blockDim.x + threadIdx.x;
	float4 rf_xyz = make_float4(0,0,0,0);
	float4 rf_ang = make_float4(0,0,0,0);
	if(p < c_par.Ntot * c_par.Ntr){
		Coord f, ri;
		if(!c_top.fixed[p%c_par.Ntot]){
			f = d_f[p];
			ri = d_r[p];
			rf_xyz = rforce(p);
			rf_ang = rforce(p + c_par.Ntot*c_par.Ntr);
			ri.x += (c_par.dt/c_par.gammaR)*f.x + c_par.varR*rf_xyz.x;
			ri.y += (c_par.dt/c_par.gammaR)*f.y + c_par.varR*rf_xyz.y;
			ri.z += (c_par.dt/c_par.gammaR)*f.z + c_par.varR*rf_xyz.z;
#ifndef REDUCE_TO_2D
            // Disallow rotations in all directions but theta
			ri.fi    += (c_par.dt/c_par.gammaTheta)*f.fi  + c_par.varTheta*rf_ang.x;
			ri.psi   += (c_par.dt/c_par.gammaTheta)*f.psi + c_par.varTheta*rf_ang.y;
#endif
			ri.theta += (c_par.dt/c_par.gammaTheta)*f.theta + c_par.varTheta*rf_ang.z;
#ifdef REDUCE_TO_2D
            // Here, we move monomers back to respective PF plane
            // It is very kludgy solution. Probably, the kludgiest of all solutions ever existed
            // But I don't actually care. It was introduced for testing purposes only
            const int pf = (p % c_par.Ntot) / (c_par.Ntot / 13); // PF index = monomer_index / monomers_per_pf
//            printf("%d %d %d\n", p, p % c_par.Ntot, pf);
            const float a = -2*3.14159 * pf / 13;
//            const float rxy = sqrtf(ri.x*ri.x + ri.y*ri.y);
            const float rcosd = ri.y * sinf(a) + ri.x * cosf(a);  // = cos(phi-a) * rxy
            ri.x = rcosd * cosf(a);
            ri.y = rcosd * sinf(a);
#endif
			d_r[p] = ri;
		}
		f.x = 0.0f;
		f.y = 0.0f;
		f.z = 0.0f;
		f.fi = 0.0f;
		f.psi = 0.0f;
		f.theta = 0.0f;
		d_f[p] = f;

#ifdef T3LS
        for(int j=0; j<c_par.Ntot; j++)
        {
           printf("%*d%*d%*f%*f%*f%*f%*f%*f\n", 5, p, 5, j, 16, d_F[p][j].x, 16, d_F[p][j].y, 16, d_F[p][j].z,
                                                             16, d_F[p][j].psi, 16, d_F[p][j].fi, 16, d_F[p][j].theta);
        }
#endif
	}
}

void compute(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies){

	Coord* d_r;
	Coord* d_f;
	Topology topGPU;

	cudaSetDevice(par.device);
	checkCUDAError("device");
#ifdef T3LS
    cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 1024*1024*1024);
#endif

	cudaMalloc((void**)&d_r, par.Ntot*par.Ntr*sizeof(Coord));
	checkCUDAError("d_r allocation");
	cudaMalloc((void**)&d_f, par.Ntot*par.Ntr*sizeof(Coord));
	checkCUDAError("d_f allocation");

	// security from jerk
	for(int i = 0; i < par.Ntot*par.Ntr; i++){
		f[i].fi = 0.0f;
		f[i].psi = 0.0f;
		f[i].theta = 0.0f;
		f[i].x = 0.0f;
		f[i].y = 0.0f;
		f[i].z = 0.0f;
	}

	cudaMemcpy(d_f, f, par.Ntot*par.Ntr*sizeof(Coord), cudaMemcpyHostToDevice);
	checkCUDAError("copy_forces");
	cudaMemcpy(d_r, r, par.Ntot*par.Ntr*sizeof(Coord), cudaMemcpyHostToDevice);
	checkCUDAError("from r to d_r copy");

	topGPU.maxHarmonicPerMonomer = top.maxHarmonicPerMonomer;
	topGPU.maxLongitudinalPerMonomer = top.maxLongitudinalPerMonomer;
	topGPU.maxLateralPerMonomer = top.maxLateralPerMonomer;

	cudaMalloc((void**)&(topGPU.harmonicCount), par.Ntot*sizeof(int));
	checkCUDAError("topGPU.harmonicCount allocation");
	cudaMalloc((void**)&(topGPU.longitudinalCount), par.Ntot*par.Ntr*sizeof(int));
	checkCUDAError("topGPU.longitudinalCount allocation");
	cudaMalloc((void**)&(topGPU.lateralCount), par.Ntot*par.Ntr*sizeof(int));
	checkCUDAError("topGPU.lateralCount allocation");
	cudaMalloc((void**)&(topGPU.harmonic), par.Ntot*sizeof(int)*topGPU.maxHarmonicPerMonomer);
	checkCUDAError("harmonic allocation");
	cudaMalloc((void**)&(topGPU.longitudinal), par.Ntot*par.Ntr*sizeof(int)*topGPU.maxLongitudinalPerMonomer);
	checkCUDAError("long allocation");
	cudaMalloc((void**)&(topGPU.lateral), par.Ntot*par.Ntr*sizeof(int)*topGPU.maxLateralPerMonomer);
	checkCUDAError("lateral allocation");

	cudaMemcpy(topGPU.harmonic, top.harmonic, par.Ntot*topGPU.maxHarmonicPerMonomer*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("harmonic copy");
	cudaMemcpy(topGPU.longitudinal, top.longitudinal, par.Ntot*par.Ntr*topGPU.maxLongitudinalPerMonomer*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("long copy");
	cudaMemcpy(topGPU.lateral, top.lateral, par.Ntot*par.Ntr*topGPU.maxLateralPerMonomer*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("lateral copy");
	cudaMemcpy(topGPU.harmonicCount, top.harmonicCount, par.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("harmonic count copy");
	cudaMemcpy(topGPU.longitudinalCount, top.longitudinalCount, par.Ntot*par.Ntr*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("long count copy");
	cudaMemcpy(topGPU.lateralCount, top.lateralCount, par.Ntot*par.Ntr*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("lat count copy");

	cudaMalloc((void**)&(topGPU.fixed), par.Ntot*sizeof(bool));
	checkCUDAError("topGPU.fixed allocation");
	cudaMalloc((void**)&(topGPU.mon_type), par.Ntot*sizeof(int));
	checkCUDAError("d_mon_type allocation");
	cudaMemcpy(topGPU.fixed, top.fixed, par.Ntot*sizeof(bool), cudaMemcpyHostToDevice);
	checkCUDAError("fixed copy");
	cudaMemcpy(topGPU.mon_type, top.mon_type, par.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("montype copy");

#ifdef LJ_on
	topGPU.maxLJPerMonomer = 256;  ///I hope it will be enough       =top.maxLJPerMonomer;
	cudaMalloc((void**)&(topGPU.LJCount), par.Ntot*sizeof(int)*par.Ntr);
	checkCUDAError("lj_count allocation");
	cudaMalloc((void**)&(topGPU.LJ), par.Ntot*sizeof(int)*par.Ntr*topGPU.maxLJPerMonomer);
	checkCUDAError("lj allocation");
	/*
	cudaMemcpy(topGPU.LJCount, top.LJCount, par.Ntot*par.Ntr*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("lj count copy");
	cudaMemcpy(topGPU.LJ, top.LJ, par.Ntot*par.Ntr*topGPU.maxLJPerMonomer*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("lj copy");
	*/
#endif
	//const memory
	cudaMemcpyToSymbol(c_top, &topGPU, sizeof(Topology), 0, cudaMemcpyHostToDevice);
	checkCUDAError("copy of topGPU pointer to const memory");

	cudaMemcpyToSymbol(c_par, &par, sizeof(Parameters), 0, cudaMemcpyHostToDevice);
	checkCUDAError("copy parameters");

#ifdef OUTPUT_EN     //energies initializing
	Energies* d_energies;
	cudaMalloc((void**)&d_energies, par.Ntot*par.Ntr * sizeof(Energies));
	checkCUDAError("d_energies allocation");
#endif
	
	initRand(par.rseed, 2*par.Ntot*par.Ntr);
																
/*
|	--------   KERNEL CALLS
*/																
	for(long long int step = 0; step < par.steps; step++){

#ifdef LJ_on
		if(step % par.ljpairsupdatefreq == 0){ //pairs update frequency 

			LJ_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r);
			checkCUDAError("lj_kernel");
#if defined(ASSEMBLY)
            pairs_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r);
            checkCUDAError("pairs_kernel");
            //printf("pairs updated\n");
#endif
		}
#endif
		if(step % par.stride == 0){ //every stride steps do energy computing and outputing DCD
#ifdef OUTPUT_EN
            energy_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r, d_energies);
            checkCUDAError("energy_kernel");
            cudaMemcpy(energies, d_energies, par.Ntr * par.Ntot * sizeof(Energies), cudaMemcpyDeviceToHost);
            checkCUDAError("energy_copy");
#endif
			cudaMemcpy(r, d_r, par.Ntr*par.Ntot*sizeof(Coord), cudaMemcpyDeviceToHost);
			checkCUDAError("r update copy");
			update(step);
		}

		compute_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r, d_f);
		checkCUDAError("compute_forces");

#if defined (TEST_OPT) || defined (OUTPUT_FORCE)
		cudaMemcpy(f, d_f, par.Ntot*par.Ntr*sizeof(Coord), cudaMemcpyDeviceToHost);
#endif

		integrate_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r, d_f);
		checkCUDAError("integrate_kernel");
	}

	cudaFree(d_r);
	cudaFree(d_f);
	cudaFree(topGPU.harmonicCount);
	cudaFree(topGPU.longitudinalCount);
	cudaFree(topGPU.lateralCount);
	cudaFree(topGPU.fixed);
	cudaFree(topGPU.mon_type);
	cudaFree(topGPU.harmonic);
	cudaFree(topGPU.longitudinal);
	cudaFree(topGPU.lateral);
#ifdef LJ_on
	cudaFree(topGPU.LJCount);
	cudaFree(topGPU.LJ);
#endif
#ifdef OUTPUT_EN
	cudaFree(d_energies);
#endif
	checkCUDAError("cleanup");
}
