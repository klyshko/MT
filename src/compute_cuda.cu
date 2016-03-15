/*
 * compute_cuda.cu
 *
 *  Created on: 04.06.2012
 *      Author: zhmurov

   Modified on: 09.01.2015
   		Author: klyshko
 */

#include <cuda.h>
#include "Cuda.h"
#include "mt.h"
#include "parameters.h"
#include "HybridTaus.cu"

#define BLOCK_SIZE 32

extern void update(long long int step, int* mt_len);
extern int change_conc(int* delta, int* mt_len);
extern void mt_length(long long int step, int* mt_len);
extern void UpdateLJPairs();
extern void UpdatePairs();

__device__ __constant__ Parameters c_par;
__device__ __constant__ Topology c_top;


__device__ real dmorse(real D, real a, real x){
    return (2*a*D*(1-exp(-a*x))*exp(-a*x));
}

__device__ real morse_en(real D, real a, real x){
    return D * (1 - exp(-a*x)) * (1 - exp(-a*x)) - D;
}

__device__ real dbarr(real a, real r, real w, real x){
    return ( - a*exp(-(x-r)*(x-r)/(2*w*w)) * (x-r)/(w*w));
}

__device__ real barr(real a, real r, real w, real x){
    return a*exp(-(x-r)*(x-r)/(2*w*w));
}


__global__ void compute_kernel(const Coord* d_r, Coord* d_f){
	const int p = blockIdx.x*blockDim.x + threadIdx.x;
	const int ind = p%c_par.Ntot;
	const int traj = p/c_par.Ntot;
	real cos_fii, cos_fij, sin_fii, sin_fij, 
		  cos_psii, cos_psij, sin_psii, sin_psij,
		  cos_thetai, cos_thetaj, sin_thetai, sin_thetaj;
	real xi, xj, yi, yj, zi, zj;
	real dUdr, dr, gradx, grady, gradz, gradfi, gradpsi, gradtheta;
	int j;
	real psiji, thetaji, fiji, psiij, thetaij, fiij;
	Coord ri, rj, fi = (Coord){0.0,0.0,0.0,0.0,0.0,0.0};

	real xp1 = xp1_def;
	real yp1 = yp1_def;
	real zp1 = zp1_def;
	real xp2 = xp2_def;
	real yp2 = yp2_def;
	real zp2 = zp2_def;
	real R_MON = r_mon;

	if(ind < c_par.Ntot && traj < c_par.Ntr){
		
		if (!c_top.extra[ind + traj * c_par.Ntot]){

			
			ri = d_r[p];
			xi = ri.x;
			yi = ri.y;
			zi = ri.z; 
			cos_fii = cosf(ri.fi); 
			sin_fii = sinf(ri.fi);
			cos_psii = cosf(ri.psi);
			sin_psii = sinf(ri.psi);
			cos_thetai = cosf(ri.theta);
			sin_thetai = sinf(ri.theta);
			
			// harmonic
			for(int k = 0; k < c_top.harmonicCount[ind]; k++){
				j = c_top.harmonic[c_top.maxHarmonicPerMonomer*ind+k];
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
	            	psiji = rj.psi - ri.psi - 2 * M_PI * (int)((rj.psi - ri.psi)/(2 * M_PI));
	                psiij = - psiji;
	            	thetaji = rj.theta - ri.theta - 2 * M_PI * (int)((rj.theta - ri.theta)/(2 * M_PI));
	            	thetaij = - thetaji;
	              	fiji = rj.fi - ri.fi - 2 * M_PI * (int)((rj.fi - ri.fi)/(2 * M_PI));
	                fiij = - fiji;
	                
	                if(R_MON > 0){

	                    fi.psi   += c_par.B_psi		*	(psiji 		- c_par.psi_0	);
	                    fi.fi	 += c_par.B_fi		*	(fiji 		- c_par.fi_0	);
	                    fi.theta += c_par.B_theta	*	(thetaji 	- c_par.theta_0	);
	                }
	                else{
	                    fi.psi   -= c_par.B_psi		*	(psiij	- 	c_par.psi_0	);
	                    fi.fi	 -= c_par.B_fi		*	(fiij 		- c_par.fi_0	);
	                    fi.theta -= c_par.B_theta	*	(thetaij 	- c_par.theta_0	);
	                }
	                
	              	
	            }
	            
			}
			
			
			
			
#if defined(MORSE)
			
			for(int k = 0; k < c_top.longitudinalCount[ind + c_par.Ntot * traj]; k++){
				j = c_top.longitudinal[c_top.maxLongitudinalPerMonomer * c_par.Ntot * traj + ind * c_top.maxLongitudinalPerMonomer + k];
				if(j < 0){
					R_MON = r_mon;
					j = abs(j);
				}
				else {
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



				if (dr == 0) dUdr = 0.0;
				else dUdr = dmorse(c_par.D_long, c_par.A_long, dr) / dr;



#if defined(BARR)
				if (dr != 0) 
	            dUdr += dbarr(c_par.a_barr_long, c_par.r_barr_long, c_par.w_barr_long, dr) / dr;
#endif
	        
				fi.x     += -dUdr*gradx;
				fi.y     += -dUdr*grady;
				fi.z     += -dUdr*gradz;
				fi.fi    += -dUdr*gradfi;
				fi.psi   += -dUdr*gradpsi;
				fi.theta += -dUdr*gradtheta;
				
	            if(dr < ANGLE_CUTOFF )
	            {
	            	
	            	psiji = rj.psi - ri.psi - 2 * M_PI * (int)((rj.psi - ri.psi)/(2 * M_PI));
	                psiij = - psiji;
	            	thetaji = rj.theta - ri.theta - 2 * M_PI * (int)((rj.theta - ri.theta)/(2 * M_PI));
	            	thetaij = - thetaji;
	              	fiji = rj.fi - ri.fi - 2 * M_PI * (int)((rj.fi - ri.fi)/(2 * M_PI));
	                fiij = - fiji;
	                
	                 if(R_MON > 0){

	                    fi.psi   += c_par.B_psi		*	(psiji 		- c_par.psi_0	);
	                    fi.fi	 += c_par.B_fi		*	(fiji 		- c_par.fi_0	);
	                    fi.theta += c_par.B_theta	*	(thetaji 	- c_par.theta_0	);
	                }
	                else{
	                    fi.psi   -= c_par.B_psi		*	(psiij	- 	c_par.psi_0	);
	                    fi.fi	 -= c_par.B_fi		*	(fiij 		- c_par.fi_0	);
	                    fi.theta -= c_par.B_theta	*	(thetaij 	- c_par.theta_0	);
	                }
	              	/*
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
	                */
	            }
	            
			} 

#endif

#if defined(MORSE)
			for(int k = 0; k < c_top.lateralCount[ind + traj * c_par.Ntot]; k++){
				j = c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + ind * c_top.maxLateralPerMonomer + k];//c_top.maxLateralPerMonomer*ind+k];
				if (abs(j) != LARGENUMBER) {

					if(j <= 0){
						j = abs(j);
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



					if (dr == 0) dUdr = 0.0;
					else if (c_top.mon_type[ind] != c_top.mon_type[j]) {
						dUdr = dmorse(c_par.D_lat / 2, c_par.A_lat, dr) / dr;
					}
		            else {
		            	dUdr = dmorse(c_par.D_lat, c_par.A_lat, dr) / dr;
		            }
		            

#if defined(BARR)
				if (dr != 0) 
	            dUdr += dbarr(c_par.a_barr_long, c_par.r_barr_long, c_par.w_barr_long, dr) / dr;
#endif

					fi.x     += -dUdr*gradx;
					fi.y     += -dUdr*grady;
					fi.z     += -dUdr*gradz;
					fi.fi    += -dUdr*gradfi;
					fi.psi   += -dUdr*gradpsi;
					fi.theta += -dUdr*gradtheta;
				}

			}
#endif


#ifdef LJ_on

						
			for(int k = 0; k < c_top.LJCount[ind + traj * c_par.Ntot]; k++){

				j = c_top.LJ[c_top.maxLJPerMonomer * c_par.Ntot * traj + ind * c_top.maxLJPerMonomer + k];
				rj = d_r[j + traj * c_par.Ntot];
				dr = sqrt(pow(ri.x-rj.x,2)+pow(ri.y-rj.y,2)+pow(ri.z-rj.z,2));
				
				if( dr < lj_cutoff )
	            {	
	            	real df = 6/pow(dr,8);

	            #ifdef REGULARIZATION
	            	if (c_par.ljscale * c_par.ljsigma6 * df * dr > c_par.gammaR * r_mon / c_par.dt ){
	            		df = c_par.gammaR * r_mon / (c_par.dt * dr * c_par.ljscale * c_par.ljsigma6);	
	            	} 
				#endif	
	            	
	            	fi.x += c_par.ljscale*c_par.ljsigma6*df*(ri.x-rj.x);
					fi.y += c_par.ljscale*c_par.ljsigma6*df*(ri.y-rj.y);
					fi.z += c_par.ljscale*c_par.ljsigma6*df*(ri.z-rj.z); 

				}
			}


#endif

#if defined(REPULSIVE)

			
		    if (ri.z < c_par.rep_leftborder){ 

		        fi.z += c_par.rep_eps * fabs(ri.z - c_par.rep_leftborder);
		    } else if (ri.z > c_par.zs[traj] + c_par.rep_leftborder) {

		    	fi.z += - c_par.rep_eps * fabs(ri.z - (c_par.zs[traj] + c_par.rep_leftborder) );
		    }

		    //fi.z += - c_par.rep_eps * fabs(ri.z + 100); // traj = 2

		    real rad2 = ri.x * ri.x + ri.y * ri.y;
		    if (rad2 > c_par.rep_r * c_par.rep_r){

		        real coeff = -c_par.rep_eps * (sqrt(rad2) - c_par.rep_r);
		        fi.x += ri.x * coeff ;
		        fi.y += ri.y * coeff;
		    }
		    
#endif

			d_f[p] = fi;
			fi = (Coord){0.0,0.0,0.0,0.0,0.0,0.0};
			
		}
	}
}

__global__ void pairs_kernel(const Coord* d_r){
	const int p = blockIdx.x * blockDim.x + threadIdx.x;
	const int i = p % c_par.Ntot;
	const int traj = p / c_par.Ntot;
	real cos_fii, cos_fij, sin_fii, sin_fij, 
		  cos_psii, cos_psij, sin_psii, sin_psij,
		  cos_thetai, cos_thetaj, sin_thetai, sin_thetaj;
	real xi, xj, yi, yj, zi, zj;
	real dr, dr2;
    Coord ri, rj;
    real xp1, yp1, zp1, xp2, yp2, zp2;
    real R_MON;
    
    if(i < c_par.Ntot && traj < c_par.Ntr){

    	c_top.lateralCount[i + traj * c_par.Ntot] = 0;
    	c_top.longitudinalCount[i + traj * c_par.Ntot] = 0;

    	if(!c_top.extra[i + traj * c_par.Ntot]){
			
			if(c_top.harmonic[c_top.maxHarmonicPerMonomer * i ] < 0) /// c_top.maxHarmonicPerMonomer*ind+k
		        R_MON = -r_mon; 
		    else
		        R_MON = r_mon;
		        
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

		    for(int j = 0; j < c_par.Ntot; j++){
		        if(c_top.mon_type[i] != c_top.mon_type[j] && abs(c_top.harmonic[c_top.maxHarmonicPerMonomer * i]) != j ){ // C_top.mon_type
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
		            dr2 = pow(-zi + zj -
		                R_MON * cos_fii * cos_thetai - R_MON * cos_fij * cos_thetaj,2) +
		                pow(-xi + xj -
		                R_MON * (sin_fii * sin_psii + cos_fii * cos_psii * sin_thetai) -
		                R_MON * (sin_fij * sin_psij + cos_fij * cos_psij * sin_thetaj),2) +
		                pow(-yi + yj -
		                R_MON * (-cos_psii * sin_fii + cos_fii * sin_psii * sin_thetai) -
		                R_MON * (-cos_psij * sin_fij + cos_fij * sin_psij * sin_thetaj),2);
		            dr  = sqrt(dr2);

		            if(dr < PAIR_CUTOFF){
		            	
		            	if(c_top.harmonic[c_top.maxHarmonicPerMonomer * i] < 0)
		                    c_top.longitudinal[c_top.maxLongitudinalPerMonomer * c_par.Ntot * traj + i * c_top.maxLongitudinalPerMonomer + c_top.longitudinalCount[i + traj * c_par.Ntot]] = j;
		                else
		                    c_top.longitudinal[c_top.maxLongitudinalPerMonomer * c_par.Ntot * traj + i * c_top.maxLongitudinalPerMonomer + c_top.longitudinalCount[i + traj * c_par.Ntot]] = -j;
		                c_top.longitudinalCount[i + traj * c_par.Ntot]++;
		            }
		            /*
		            if(dr < curMinDist)
		            {
		                curMinDist = dr;
		                if(c_top.harmonic[i] < 0)
		                    c_top.longitudinal[c_top.maxLongitudinalPerMonomer * c_par.Ntot * traj + i * c_top.maxLongitudinalPerMonomer] = j;
		                else
		                    c_top.longitudinal[c_top.maxLongitudinalPerMonomer * c_par.Ntot * traj + i * c_top.maxLongitudinalPerMonomer] = -j;

		                c_top.longitudinalCount[i + traj * c_par.Ntot] = 1;
		            }
		            */
		        }
		    }

		    //float curMinDistArr[2] = {PAIR_CUTOFF, PAIR_CUTOFF};
		    
		    for(int j = 0; j < c_par.Ntot; j++){
		        if (i != j && abs(c_top.harmonic[c_top.maxHarmonicPerMonomer * i]) != j){// && c_top.mon_type[i] == c_top.mon_type[j]) {

		        	rj = d_r[j + traj * c_par.Ntot];
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
						

						if (dr < PAIR_CUTOFF - 0.5){
							if (ind == 0){ 
								c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + i * c_top.maxLateralPerMonomer + c_top.lateralCount[i + traj * c_par.Ntot]] = -j;
								c_top.lateralCount[i + traj * c_par.Ntot]++;
							} else if (c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + i * c_top.maxLateralPerMonomer + c_top.lateralCount[i + traj * c_par.Ntot] - 1] + j != 0){
								c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + i * c_top.maxLateralPerMonomer + c_top.lateralCount[i + traj * c_par.Ntot]] = j;
								c_top.lateralCount[i + traj * c_par.Ntot]++;
							}
		 	
						}
		                /*
		                if (ind == 0) {
		                    if (dr < curMinDistArr[0]) {
		                        curMinDistArr[0] = dr;
		                        c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + i * c_top.maxLateralPerMonomer + 0] = -j;
		                    }
		                } else {
		                    if ((dr < curMinDistArr[1]) && (c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + i * c_top.maxLateralPerMonomer + 0] + j != 0)) {
		                        curMinDistArr[1] = dr;
		                        c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + i * c_top.maxLateralPerMonomer + 1] = j;
		                    }
		                }
		                */

					}
					/*
                    if (curMinDistArr[0] < curMinDistArr[1] && curMinDistArr[0] != PAIR_CUTOFF){
                    	c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + i * c_top.maxLateralPerMonomer + 0] = -j;
                    } else if (curMinDistArr[1] < curMinDistArr[0] && curMinDistArr[1] != PAIR_CUTOFF){
                    	c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + i * c_top.maxLateralPerMonomer + 1] = j;
                    } 

		           */
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
	real psiji, thetaji, fiji, psiij, thetaij, fiij;
	int j;
	Coord ri, rj;
	real dr, dr2;
	Energies en;
	real U_lat = 0.0, U_long = 0.0, U_harm = 0.0, U_fi = 0.0, U_psi = 0.0, U_teta = 0.0, U_lj = 0.0;

	real xp1 = xp1_def;
	real yp1 = yp1_def;
	real zp1 = zp1_def;
	real xp2 = xp2_def;
	real yp2 = yp2_def;
	real zp2 = zp2_def;
	real R_MON = r_mon;

	if (ind < c_par.Ntot && traj < c_par.Ntr){
		if(!c_top.extra[ind + traj * c_par.Ntot]){
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
				j = c_top.harmonic[c_top.maxHarmonicPerMonomer * ind + k];
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

				U_harm += (c_par.C / 2) * pow(dr,2);

	            if(dr < ANGLE_CUTOFF){

	            	psiji = rj.psi - ri.psi - 2 * M_PI * (int)((rj.psi - ri.psi)/(2 * M_PI));
	                psiij = - psiji;
	            	thetaji = rj.theta - ri.theta - 2 * M_PI * (int)((rj.theta - ri.theta)/(2 * M_PI));
	            	thetaij = - thetaji;
	              	fiji = rj.fi - ri.fi - 2 * M_PI * (int)((rj.fi - ri.fi)/(2 * M_PI));
	                fiij = - fiji;

	              	if(R_MON > 0){
	                    U_psi  	 += c_par.B_psi		/2	*	pow(psiji 		- c_par.psi_0		,2);
	                    U_fi	 += c_par.B_fi		/2	*	pow(fiji 		- c_par.fi_0		,2);
	                    U_teta	 += c_par.B_theta	/2	*	pow(thetaji 	- c_par.theta_0		,2);
	                }
	                else{
	                    U_psi  	 += c_par.B_psi	  	/2	*	pow(psiij 		- c_par.psi_0		,2);
	                    U_fi	 += c_par.B_fi	  	/2	*	pow(fiij 		- c_par.fi_0		,2);
	                    U_teta	 += c_par.B_theta	/2	*	pow(thetaij 	- c_par.theta_0		,2);
	                }
	              	/*
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
	                */
	            }
			}

#if defined(MORSE)		
			for(int k = 0; k < c_top.longitudinalCount[ind + traj * c_par.Ntot]; k++){
				j = c_top.longitudinal[c_top.maxLongitudinalPerMonomer * c_par.Ntot * traj + c_top.maxLongitudinalPerMonomer * ind + k];
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


	            U_long += morse_en(c_par.D_long, c_par.A_long, dr);
				//U_long += (c_par.A_long*(c_par.b_long * dr2 * exp(-dr / c_par.r0_long) - c_par.c_long*exp(-dr2/(c_par.d_long*c_par.r0_long)))); 


#if defined(BARR)
	            U_long += barr(c_par.a_barr_long, c_par.r_barr_long, c_par.w_barr_long, dr);
#endif
				if(dr < ANGLE_CUTOFF){

	            	psiji = rj.psi - ri.psi - 2 * M_PI * (int)((rj.psi - ri.psi)/(2 * M_PI));
	                psiij = - psiji;
	            	thetaji = rj.theta - ri.theta - 2 * M_PI * (int)((rj.theta - ri.theta)/(2 * M_PI));
	            	thetaij = - thetaji;
	              	fiji = rj.fi - ri.fi - 2 * M_PI * (int)((rj.fi - ri.fi)/(2 * M_PI));
	                fiij = - fiji;

	              	if(R_MON > 0){
	                    U_psi  	 += c_par.B_psi		/2	*	pow(psiji 		- c_par.psi_0		,2);
	                    U_fi	 += c_par.B_fi		/2	*	pow(fiji 		- c_par.fi_0		,2);
	                    U_teta	 += c_par.B_theta	/2	*	pow(thetaji 	- c_par.theta_0		,2);
	                }
	                else{
	                    U_psi  	 += c_par.B_psi	  	/2	*	pow(psiij 		- c_par.psi_0		,2);
	                    U_fi	 += c_par.B_fi	  	/2	*	pow(fiij 		- c_par.fi_0		,2);
	                    U_teta	 += c_par.B_theta	/2	*	pow(thetaij 	- c_par.theta_0		,2);
	                }
	              	/*
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
	                */
	            }
			}
#endif

#if defined(MORSE)
			for(int k = 0; k < c_top.lateralCount[ind + traj * c_par.Ntot]; k++){
				j = c_top.lateral[c_top.maxLateralPerMonomer * c_par.Ntot * traj + c_top.maxLateralPerMonomer * ind + k];
				
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
			

	            if (c_top.mon_type[ind] != c_top.mon_type[j]) {
					U_lat += morse_en(c_par.D_lat / 2, c_par.A_lat, dr); 	
				}
	            else {
	            	U_lat += morse_en(c_par.D_lat, c_par.A_lat, dr);
	            }
				


#if defined(BARR)
	            U_lat += barr(c_par.a_barr_lat, c_par.r_barr_lat, c_par.w_barr_lat, dr);
#endif
			}
#endif
		
#ifdef LJ_on

	        for(int k = 0; k < c_top.LJCount[ind + traj * c_par.Ntot]; k++){
	        	j = c_top.LJ[c_top.maxLJPerMonomer * c_par.Ntot * traj + ind * c_top.maxLJPerMonomer + k];
				rj = d_r[j + traj * c_par.Ntot];
	            dr = sqrt(pow(ri.x - rj.x, 2) + pow(ri.y - rj.y, 2) + pow(ri.z - rj.z, 2));
	            if(dr < lj_cutoff){
	                U_lj += c_par.ljscale * c_par.ljsigma6 / pow(dr,6);
	            }
	        }
#endif       

		}

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
    Coord ri, rj;
    
    if(i < c_par.Ntot && tr < c_par.Ntr){

    	ri = r[i + tr * c_par.Ntot];
	    c_top.LJCount[i + tr * c_par.Ntot] = 0;

    	if(!c_top.extra[i + tr * c_par.Ntot]){

	        for(int j = 0; j < c_par.Ntot; j++){
	            rj = r[j + tr * c_par.Ntot];
	            real dr = sqrt(pow(ri.x - rj.x,2)+
	                            pow(ri.y - rj.y,2)+
	                            pow(ri.z - rj.z,2));

	            if((dr < c_par.ljpairscutoff) && (i != j)){
	                c_top.LJCount[i + tr * c_par.Ntot]++;
	                c_top.LJ[c_top.maxLJPerMonomer * c_par.Ntot * tr + i * c_top.maxLJPerMonomer + c_top.LJCount[i + tr * c_par.Ntot] - 1] = j;
	            }
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
		if(!c_top.fixed[p % c_par.Ntot] && !(c_top.extra[p])){
			f = d_f[p];
			ri = d_r[p];
			rf_xyz = rforce(p);
			rf_ang = rforce(p + c_par.Ntot*c_par.Ntr);

			ri.x += (c_par.dt/c_par.gammaR)*f.x + c_par.varR*rf_xyz.x;
			ri.y += (c_par.dt/c_par.gammaR)*f.y + c_par.varR*rf_xyz.y;
			ri.z += (c_par.dt/c_par.gammaR)*f.z + c_par.varR*rf_xyz.z;

            // Disallow rotations in all directions but theta
			ri.fi    += (c_par.dt/(c_par.gammaTheta * c_par.alpha))*f.fi  + (c_par.varTheta * sqrt(c_par.freeze_temp / c_par.alpha))*rf_ang.x;
			ri.psi   += (c_par.dt/(c_par.gammaTheta * c_par.alpha))*f.psi + (c_par.varTheta * sqrt(c_par.freeze_temp / c_par.alpha))*rf_ang.y;
			ri.theta += (c_par.dt/c_par.gammaTheta)*f.theta + c_par.varTheta*rf_ang.z;

            
			d_r[p] = ri;
		}
		
		/*
		f.x = 0.0f;
		f.y = 0.0f;
		f.z = 0.0f;
		f.fi = 0.0f;
		f.psi = 0.0f;
		f.theta = 0.0f;
		d_f[p] = f;
		*/
	}
}

void compute(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies){

	Coord* d_r;
	Coord* d_f;
	Topology topGPU;
	
	cudaSetDevice(par.device);
	checkCUDAError("device");

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
	cudaMalloc((void**)&(topGPU.extra), par.Ntr * par.Ntot*sizeof(bool));
	checkCUDAError("topGPU.extra allocation");
	cudaMemcpy(topGPU.fixed, top.fixed, par.Ntot*sizeof(bool), cudaMemcpyHostToDevice);
	checkCUDAError("fixed copy");
	cudaMemcpy(topGPU.mon_type, top.mon_type, par.Ntot*sizeof(int), cudaMemcpyHostToDevice);
	checkCUDAError("montype copy");
	cudaMemcpy(topGPU.extra, top.extra, par.Ntr * par.Ntot*sizeof(bool), cudaMemcpyHostToDevice);
	checkCUDAError("extra copy");

#ifdef LJ_on
	topGPU.maxLJPerMonomer = 256;  ///I hope it will be enough       =top.maxLJPerMonomer;
	cudaMalloc((void**)&(topGPU.LJCount), par.Ntot*sizeof(int)*par.Ntr);
	checkCUDAError("lj_count allocation");
	cudaMalloc((void**)&(topGPU.LJ), par.Ntot*sizeof(int)*par.Ntr*topGPU.maxLJPerMonomer);
	checkCUDAError("lj allocation");	
#endif

	//const memory
	cudaMemcpyToSymbol(c_top, &topGPU, sizeof(Topology), 0, cudaMemcpyHostToDevice);
	checkCUDAError("copy of topGPU pointer to const memory");

	cudaMemcpyToSymbol(c_par, &par, sizeof(Parameters), 0, cudaMemcpyHostToDevice);
	checkCUDAError("copy parameters to const memory");

#ifdef OUTPUT_EN     //energies initializing
	Energies* d_energies;
	cudaMalloc((void**)&d_energies, par.Ntot*par.Ntr * sizeof(Energies));
	checkCUDAError("d_energies allocation");
#endif

	initRand(par.rseed, 2*par.Ntot*par.Ntr);

    int* mt_len = (int*)malloc(par.Ntr * sizeof(int));
	int* mt_len_prev = (int*)malloc(par.Ntr * sizeof(int));
	//int* delta = (int*)malloc(par.Ntr * sizeof(int))

	for(long long int step = 0; step < par.steps; step++){


		if(step % par.ljpairsupdatefreq == 0){ //pairs update frequency 
	#ifdef LJ_on
			//printf("LJPairs are updated");
			LJ_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r);
			checkCUDAError("lj_kernel");
	#endif

	#if defined(ASSEMBLY)

			pairs_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r);
            checkCUDAError("pairs_kernel");

            
/*
            cudaMemcpy(top.lateral, topGPU.lateral, par.Ntot*par.Ntr*topGPU.maxLateralPerMonomer*sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(top.lateralCount, topGPU.lateralCount, par.Ntot*sizeof(int), cudaMemcpyDeviceToHost);
			checkCUDAError("harmonic count copy");
            for (int i = 0; i < par.Ntot; i ++ ){
            	//printf("%d:\t", i);
            	for (int k = 0; k < top.lateralCount[i]; k++){
            		//printf("%d ", top.longitudinal[top.maxLongitudinalPerMonomer * i + k]);
            		if (k > 3) printf("=======================================================");
            	}
            	//printf("\n");
            }
            */

	#endif
		}

		compute_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r, d_f);
		checkCUDAError("compute_kernel");

		integrate_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r, d_f);
		checkCUDAError("integrate_kernel");

		if(step % par.stride == 0){ //every stride steps do energy computing and outputing DCD
			/*	////////////////////*
			cudaMemcpy(top.longitudinal, topGPU.longitudinal, par.Ntot*par.Ntr*topGPU.maxLongitudinalPerMonomer*sizeof(int), cudaMemcpyDeviceToHost);
            cudaMemcpy(top.longitudinalCount, topGPU.longitudinalCount, par.Ntot*sizeof(int), cudaMemcpyDeviceToHost);
			checkCUDAError("harmonic count copy");
            for (int i = 0; i < par.Ntot; i ++ ){
            	if (top.longitudinalCount[i] > 0 && !top.fixed[i]){
            		printf("%d:\t", i);
	            	for (int k = 0; k < top.longitudinalCount[i]; k++){
	            		printf("%d ", top.longitudinal[top.maxLongitudinalPerMonomer * i + k]);
	            		if (k == 2) printf("=======================================================");
	            	}
	            	printf("\n");
            	}
            	
            }//////
*/
	#ifdef OUTPUT_EN
            energy_kernel<<<par.Ntot*par.Ntr/BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r, d_energies);
            checkCUDAError("energy_kernel");
            cudaMemcpy(energies, d_energies, par.Ntr * par.Ntot * sizeof(Energies), cudaMemcpyDeviceToHost);
            checkCUDAError("energy_copy");
	#endif
			cudaMemcpy(r, d_r, par.Ntr*par.Ntot*sizeof(Coord), cudaMemcpyDeviceToHost);
			checkCUDAError("r copy");

	#if defined (OUTPUT_FORCE)
			cudaMemcpy(f, d_f, par.Ntot*par.Ntr*sizeof(Coord), cudaMemcpyDeviceToHost);
	#endif

#ifdef MT_LENGTH
			
			if (step != 0){

				memcpy(mt_len_prev, mt_len, par.Ntr*sizeof(float));
				mt_length(step, mt_len);
	#ifdef CONCENTRATION
				
				for (int i = 0; i < par.Ntr; i++){
					mt_len_prev[i] = mt_len[i] - mt_len_prev[i];
				}

				if (change_conc(mt_len_prev, mt_len)){
					cudaMemcpy(topGPU.extra, top.extra, par.Ntr * par.Ntot*sizeof(bool), cudaMemcpyHostToDevice);
					checkCUDAError("extra copy to device");
					cudaMemcpy(d_r, r, par.Ntot*par.Ntr*sizeof(Coord), cudaMemcpyHostToDevice);
					checkCUDAError("from r to d_r copy");
					cudaMemcpyToSymbol(c_par, &par, sizeof(Parameters), 0, cudaMemcpyHostToDevice);
					checkCUDAError("copy parameters to const memory");

					cudaMemcpyToSymbol(c_top, &topGPU, sizeof(Topology), 0, cudaMemcpyHostToDevice);
					checkCUDAError("copy of topGPU pointer to const memory");
				}
	#endif
				update(step, mt_len);	
				
			} 
			else{
				update(step, mt_len);
				mt_length(step, mt_len);
			}

#endif
		
#ifndef MT_LENGTH			
			update(step, mt_len);
#endif			
		}
	}

	cudaFree(d_r);
	cudaFree(d_f);

	cudaFree(topGPU.harmonicCount);
	cudaFree(topGPU.longitudinalCount);
	cudaFree(topGPU.lateralCount);
	cudaFree(topGPU.fixed);
	cudaFree(topGPU.mon_type);
	cudaFree(topGPU.extra);
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