#include "mt.h"
#include "output.h"
#include "parameters.h"
extern void compute(Coord* r, Coord* f, Parameters &par, Topology &top, Energies* energies);

double barr(double a, double r, double w, double x)
{
    return a*expf(-(x-r)*(x-r)/(2*w*w));
}

double exp_en(double g0, double k, double x)
{
    return -g0*exp(-k*x);
}

double morse_en(double D, double A, double x)
{
    return D*(1-exp(-A*x))*(1-exp(-A*x)) - D;
}

double U(Coord* r, Parameters &par, Topology &top){
	int ind, traj;
	double cos_fii, cos_fij, sin_fii, sin_fij, 
		  cos_psii, cos_psij, sin_psii, sin_psij,
		  cos_thetai, cos_thetaj, sin_thetai, sin_thetaj;
	double xi, xj, yi, yj, zi, zj;
	double dUdr, dr, dr2, gradx, grady, gradz, gradfi, gradpsi, gradtheta, expon, expon2, df;
	int i,j;
	Coord ri, rj;
	double xp1 = xp1_def;
	double yp1 = yp1_def;
	double zp1 = zp1_def;
	double xp2 = xp2_def;
	double yp2 = yp2_def;
	double zp2 = zp2_def;
	double R_MON = r_mon;
	double U = 0;
    double U_lat = 0, U_long = 0, U_harm = 0, U_fi = 0, U_psi = 0, U_teta = 0, U_lj = 0;

	for(int p = 0; p < par.Ntot; p++){
		ind = p%par.Ntot;
		traj = p/par.Ntot;
		if(ind < par.Ntot && traj < par.Ntr){
			ri = r[p];
			cos_fii = cosf(ri.fi); 
			sin_fii = sinf(ri.fi);
			cos_psii = cosf(ri.psi);
			sin_psii = sinf(ri.psi);
			cos_thetai = cosf(ri.theta);
			sin_thetai = sinf(ri.theta);
			xi = ri.x;
			yi = ri.y;
			zi = ri.z;
			for(int k = 0; k < top.harmonicCount[ind]; k++){
				j = top.harmonic[top.maxHarmonicPerMonomer*ind+k];
				if(j<0){
					R_MON = r_mon;
					j *= -1;
				}
				else{
					R_MON = -r_mon;
				}
				rj = r[j + traj*par.Ntot];
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
                U_harm += par.A_xx+par.b_xx/(par.A_xx+exp(-par.c_xx*(dr-par.d_xx)));
#else
				U_harm += (par.C/2)*pow(dr,2);
#endif
                if(dr < ANGLE_CUTOFF)
                {
                    if(R_MON > 0){
                        U_psi  	 += par.B_psi		/2	*	pow(rj.psi 		-	ri.psi 		- par.psi_0		,2);
                        U_fi	 += par.B_fi		/2	*	pow(rj.fi 		-	ri.fi 		- par.fi_0		,2);
                        U_teta	 += par.B_theta		/2	*	pow(rj.theta 	- 	ri.theta 	- par.theta_0	,2);
                    }
                    else{
                        U_psi  	 += par.B_psi	  	/2	*	pow(ri.psi 		- 	rj.psi 		- par.psi_0		,2);
                        U_fi	 += par.B_fi	  	/2	*	pow(ri.fi 		- 	rj.fi 		- par.fi_0		,2);
                        U_teta	 += par.B_theta		/2	*	pow(ri.theta 	- 	rj.theta 	- par.theta_0	,2);
                    }
                }
			}
			for(int k = 0; k < top.longitudinalCount[ind]; k++){
				j = top.longitudinal[top.maxLongitudinalPerMonomer*ind+k];
				if(j<0){
					R_MON = r_mon;
					j *= -1;
				}
				else{
					R_MON = -r_mon;
				}
				rj = r[j + traj*par.Ntot];
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
#if defined(SIGMOID)
                U_long += par.A_xy+par.b_xy/(par.A_xy+exp(-par.c_xy*(dr-par.d_xy)));
#elif defined(EXP)
                U_long += exp_en(par.g0_long, par.k_long, dr);
#elif defined(MORSE)
                U_long += morse_en(par.D_long, par.A_long, dr);
#else
				U_long += (par.A_long*(par.b_long*dr2*exp(-dr/par.r0_long) - par.c_long*exp(-dr2/(par.d_long*par.r0_long)))); 
#endif
#if defined(BARR)
                U_long += barr(par.a_barr_long, par.r_barr_long, par.w_barr_long, dr);
#endif
				if(dr < ANGLE_CUTOFF)
                {
                    if(R_MON > 0){
                        U_psi  	 += par.B_psi		/2	*	pow(rj.psi 		-	ri.psi 		- par.psi_0		,2);
                        U_fi	 += par.B_fi		/2	*	pow(rj.fi 		-	ri.fi 		- par.fi_0		,2);
                        U_teta	 += par.B_theta		/2	*	pow(rj.theta 	- 	ri.theta 	- par.theta_0	,2);
                    }
                    else{
                        U_psi  	 += par.B_psi	  	/2	*	pow(ri.psi 		- 	rj.psi 		- par.psi_0		,2);
                        U_fi	 += par.B_fi	  	/2	*	pow(ri.fi 		- 	rj.fi 		- par.fi_0		,2);
                        U_teta	 += par.B_theta		/2	*	pow(ri.theta 	- 	rj.theta 	- par.theta_0	,2);
                    }
                }
			}
			for(int k = 0; k < top.lateralCount[ind]; k++){
				j = top.lateral[top.maxLateralPerMonomer*ind+k];
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

				rj = r[j + traj*par.Ntot];
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
#if defined(SIGMOID)
                float c_lat, d_lat, b_lat, A_lat;
                int type = (top.mon_type[p]==0 && top.mon_type[j]==0) ? 0: //AA
                            ((top.mon_type[p]==0 && top.mon_type[j]==1)? 1 : //AB
                            ((top.mon_type[p]==1 && top.mon_type[j]==0)? 2 : 3)); //BA:BB
                A_lat = par.A_lat[type];
                b_lat = par.b_lat[type];
                c_lat = par.c_lat[type];
                d_lat = par.d_lat[type];
                U_lat += A_lat+b_lat/(A_lat+exp(-c_lat*(dr-d_lat)));
#elif defined(EXP)
                U_lat += exp_en(par.g0_lat, par.k_lat, dr);
#elif defined(MORSE)
                U_lat += morse_en(par.D_lat, par.A_lat, dr);
#else
				U_lat += (par.A_lat*(par.b_lat*dr2*exp(-dr/par.r0_lat) - par.c_lat*exp(-dr2/(par.d_lat*par.r0_lat))));
#endif
#if defined(BARR)
                U_lat += barr(par.a_barr_lat, par.r_barr_lat, par.w_barr_lat, dr);
#endif
			}
#ifdef LJ_on
            for(int k=0; k < top.LJCount[p]; k++)
            {
                rj = r[top.LJ[k*par.Ntr*par.Ntot+p]];
                dr = sqrt(pow(ri.x-rj.x,2)+pow(ri.y-rj.y,2)+pow(ri.z-rj.z,2));
                if(dr < lj_cutoff){
                    U_lj += par.ljscale*par.ljsigma6/pow(dr,6);
                }
            }
#endif
		}
	}
    U = U_harm + U_long + U_lat + U_lj + U_psi + U_fi + U_teta;

    Energies en;
    en.U_harm = U_harm/2;
    en.U_long = U_long/2;
    en.U_lat = U_lat/2;
    en.U_psi = U_psi/2;
    en.U_fi = U_fi/2;
    en.U_teta = U_teta/2;
    en.U_lj = U_lj/2;
    OutputEn(en);


    for (int p = 0; p < par.Ntot; p++){
    	if (!(top.lateral[2*p] == top.lateral[2*p + 1] && top.lateral[2*p] == 999999))
    		printf("Lateral[%d] -  %d %d\n", p, top.lateral[2*p], top.lateral[2*p + 1]);
    }
    for (int p = 0; p < par.Ntot; p++){
    	if (top.longitudinalCount[p] != 0)
    		printf("Longitudial[%d] -  %d, count = %d\n", p, top.longitudinal[p], top.longitudinalCount[p]);
    }
	return U;
}



void num_test(Coord* r, Coord* f, Parameters &par, Topology &top){
	Coord* r_old = (Coord*)calloc(par.Ntot, sizeof(Coord));
	Coord* f_num = (Coord*)calloc(par.Ntot, sizeof(Coord));
	Coord* f_sim = (Coord*)calloc(par.Ntot, sizeof(Coord));
	real U0 = U(r, par, top);
	real U1=0, U2=0;
	Energies* energies = NULL;
	for(int i = 0; i < par.Ntot; i++){
		f_num[i].x = 0;
		f_num[i].y = 0;
		f_num[i].z = 0;
		f_num[i].fi = 0;
		f_num[i].psi = 0;
		f_num[i].theta = 0;
	}
	for(int i = 0; i < par.Ntot; i++){
		r_old[i].x = r[i].x;
		r_old[i].y = r[i].y;
		r_old[i].z = r[i].z;
		r_old[i].fi = r[i].fi;
		r_old[i].psi = r[i].psi;
		r_old[i].theta = r[i].theta;
	}
	for(int i = 0; i < par.Ntot; i++){	
		r[i].x = r_old[i].x + delta_x;
		U1 = U(r, par, top);
		r[i].x = r_old[i].x - delta_x;
		U2 = U(r, par, top);
		f_num[i].x = (U2-U1)/(2*delta_x)/2;
		r[i].x = r_old[i].x;
		compute(r, f, par, top, energies);
		f_sim[i].x = f[i].x;
		U1=0; U2=0;
		for(int j = 0; j < par.Ntot; j++){
			r[j].x = r_old[j].x;		f[j].x=0;
			r[j].y = r_old[j].y;		f[j].y=0;
			r[j].z = r_old[j].z;		f[j].z=0;
			r[j].fi = r_old[j].fi;		f[j].fi=0;
			r[j].psi = r_old[j].psi;	f[j].psi=0;
			r[j].theta = r_old[j].theta;	f[j].theta=0;
		}

		r[i].y = r_old[i].y + delta_y;
		U1 = U(r, par, top);
		r[i].y = r_old[i].y - delta_y;
		U2 = U(r, par, top);
		f_num[i].y = (U2-U1)/(2*delta_y)/2;
		r[i].y = r_old[i].y;
		compute(r, f, par, top, energies);
		f_sim[i].y = f[i].y;
		U1=0; U2=0;
		for(int j = 0; j < par.Ntot; j++){
			r[j].x = r_old[j].x;		f[j].x=0;
			r[j].y = r_old[j].y;		f[j].y=0;
			r[j].z = r_old[j].z;		f[j].z=0;
			r[j].fi = r_old[j].fi;		f[j].fi=0;
			r[j].psi = r_old[j].psi;	f[j].psi=0;
			r[j].theta = r_old[j].theta;	f[j].theta=0;
		}

		r[i].z = r_old[i].z + delta_z;
		U1 = U(r, par, top);
		r[i].z = r_old[i].z - delta_z;
		U2 = U(r, par, top);
		f_num[i].z = (U2-U1)/(2*delta_z)/2;
		r[i].z = r_old[i].z;
		compute(r, f, par, top, energies);
		f_sim[i].z = f[i].z;
		U1=0; U2=0;
		for(int j = 0; j < par.Ntot; j++){
			r[j].x = r_old[j].x;		f[j].x=0;
			r[j].y = r_old[j].y;		f[j].y=0;
			r[j].z = r_old[j].z;		f[j].z=0;
			r[j].fi = r_old[j].fi;		f[j].fi=0;
			r[j].psi = r_old[j].psi;	f[j].psi=0;
			r[j].theta = r_old[j].theta;	f[j].theta=0;
		}

		r[i].fi = r_old[i].fi + delta_fi;
		U1 = U(r, par, top);
		r[i].fi = r_old[i].fi - delta_fi;
		U2 = U(r, par, top);
		f_num[i].fi = (U2-U1)/(2*delta_fi)/2;
		r[i].fi = r_old[i].fi;
		compute(r, f, par, top, energies);
		f_sim[i].fi = f[i].fi;
		U1=0; U2=0;
		for(int j = 0; j < par.Ntot; j++){
			r[j].x = r_old[j].x;		f[j].x=0;
			r[j].y = r_old[j].y;		f[j].y=0;
			r[j].z = r_old[j].z;		f[j].z=0;
			r[j].fi = r_old[j].fi;		f[j].fi=0;
			r[j].psi = r_old[j].psi;	f[j].psi=0;
			r[j].theta = r_old[j].theta;	f[j].theta=0;
		}

		r[i].psi = r_old[i].psi + delta_psi;
		U1 = U(r, par, top);
		r[i].psi = r_old[i].psi - delta_psi;
		U2 = U(r, par, top);
		f_num[i].psi = (U2-U1)/(2*delta_psi)/2;
		r[i].psi = r_old[i].psi;
		compute(r, f, par, top, energies);
		f_sim[i].psi = f[i].psi;
		U1=0; U2=0;
		for(int j = 0; j < par.Ntot; j++){
			r[j].x = r_old[j].x;		f[j].x=0;
			r[j].y = r_old[j].y;		f[j].y=0;
			r[j].z = r_old[j].z;		f[j].z=0;
			r[j].fi = r_old[j].fi;		f[j].fi=0;
			r[j].psi = r_old[j].psi;	f[j].psi=0;
			r[j].theta = r_old[j].theta;	f[j].theta=0;
		}

		r[i].theta = r_old[i].theta + delta_theta;
		U1 = U(r, par, top);
		r[i].theta = r_old[i].theta - delta_theta;
		U2 = U(r, par, top);
		f_num[i].theta = (U2-U1)/(2*delta_theta)/2;
		r[i].theta = r_old[i].theta;
		compute(r, f, par, top, energies);
		f_sim[i].theta = f[i].theta;
		U1=0; U2=0;
		for(int j = 0; j < par.Ntot; j++){
			r[j].x = r_old[j].x;		f[j].x=0;
			r[j].y = r_old[j].y;		f[j].y=0;
			r[j].z = r_old[j].z;		f[j].z=0;
			r[j].fi = r_old[j].fi;		f[j].fi=0;
			r[j].psi = r_old[j].psi;	f[j].psi=0;
			r[j].theta = r_old[j].theta;	f[j].theta=0;
		}
//		printf("%d %d %d\n",i, top.lateralCount[i], top.lateral[i]);
	}
	for(int k = 0; k < par.Ntot; k++){
		if(sqrt(pow(f_num[k].x-f_sim[k].x,2) + pow(f_num[k].y-f_sim[k].y,2) + pow(f_num[k].z-f_sim[k].z,2) + pow(f_num[k].fi-f_sim[k].fi,2) + pow(f_num[k].psi-f_sim[k].psi,2) + pow(f_num[k].theta-f_sim[k].theta,2)) > 0.1){
		//	printf("%f %f : %f %f : %f %f : %f %f : %f %f : %f %f\n", f_num[k].x, f_sim[k].x, f_num[k].y, f_sim[k].y, f_num[k].z , f_sim[k].z, f_num[k].fi, f_sim[k].fi, f_num[k].psi, f_sim[k].psi, f_num[k].theta , f_sim[k].theta);
			printf("%-*d NUM %*f%*f%*f%*f%*f%*f\n", 16, k, 16, f_num[k].x, 16, f_num[k].y, 16, f_num[k].z, 16, f_num[k].fi, 16, f_num[k].psi, 16, f_num[k].theta);
			printf("%-*d SIM %*f%*f%*f%*f%*f%*f\n", 16, k, 16, f_sim[k].x, 16, f_sim[k].y, 16, f_sim[k].z, 16, f_sim[k].fi, 16, f_sim[k].psi, 16, f_sim[k].theta);
		}
	}
	return;
}

