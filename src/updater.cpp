#include "updater.h"

void OutputAllEnergies(long long int step){

    Energies* fullEnergy = (Energies*)malloc(par.Ntr * sizeof(Energies));

    for (int tr = 0; tr < par.Ntr; tr++){
        fullEnergy[tr].U_harm = 0;
        fullEnergy[tr].U_long = 0;
        fullEnergy[tr].U_lat = 0;
        fullEnergy[tr].U_psi = 0;
        fullEnergy[tr].U_fi = 0;
        fullEnergy[tr].U_teta = 0;
        fullEnergy[tr].U_lj = 0;

        for (int i = 0; i < par.Ntot; i++){
            fullEnergy[tr].U_harm += energies[i + par.Ntot * tr].U_harm;
            fullEnergy[tr].U_long += energies[i + par.Ntot * tr].U_long;
            fullEnergy[tr].U_lat += energies[i + par.Ntot * tr].U_lat;
            fullEnergy[tr].U_psi += energies[i + par.Ntot * tr].U_psi;
            fullEnergy[tr].U_fi += energies[i + par.Ntot * tr].U_fi;
            fullEnergy[tr].U_teta += energies[i + par.Ntot * tr].U_teta;
            fullEnergy[tr].U_lj += energies[i + par.Ntot * tr].U_lj;
            //printf("ENERGY[%d] long = %f\n", i, energies[i + par.Ntot * tr].U_long);
        }

        if (!isfinite(fullEnergy[tr].U_harm) || !isfinite(fullEnergy[tr].U_long) || !isfinite(fullEnergy[tr].U_lat) || !isfinite(fullEnergy[tr].U_lj)) {
            printf("Some energy in %d trajectory is NaN. NOT Exit program\n", tr);
            //exit(0);
        }

        //char fileName[64];
        //sprintf(fileName, "energies%d.dat", tr);
        //FILE* energyFile = fopen(fileName, "a");
        printf("Energies[%d]:\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", tr, 
            fullEnergy[tr].U_harm, fullEnergy[tr].U_long, fullEnergy[tr].U_lat, fullEnergy[tr].U_psi, fullEnergy[tr].U_fi, fullEnergy[tr].U_teta, fullEnergy[tr].U_lj);

        //fprintf(energyFile, "%lld\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\t\t%f\n",
        // step, fullEnergy[tr].U_harm, fullEnergy[tr].U_long, fullEnergy[tr].U_lat, fullEnergy[tr].U_psi, fullEnergy[tr].U_fi, fullEnergy[tr].U_teta, fullEnergy[tr].U_lj);
        //fclose(energyFile);
    }
    
}

void OutputSumForce(){
    for(int tr = 0; tr < par.Ntr; tr++){
        Coord f_sum;
        f_sum.x = 0; f_sum.y = 0; f_sum.z = 0;
        for(int i = tr * par.Ntot; i < (tr + 1) * par.Ntot; i++){
            f_sum.x += f[i].x;
            f_sum.y += f[i].y;
            f_sum.z += f[i].z;
        }
        printf("SF for %d traj %.16f\n", tr, sqrt(f_sum.x * f_sum.x + f_sum.y * f_sum.y + f_sum.z * f_sum.z));
       // printf("SF for %d traj: %f %f %f \n", tr, f_sum.x, f_sum.y, f_sum.z);// * f_sum.x + f_sum.y * f_sum.y + f_sum.z * f_sum.z));
    } 
}

void OutputForces(){

    for (int i = 0; i < par.Ntr; i++){
         printf("force%d\t%f\t%f\t%f\n",i, f[sphere + i * par.Ntot].x, f[sphere + i * par.Ntot].y, f[sphere + i * par.Ntot].z);
    }
   

    //for (int i = 0; i < par.Ntot * par.Ntr; i++){
        //printf("Force[%d].x = %f, y = %f, z = %f\n", i, f[i].x, f[i].y, f[i].z);
        //if (f[i].theta)
        
        //printf("Angle[%d].theta = %f, fi = %f, psi = %f\n", i, r[i].theta, r[i].fi, r[i].psi);
        //printf("Force[%d].x = %f, y = %f, z = %f\n", i, f[i].x, f[i].y, f[i].z);
        //printf("Coord[%d].x = %f, y = %f, z = %f\n", i, r[i].x, r[i].y, r[i].z);
        //printf("type[%d]  %d\n", i, top.mon_type[i]);

    //}
    
}

void update(long long int step, int* mt_len){
    printf("Saving coordinates at step %lld\n", step);
    printTime(step);
    printEstimatedTimeleft((real)step/(real)par.steps);
    saveCoordDCD();
    if (par.hydrolysis && (step % (par.stride * 10) == 0) ){
        if (step == 0){
            FILE* firsttime = fopen("dcd/hydrolysis.pdb", "w");
            fclose(firsttime);
        }
        appendCoordPDB();
    }

    if (par.out_force){
    	//OutputSumForce();
    	OutputForces();
    }
    if (par.out_energy){
    	OutputAllEnergies(step);
    }

}

int change_conc(int* delta, int* mt_len){
    int flag = 0;

    for(int tr = 0; tr < par.Ntr; tr++){

        int num_of_extra = 0;
        for (int i = 0; i < par.Ntot; i++){
            if (top.extra[i + tr * par.Ntot]){
                num_of_extra++;
            } 
        }
        float Vol = float(3.14 * par.rep_r * par.rep_r * par.zs[tr]);
        int NFreeDimers = (par.Ntot - mt_len[tr] - num_of_extra) / 2;

        while (1.0e7 * NFreeDimers / 6.0 < par.conc * Vol){
            for(int i = 0; i < par.Ntot; i+=2){
                if (top.extra[i + tr * par.Ntot] && top.mon_type[i] == 0){
                    top.extra[i + tr * par.Ntot] = false;
                    top.extra[i + tr * par.Ntot+1] = false;
                    float x,y;
                    for(int j = 0; j > -1 ; j++) {
                        x = par.rep_r - 2*(rand() % int(par.rep_r));
                        y = par.rep_r - 2*(rand() % int(par.rep_r));
                        if (x*x + y*y <= par.rep_r * par.rep_r){   //// check if its in cylinder
                            printf("New x,y coordinates for extra particle: %f  %f index: %d\n", x, y, i + tr * par.Ntot);
                            num_of_extra -= 2;
                            break;
                        }
                    }
                    float z = par.zs[tr] + par.rep_leftborder + 3*2*r_mon;
                    r[i + tr * par.Ntot].x = x; //
                    r[i + tr * par.Ntot].y = y;
                    r[i + tr * par.Ntot].z = z; ///fix
                    r[i + tr * par.Ntot +1].x = x;
                    r[i + tr * par.Ntot +1].y = y;
                    r[i + tr * par.Ntot +1].z = z + 2*r_mon;
                    flag++;
                    break;
                }
            }

            NFreeDimers += 1;
            if (num_of_extra == 0) {
                printf("No more extra particles for trajectory[%d]!\n",tr); 
                break;   
            }
        }
            
           // par.zs[tr] += 2.0 * r_mon * delt / 13.0;
               
        //Vol = float(3.14 * par.rep_r * par.rep_r * par.zs[tr]);
        //Nfree = par.Ntot - mt_len[tr] - num_of_extra;
        printf("Concentration for tajectory[%d]: %f [muMole / L],\t %f [Dimers / nm^3],\t %d [Dimers / Volume],\t  Volume: %f [nm^3]\n", tr, 1.0e7 * NFreeDimers / (6.0 * Vol), NFreeDimers / Vol, NFreeDimers, Vol);
    }
    return flag;
}

void mt_length(long long int step, int* mt_len){
    if (step == 0){
        FILE* first = fopen("mt_len.dat", "w");
        fclose(first);
    }

    for (int traj = 0; traj < par.Ntr; traj++){
        int sum = 0;

        //int protofilamentsLength[PF_NUMBER] = {0,0,0,0,0,0,0,0,0,0,0,0,0};
        for (int i = traj * par.Ntot; i < (traj + 1) * par.Ntot; i++) {

            float rad = sqrt(r[i].x * r[i].x + r[i].y * r[i].y);
            if ((rad < R_MT + R_THRES) && (rad > R_MT - R_THRES) && (cosf(r[i].theta) > cosf(ANG_THRES))) {
               
            // detect particles inside microtubule 
                sum++;
                top.on_tubule[i] = 1;

            /*  
                for (int j = 0; j < PF_NUMBER; j++){

                    float angle = j * 2 * M_PI / PF_NUMBER;
                    float psi = 0.0;

                    if (r[i].psi >= -ANG_THRES) {
                        psi = r[i].psi;
                    }
                    else {
                        psi = 2 * M_PI + r[i].psi;
                    }

                    if (psi < angle + ANG_THRES && psi > angle - ANG_THRES ) {

                        if ( (r[i].theta < ANG_THRES && r[i].theta > -ANG_THRES) 
                           ||    (( 2 * M_PI - fabs(r[i].theta)) < ANG_THRES &&  (2 * M_PI - fabs(r[i].theta)) > -ANG_THRES) ) {
                            if (r[i].fi < ANG_THRES && r[i].fi > - ANG_THRES){
                                protofilamentsLength[j]++;
                                //count++;
                                //printf("%d\t%f\t%f\t%f\t\t%f\t%f\t%f\t\t%d%s\t%d\t%d\n", j, r[i].x, r[i].y,r[i].z,r[i].fi,r[i].psi,r[i].theta, pdb.atoms[i].resid, pdb.atoms[i].name, count, count2);
                                break;
                            }
                            
                        }
                    }

                }
                */
            } else {
                top.on_tubule[i] = 0;
            }
        }
        /*
        int sum = 0;
        for (int i = 0; i < PF_NUMBER; i++){
            //if (protofilamentsLength[i] % 2) protofilamentsLength[i]--;
            //printf("PF[%d] = %d\n", i, protofilamentsLength[i]);
            sum += protofilamentsLength[i];

        }*/
        mt_len[traj] = sum;
        printf("tubule[%d]: %d\n", traj, sum);
        /*char fileName[64];
        sprintf(fileName, "mtlength%d.dat", traj);
        */
    }
    FILE* mtLenFile = fopen("mt_len.dat", "a");
    fprintf(mtLenFile, "%lld\t" , step);
    for (int traj = 0; traj < par.Ntr; traj++){
        fprintf(mtLenFile, "%f\t" , 4.0 * (float)mt_len[traj] / PF_NUMBER);
    }
    fprintf(mtLenFile, "\n");
    fclose(mtLenFile);
        
}

void hydrolyse(){
    
    double prob;
    for(int i = 0; i < par.Ntot; i+=2){
        for (int tr = 0; tr < par.Ntr; tr++){
            if (top.gtp[i + tr * par.Ntot] == 1 && !top.extra[i + tr * par.Ntot]){
                prob = rand() / (double)RAND_MAX;
                if (prob < 0.02){
                    top.gtp[i     + tr * par.Ntot] = 0;
                    top.gtp[i + 1 + tr * par.Ntot] = 0;

                    printf("*** Hydrolysis occured to dimer # %d trajectory #%d ***\n", i/2, tr);
                }
            }
        }     
    }
}