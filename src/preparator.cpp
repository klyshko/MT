#include "preparator.h"


void initParameters(int argc, char* argv[]){
#ifdef CUDA
#ifdef USE_MPI
    if (getYesNoParameter(PARAMETER_MPI_DEV_AUTO, 0)) {
        int mpi_dpn = getIntegerParameter(PARAMETER_MPI_DEVPERNODE);
        mpi_dev_cur = mpi_rank % mpi_dpn;
    } else if (mpi_size > 1) {
        mpi_dev_cur = getIntegerParameter(PARAMETER_MPI_DEVICE(mpi_rank));
    } else {
        mpi_dev_cur = getIntegerParameter(PARAMETER_MPI_DEVICE(mpi_rank),getIntegerParameter(PARAMETER_DEVICE));
    }
    par.device = mpi_dev_cur;
#else
    par.device = getIntegerParameter(PARAMETER_DEVICE);
#endif
#endif
    par.dt = getFloatParameter(PARAMETER_TIMESTEP);
    par.rseed = getIntegerParameter(PARAMETER_RANDOM_SEED);
    par.steps = getLongIntegerParameter(PARAMETER_NUMSTEPS,-1);
    par.stride = getLongIntegerParameter(PARAMETER_STRIDE,-1);
    par.ljpairscutoff = getFloatParameter(PARAMETER_LJPAIRSCUTOFF);
    par.ljpairsupdatefreq = getIntegerParameter(PARAMETER_LJPAIRSUPDATEFREQ);

    getMaskedParameter(par.coordFilename_ang, PARAMETER_COORD_FILE_ANG);
    getMaskedParameter(par.coordFilename_xyz, PARAMETER_COORD_FILE_XYZ);
    getMaskedParameter(par.ffFilename, PARAMETER_FORCEFIELD_FILE);
    getMaskedParameter(par.condFilename, PARAMETER_CONDITIONS_FILE);
    par.fix = getIntegerParameter(PARAMETER_FIX,1);
    par.Ntr = getIntegerParameter(PARAMETER_RUNNUM,1);

    if (getYesNoParameter(PAR_ASSEMBLY,1)){
    	par.is_assembly = true;
    } else {
    	par.is_assembly = false;
    }
    if (getYesNoParameter(TUB_LENGTH,1)){
    	par.tub_length = true;
    } else {
    	par.tub_length = false;
    }
    if (getYesNoParameter(PAR_OUTPUT_ENERGY,1)){
    	par.out_energy = true;
    } else {
    	par.out_energy = false;
    }
    if (getYesNoParameter(PAR_OUTPUT_FORCE,0)){
    	par.out_force = true;
    } else {
    	par.out_force = false;
    }

    if(getYesNoParameter(BDHITEA_ON_STRING,0))
        par.hdi_on = true;
    else
        par.hdi_on = false;

    if (par.hdi_on) {
        tea.a = getFloatParameter(BDHITEA_A_STRING, r_mon, 1); // in Angstroms
        tea.capricious = getYesNoParameter(BDHITEA_CAPRICIOUS_STRING, 1, 1); // Be paranoid about tensor values?
        tea.epsilon_freq = getIntegerParameter(BDHITEA_EPSILONFREQ_STRING, 0, 0); // How often recalculate epsilon?
        tea.epsmax = getFloatParameter(BDHITEA_EPSMAX_STRING, 999.f, 1); // Epsilon will never exceed 1, so epsmax=999 will never trigger halt by itself; used in capricious mode
        tea.unlisted = 1; 
    } 
    

  
    getMaskedParameter(par.restartkey, PARAMETER_RESTARTKEY);

    if(getYesNoParameter(PARAMETER_ISRESTART,0))
        par.is_restart = true;
    else
        par.is_restart = false;
    if(!par.is_restart)
    {
        par.firststep = getLongIntegerParameter(PARAMETER_FIRSTSTEP, 0);
    }
    else
    {
        FILE *keyf = safe_fopen(par.restartkey, "r");
        if (fscanf(keyf,"%lld",&(par.firststep)) != 1)
            DIE("Reading restartkey %s: unable to get firststep", &(par.restartkey));
        fclose(keyf);
    }   
    read_PDB(par.coordFilename_xyz, par.coordFilename_ang);   ////// <----- read topology from pdb files

    par.dcdFilename_xyz = (char**)calloc(par.Ntr, sizeof(char*));
    par.dcdFilename_ang = (char**)calloc(par.Ntr, sizeof(char*));
    par.restart_xyzFilename = (char**)calloc(par.Ntr, sizeof(char*));
    par.restart_angFilename = (char**)calloc(par.Ntr, sizeof(char*));

#if not defined(READ_FROM_DCD)
    createDCD(&dcd, par.Ntot, par.steps/par.stride, 1, par.dt, par.stride, 0, 0, 0, 0);
    char trajnum[10];
    for(int traj = 0; traj < par.Ntr; traj++){
            sprintf(trajnum, "%d", traj);
            par.dcdFilename_xyz[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
            par.dcdFilename_ang[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
            getMaskedParameterWithReplacement(par.dcdFilename_ang[traj], PARAMETER_DCD_FILE_ANG, trajnum, "<run>");
            getMaskedParameterWithReplacement(par.dcdFilename_xyz[traj], PARAMETER_DCD_FILE_XYZ, trajnum, "<run>");
            dcdOpenWrite(&dcd, par.dcdFilename_xyz[traj]);
            dcdWriteHeader(dcd);
            dcdClose(dcd);
            dcdOpenWrite(&dcd, par.dcdFilename_ang[traj]);
            dcdWriteHeader(dcd);
            dcdClose(dcd);

            par.restart_xyzFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
            par.restart_angFilename[traj] = (char*)calloc(PARAMETER_LENGTH, sizeof(char));
            getMaskedParameterWithReplacement(par.restart_xyzFilename[traj], PARAMETER_RESTART_XYZ_FILE, trajnum, "<run>");
            getMaskedParameterWithReplacement(par.restart_angFilename[traj], PARAMETER_RESTART_ANG_FILE, trajnum, "<run>");

    }
#endif
    if(par.is_restart)
    {
        readRestart();
    }
    parseParametersFile(par.ffFilename, argc, argv);

#if defined(MORSE)
    par.A_long = getFloatParameter(A_LONG);
    par.A_lat  = getFloatParameter(A_LAT);
    par.D_long  = getFloatParameter(D_LONG);
    par.D_lat   = getFloatParameter(D_LAT);
#else
   
    printf("No long/lat potentials. Continue to compute\n");
    //exit(0);
#endif
    
#if defined(BARR)
    par.a_barr_long = getFloatParameter(A_BARR_LONG);
    par.r_barr_long = getFloatParameter(R_BARR_LONG);
    par.w_barr_long = getFloatParameter(W_BARR_LONG)/( 2*(sqrt(-2*log(0.5))) );
    par.a_barr_lat = getFloatParameter(A_BARR_LAT);
    par.r_barr_lat = getFloatParameter(R_BARR_LAT);
    par.w_barr_lat = getFloatParameter(W_BARR_LAT)/( 2*(sqrt(-2*log(0.5))) );
#endif

 //ANGLES_HARMONIC  
    par.C = getFloatParameter(PARAMETER_HARMONIC_C);

    par.B_fi = getFloatParameter(PARAMETER_BENDING_B_FI);
    par.B_psi = getFloatParameter(PARAMETER_BENDING_B_PSI);
    par.B_theta = getFloatParameter(PARAMETER_BENDING_B_THETA);
    par.fi_0 = getFloatParameter(PARAMETER_BENDING_FI_0);
    par.psi_0 = getFloatParameter(PARAMETER_BENDING_PSI_0);
    par.theta0_gdp = getFloatParameter(PARAMETER_BENDING_THETA0_GDP);
    par.theta0_gtp = getFloatParameter(PARAMETER_BENDING_THETA0_GTP);

    if (getYesNoParameter(PARAMETER_LJ_ON, 1, 0)){
    	par.lj_on = true;
    	par.ljscale = getFloatParameter(PARAMETER_LJSCALE);
    	par.ljsigma6 = pow(getFloatParameter(PARAMETER_LJSIGMA),6);
    } else {
    	par.lj_on = false;
    }	

    if (getYesNoParameter(REPULSIVE_WALLS,1,0)){
    	par.is_wall = true;	
    } else {
    	par.is_wall = false;
    }
    
    par.rep_leftborder = getFloatParameter(REP_LEFTBORDER);
    par.rep_h = getFloatParameter(REP_H);
    par.rep_r = getFloatParameter(REP_R);
    par.rep_eps = getFloatParameter(REP_EPS); 
    for (int i = 0; i < par.Ntr; i++){
        par.zs[i] = par.rep_h;
    }   

    parseParametersFile(par.condFilename, argc, argv);

    if (getYesNoParameter(CONST_CONC, 1, 0)){
    	par.is_const_conc = true;
    	par.conc = getFloatParameter(CONC);	
    } else {
    	par.is_const_conc = false;
    }

    if (getYesNoParameter(HYDROLYSIS, 1, 0)){
        par.hydrolysis = true;
        par.khydro = getFloatParameter(KHYDRO);
        par.hydrostep = (long int) (0.02 * 1000000000000 / (par.dt * par.khydro));
    } else {
        par.hydrolysis = false;
    }

    for(int traj = 0; traj < par.Ntr; traj++){
        if (par.hydrolysis){
           for(int i = 0; i < par.Ntot; i++){
                top.gtp[i + traj*par.Ntot] = 1;     
            } 
        } else {
            for(int i = 0; i < par.Ntot; i++){
                top.gtp[i + traj*par.Ntot] = 0;     
            }
        }   
    }

    par.Temp = getFloatParameter(PARAMETER_TEMPERATURE);
    par.viscosity = getFloatParameter(PARAMETER_VISCOSITY);
    if (par.hdi_on) {
        par.gammaR = 6 * M_PI * par.viscosity * tea.a;
    } else {
        par.gammaR = 6 * M_PI * par.viscosity * r_mon;
    }
    par.gammaTheta = 8 * M_PI * par.viscosity * pow(r_mon, 3); 
    
    par.varR = sqrtf(2.0f*KB*par.Temp*par.dt/par.gammaR);
    par.varTheta = sqrtf(2.0f*KB*par.Temp*par.dt/par.gammaTheta);
    par.alpha = getFloatParameter(ALPHA_GEOMETRY);
    par.freeze_temp = getFloatParameter(FREEZE_TEMP);

    int extra_counter = 0;
    for (int i = 0; i < par.Ntot; i++ ){
    	if (top.extra[i]) {
    		extra_counter++;
    	}
    }
    if(extra_counter == 0 && par.is_const_conc) {
    	printf("Error! You want constant concentration! Load structure with extra particles (chain X) to support constant concentration!\n");
    	exit(-1);
    } 

    if (par.hydrolysis) {
        coordspdb.atomCount = par.Ntot * par.Ntr;
        coordspdb.atoms = (PDBAtom*)malloc(coordspdb.atomCount * sizeof(PDBAtom));
        for (int tr = 0; tr < par.Ntr; tr++ ){
            for (int i = 0; i < par.Ntot; i++){
                memcpy(&(coordspdb.atoms[i + tr * par.Ntot]), &(pdb.atoms[i]), sizeof(PDBAtom));
            }
        }

    }
    
    
}


void read_PDB(const char* filename_xyz, const char* filename_ang){
    int i, j;
    readPDB(filename_xyz, &pdb);
    printf("Building topology....\n");

    par.Ntot = pdb.atomCount;

    r = (Coord*)calloc(par.Ntot*par.Ntr, sizeof(Coord));
    f = (Coord*)calloc(par.Ntot*par.Ntr, sizeof(Coord));

    top.mon_type = (int*)calloc(par.Ntot, sizeof(int));
    for(i = 0; i < par.Ntot; i++){
        if(pdb.atoms[i].name[1] == 'A')
            top.mon_type[i] = 0;
        else if(pdb.atoms[i].name[1] == 'B')
            top.mon_type[i] = 1;
    }

    top.gtp = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));

    top.on_tubule = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));

    for(int traj = 0; traj < par.Ntr; traj++){
        for(i = 0; i < par.Ntot; i++){
            top.on_tubule[i + traj*par.Ntot] = 0;     
        } 
    } 

     //initialising fixed atoms list
    top.fixed = (bool*)calloc(par.Ntot, sizeof(bool));
    for(i = 0; i < par.Ntot; i++){
            if( pdb.atoms[i].resid <= par.fix)
                    top.fixed[i] = true;
            else
                    top.fixed[i] = false;
    }

     //initialising extra atoms list
    top.extra = (bool*)calloc(par.Ntot * par.Ntr, sizeof(bool));
   
    for(int traj = 0; traj < par.Ntr; traj++){
        for(i = 0; i < par.Ntot; i++){
            if(pdb.atoms[i].chain == 'X'){
            	top.extra[i + traj*par.Ntot] = true;
            }        
            else {
            	 top.extra[i + traj*par.Ntot] = false;
            }
                   
        }
    }

    for(int traj = 0; traj < par.Ntr; traj++){
            for(i = 0; i < par.Ntot; i++){
                r[i+traj*par.Ntot].x = pdb.atoms[i].x;
                r[i+traj*par.Ntot].y = pdb.atoms[i].y;
                r[i+traj*par.Ntot].z = pdb.atoms[i].z;
            }
    }
    readPDB(filename_ang, &pdb_ang);
    for(int traj = 0; traj < par.Ntr; traj++){
            for(i = 0; i < par.Ntot; i++){
                r[i+traj*par.Ntot].fi = pdb_ang.atoms[i].x;
                r[i+traj*par.Ntot].psi = pdb_ang.atoms[i].y;
                r[i+traj*par.Ntot].theta = pdb_ang.atoms[i].z;
            }
    }


    top.LJCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
    top.harmonicCount = (int*)calloc(par.Ntot, sizeof(int));
    top.longitudinalCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
    top.lateralCount = (int*)calloc(par.Ntot*par.Ntr, sizeof(int));
    
    // Covalent bonds map
    for(i = 0; i < par.Ntot; i++) top.harmonicCount[i] = 0;
    for(i = 0; i < par.Ntot; i++){
        for(j = 0; j < par.Ntot; j++){
            if(pdb.atoms[i].resid == pdb.atoms[j].resid &&
                    pdb.atoms[i].chain == pdb.atoms[j].chain
                    && (i!=j) && (abs(i-j) < 2) ){
                top.harmonicCount[i]++;
                if(top.maxHarmonicPerMonomer < top.harmonicCount[i])
                        top.maxHarmonicPerMonomer = top.harmonicCount[i];
            }
        }
    }
    top.harmonic = (int*)calloc(top.maxHarmonicPerMonomer*par.Ntot, sizeof(int));
    for(i = 0; i < par.Ntot; i++) top.harmonicCount[i] = 0;
    for(i = 0; i < top.maxHarmonicPerMonomer*par.Ntot; i++) top.harmonic[i] = -1;
    for(i = 0; i < par.Ntot; i++){
        for(j = 0; j < par.Ntot; j++){
            if(pdb.atoms[i].resid == pdb.atoms[j].resid &&
                    pdb.atoms[i].chain == pdb.atoms[j].chain &&
                    (i!=j)  && (abs(i-j) < 2)){
                if(pdb.atoms[i].id > pdb.atoms[j].id){
                    top.harmonic[top.maxHarmonicPerMonomer*i + top.harmonicCount[i]] = j; 
                }
                else{
                    top.harmonic[top.maxHarmonicPerMonomer*i + top.harmonicCount[i]] = -j; 
                }
                top.harmonicCount[i]++;
            }
        }
    }

	if (!par.is_assembly){
		 // Longitudal exponential
	    for(i = 0; i < par.Ntot * par.Ntr; i++) top.longitudinalCount[i] = 0;

	    for (int traj = 0; traj < par.Ntr; traj++){
	        for(i = 0; i < par.Ntot; i++){
	            if (top.extra[i + traj * par.Ntot]) { //////////////////////////////////////////////////
	                continue;
	            }
	            for(j = 0; j < par.Ntot; j++){
	                if( (pdb.atoms[i].resid == (pdb.atoms[j].resid + 1)) &&
	                    (pdb.atoms[i].chain == pdb.atoms[j].chain) &&
	                    (pdb.atoms[i].id == (pdb.atoms[j].id + 1)) && (i != j) ){
	                    top.longitudinalCount[i + par.Ntot * traj]++;
	                    if(top.maxLongitudinalPerMonomer < top.longitudinalCount[i + par.Ntot * traj])
	                            top.maxLongitudinalPerMonomer = top.longitudinalCount[i + par.Ntot * traj];
	                }
	            }
	        }
	    }
	    top.longitudinal = (int*)calloc(top.maxLongitudinalPerMonomer * par.Ntot * par.Ntr, sizeof(int));
	    
	    for(i = 0; i < par.Ntot * par.Ntr; i++) top.longitudinalCount[i] = 0;

	    for(int traj = 0; traj < par.Ntr; traj++){
	        for(i = 0; i < par.Ntot; i++){
	            if (top.extra[i + traj * par.Ntot]) {
	                continue;                              /////////////////////////////////////////////////
	            }
	            for(j = 0; j < par.Ntot; j++){
	                if(abs(pdb.atoms[i].resid - pdb.atoms[j].resid) == 1 &&
	                    pdb.atoms[i].chain == pdb.atoms[j].chain &&
	                    abs(pdb.atoms[i].id -  pdb.atoms[j].id) == 1){
	                        if(pdb.atoms[i].id > pdb.atoms[j].id){
	                            top.longitudinal[top.maxLongitudinalPerMonomer * par.Ntot * traj + i * top.maxLongitudinalPerMonomer + top.longitudinalCount[i + par.Ntot * traj]] = j;
	                            //top.longitudinal[top.maxLongitudinalPerMonomer*i + top.longitudinalCount[i]] = j; // << BUG: indexation
	                        }
	                        else{
	                            top.longitudinal[top.maxLongitudinalPerMonomer * par.Ntot * traj + i * top.maxLongitudinalPerMonomer + top.longitudinalCount[i + par.Ntot * traj]] = -j;
	                            //top.longitudinal[top.maxLongitudinalPerMonomer*i + top.longitudinalCount[i]] = -j; /* BUG: same here */
	                        }
	                        top.longitudinalCount[i + par.Ntot * traj]++;
	                }
	            }
	        }

	    }
	    
	    // Lateral exponential
	    for(i = 0; i < par.Ntot * par.Ntr; i++) top.lateralCount[i] = 0;
	    
	    real dr, xi, xj, yi, yj, zi, zj;
	    real cos_fii, cos_fij, sin_fii, sin_fij,
	          cos_psii, cos_psij, sin_psii, sin_psij,
	          cos_thetai, cos_thetaj, sin_thetai, sin_thetaj;
	    real xp1 = xp1_def;
	    real yp1 = yp1_def;
	    real zp1 = zp1_def;
	    real xp2 = xp2_def;
	    real yp2 = yp2_def;
	    real zp2 = zp2_def;
	    for(int traj = 0; traj < par.Ntr; traj++){

	        for(i = 0; i < par.Ntot; i++){
	            if (top.extra[i + traj * par.Ntot]) {
	                continue;
	            }
	            xi = r[i].x;    
	            yi = r[i].y;    
	            zi = r[i].z;
	            sin_fii = sinf(r[i].fi);
	            cos_fii = cosf(r[i].fi);
	            sin_psii = sinf(r[i].psi);
	            cos_psii = cosf(r[i].psi);      
	            sin_thetai = sinf(r[i].theta);
	            cos_thetai = cosf(r[i].theta);  
	            for(j = 0; j < par.Ntot; j++){
	                xj = r[j].x;    
	                yj = r[j].y;    
	                zj = r[j].z;
	                sin_fij = sinf(r[j].fi);
	                cos_fij = cosf(r[j].fi);
	                sin_psij = sinf(r[j].psi);
	                cos_psij = cosf(r[j].psi);      
	                sin_thetaj = sinf(r[j].theta);
	                cos_thetaj = cosf(r[j].theta);  
	                for(int ind = 0; ind < 2; ind++)
	                {
	                    if(ind == 0)
	                    {
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
	                    
	                    if(dr < r_mon){
	                        top.lateralCount[i + traj * par.Ntot]++;
	                        if(top.maxLateralPerMonomer < top.lateralCount[i + traj * par.Ntot])
	                                top.maxLateralPerMonomer = top.lateralCount[i + traj * par.Ntot];
	                    }
	                }
	            }   
	        }
	    }
	    
	    top.lateral = (int*)calloc(top.maxLateralPerMonomer * par.Ntot * par.Ntr, sizeof(int));
	    for(i = 0; i < par.Ntot * par.Ntr; i++) top.lateralCount[i] = 0;

	    for(int traj = 0; traj < par.Ntr; traj++){
	        for(i = 0; i < par.Ntot; i++){
	             if (top.extra[i + traj * par.Ntot]) {
	                continue;                                  /////////////////////////////////////////////////
	            }
	            xi = r[i].x;    
	            yi = r[i].y;    
	            zi = r[i].z;
	            sin_fii = sinf(r[i].fi);
	            cos_fii = cosf(r[i].fi);
	            sin_psii = sinf(r[i].psi);
	            cos_psii = cosf(r[i].psi);      
	            sin_thetai = sinf(r[i].theta);
	            cos_thetai = cosf(r[i].theta);  

	            for(j = 0; j < par.Ntot; j++){
	                xj = r[j].x;    
	                yj = r[j].y;    
	                zj = r[j].z;
	                sin_fij = sinf(r[j].fi);
	                cos_fij = cosf(r[j].fi);
	                sin_psij = sinf(r[j].psi);
	                cos_psij = cosf(r[j].psi);      
	                sin_thetaj = sinf(r[j].theta);
	                cos_thetaj = cosf(r[j].theta);  
	                for(int ind = 0; ind < 2; ind++)
	                {   
	                    if(ind == 0)
	                    {
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
	                    if(dr < r_mon){
	                        if(ind == 1){
	                            top.lateral[top.maxLateralPerMonomer * par.Ntot * traj + i * top.maxLateralPerMonomer + top.lateralCount[i + par.Ntot * traj]] = j;
	                        }else{
	                            top.lateral[top.maxLateralPerMonomer * par.Ntot * traj + i * top.maxLateralPerMonomer + top.lateralCount[i + par.Ntot * traj]] = -j;
	                        }
	                        top.lateralCount[i + par.Ntot * traj]++;
	                    }
	                }
	            }
	        }
	    }
	}   
   
    printf("done building topology without LJ.\n");
}

void saveCoordPDB(const char* pdbfilename_xyz, const char* pdbfilename_ang){
    int i;
    for(i = 0; i < par.Ntot; i++){
        pdb.atoms[i].x = r[i].x;
        pdb.atoms[i].y = r[i].y;
        pdb.atoms[i].z = r[i].z;
    }
    writePDB(pdbfilename_xyz, &pdb);
    for(i = 0; i < par.Ntot; i++){
        pdb.atoms[i].x = r[i].fi;
        pdb.atoms[i].y = r[i].psi;
        pdb.atoms[i].z = r[i].theta;
    }
    writePDB(pdbfilename_ang, &pdb);
}

void saveCoordDCD(){
    int i;
    for(int j=0; j < par.Ntr; j++){
        for(i = 0; i < par.Ntot; i++){
            dcd.frame.X[i] = r[i+par.Ntot*j].x;
            dcd.frame.Y[i] = r[i+par.Ntot*j].y;
            dcd.frame.Z[i] = r[i+par.Ntot*j].z;
        }
        dcdOpenAppend(&dcd, par.dcdFilename_xyz[j]);
        dcdWriteFrame(dcd);
        dcdClose(dcd);
        for(i = 0; i < par.Ntot; i++){
            dcd.frame.X[i] = r[i+par.Ntot*j].fi;
            dcd.frame.Y[i] = r[i+par.Ntot*j].psi;
            dcd.frame.Z[i] = r[i+par.Ntot*j].theta;
        }
        dcdOpenAppend(&dcd, par.dcdFilename_ang[j]);
        dcdWriteFrame(dcd);
        dcdClose(dcd);
    }
}

void appendCoordPDB(){

    /*
typedef struct {

  int id;
  char   name[5], chain, resName[4], altLoc;
  int    resid;
  double x, y, z;

  double occupancy;
  double beta;

} PDBAtom;

    */
    int i;
    for(int j = 0; j < par.Ntr; j++){
        for(i = 0; i < par.Ntot; i++){
            coordspdb.atoms[i+par.Ntot*j].x = (float)r[i+par.Ntot*j].x;
            coordspdb.atoms[i+par.Ntot*j].y = (float)r[i+par.Ntot*j].y;
            coordspdb.atoms[i+par.Ntot*j].z = (float)r[i+par.Ntot*j].z;
            coordspdb.atoms[i+par.Ntot*j].id = i + par.Ntot * j;
            coordspdb.atoms[i+par.Ntot*j].beta = (double)j;
            coordspdb.atoms[i+par.Ntot*j].occupancy = (double)top.gtp[i+par.Ntot*j];
        }
    }    
    appendPDB("dcd/hydrolysis.pdb", &coordspdb);
}


void ReadFromDCD(Parameters par, Topology top, char* dcdfilename_xyz, char* dcdfilename_ang)
{

    DCD dcd_xyz, dcd_ang;
    dcd_xyz.frame.X = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_xyz.frame.Y = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_xyz.frame.Z = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_ang.frame.X = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_ang.frame.Y = (float*)calloc(pdb.atomCount, sizeof(float));
    dcd_ang.frame.Z = (float*)calloc(pdb.atomCount, sizeof(float));
    
    dcdOpenRead(&dcd_xyz, dcdfilename_xyz);
    dcdOpenRead(&dcd_ang, dcdfilename_ang);
    dcdReadHeader(&dcd_xyz);
    dcdReadHeader(&dcd_ang);
    while( ( dcdReadFrame(&dcd_xyz) == 0 ) && ( dcdReadFrame(&dcd_ang) == 0 ) )
    {
        for(int i = 0; i < pdb.atomCount; i++)
        {
            r[i].x = dcd_xyz.frame.X[i];
            r[i].y = dcd_xyz.frame.Y[i];
            r[i].z = dcd_xyz.frame.Z[i];
            r[i].fi    = dcd_ang.frame.X[i];
            r[i].psi   = dcd_ang.frame.Y[i];
            r[i].theta = dcd_ang.frame.Z[i];
        }
    }
}

void readRestart()
{
    for(int traj = 0; traj < par.Ntr; traj++)
    {
        readXYZ(par.restart_xyzFilename[traj], &xyz);
        for(int i = 0; i < par.Ntot; i++)
        {
            r[i + par.Ntot*traj].x = xyz.atoms[i].x;
            r[i + par.Ntot*traj].y = xyz.atoms[i].y;
            r[i + par.Ntot*traj].z = xyz.atoms[i].z;
        }
        readXYZ(par.restart_angFilename[traj], &xyz);
        for(int i = 0; i < par.Ntot; i++)
        {
            r[i + par.Ntot*traj].fi = xyz.atoms[i].x;
            r[i + par.Ntot*traj].psi = xyz.atoms[i].y;
            r[i + par.Ntot*traj].theta = xyz.atoms[i].z;
        }
    }
    free(xyz.atoms);
}

void writeRestart(long long int step)
{
    xyz.atomCount = par.Ntot;
    xyz.atoms = (XYZAtom*)calloc(xyz.atomCount, sizeof(XYZAtom));
    for(int traj = 0; traj < par.Ntr; traj++)
    {
        for(int i = 0; i < par.Ntot; i++)
        {
            xyz.atoms[i].x = r[i + par.Ntot*traj].x;
            xyz.atoms[i].y = r[i + par.Ntot*traj].y;
            xyz.atoms[i].z = r[i + par.Ntot*traj].z;
            xyz.atoms[i].name = pdb.atoms[i].name[0];
        }
        writeXYZ(par.restart_xyzFilename[traj], &xyz);
        for(int i = 0; i < par.Ntot; i++)
        {
            xyz.atoms[i].x = r[i + par.Ntot*traj].fi;
            xyz.atoms[i].y = r[i + par.Ntot*traj].psi;
            xyz.atoms[i].z = r[i + par.Ntot*traj].theta;
            xyz.atoms[i].name = pdb.atoms[i].name[0];
        }
        writeXYZ(par.restart_angFilename[traj], &xyz);
    }
    free(xyz.atoms);
    FILE *keyf = safe_fopen(par.restartkey, "w");
    if (fprintf(keyf,"%lld", step) == 0)
        DIE("Reading restartkey %s: unable to write firststep", &(par.restartkey));
    fclose(keyf);
}


void AssemblyInit()
{
    free(top.lateral);
    free(top.longitudinal);
    top.maxLateralPerMonomer = 16;
    top.maxLongitudinalPerMonomer = 8;

    top.lateral = (int*)calloc(top.maxLateralPerMonomer * par.Ntot * par.Ntr, sizeof(int));
    top.longitudinal = (int*)calloc(top.maxLongitudinalPerMonomer * par.Ntot * par.Ntr, sizeof(int));
    for(int i = 0; i < par.Ntot * par.Ntr; i++)
    {
        top.longitudinalCount[i] = 0;
        top.lateralCount[i] = 0;

    }
}