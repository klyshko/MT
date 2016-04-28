/*
 * bdhitea.cu
 *
 *  Created on: Mar 01, 2013
 *	  Author: alekseenko
 * 
 * The description of algorithm is available in Geyer & Winter, 2009 [doi:10.1063/1.3089668]
 * All the equations referenced here are from the very same paper
 */

#include "bdhitea.cuh"

void initTeaIntegrator(){
	
	//initLangevinIntegrator();
	//teaIntegrator.h = par.dt;
	tea.Ntot = par.Ntot;
	 
	cudaMalloc(&tea.rforce, par.Ntr * par.Ntot * sizeof(float4));
	cudaMalloc(&tea.mforce, par.Ntr * par.Ntot * sizeof(float4));
	cudaMalloc(&tea.coords, par.Ntr * par.Ntot * sizeof(float4));

	cudaMalloc(&tea.d_epsilon, par.Ntot * par.Ntr * sizeof(float));
	cudaMalloc(&tea.d_ci, par.Ntot * par.Ntr * sizeof(float4));
	cudaMalloc(&tea.d_beta_ij, par.Ntr * sizeof(float));
    tea.h_epsilon = (float*) malloc(par.Ntot * par.Ntr * sizeof(float));
	tea.h_beta_ij = (float*) malloc( par.Ntr * sizeof(float));

	checkCUDAError("tea memory allocation");
	cudaMemcpyToSymbol(c_tea, &tea, sizeof(Tea), 0, cudaMemcpyHostToDevice);
	checkCUDAError("tea constant copy");

    //createTeaUpdater();
    printf("TEA integrator initialized, a = %f A, freq = %d steps; capricious mode: %s; using pairlist: %s; block size: %d\n", tea.a, tea.epsilon_freq, (tea.capricious ? "on" : "off"), (tea.unlisted ? "no" : "yes"), BLOCK_SIZE);   
}

void integrateTea(){
	// Pregenerate random forces
	integrateTea_prepare<<<(par. Ntr * par.Ntot) / BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_f, d_r);
	integrateTea_kernel_unlisted<<<(par. Ntr * par.Ntot) / BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_f, d_r);
	checkCUDAError("integrate TEA prepare and kernel");
}

void deleteTeaIntegrator(){
	//deleteLangevinIntegrator();
	cudaFree(tea.rforce);
	cudaFree(tea.mforce);
	cudaFree(tea.coords);
   	cudaFree(tea.d_epsilon);
	cudaFree(tea.d_beta_ij);
	cudaFree(tea.d_ci);
    free(tea.h_epsilon);
	free(tea.h_beta_ij);

}

void updateTea(long long int step){
	const int update_epsilon = (step % tea.epsilon_freq) == 0;
	const int N = par.Ntot;

	int* Nnoextra = (int*)calloc(par.Ntr, sizeof(int));
	for(int tr = 0; tr < par.Ntr; tr++){
		for(int j = 0; j < par.Ntot; j++){
			if (!top.extra[j + par.Ntot * tr]){
				Nnoextra[tr]++;
			}
		}
		//printf("Nnoextra for [%d] is %d\n", tr, Nnoextra[tr]);
	}

	if (update_epsilon){
		// Calculate relative coupling
		integrateTea_epsilon_unlisted<<<(par.Ntr * par.Ntot) / BLOCK_SIZE + 1, BLOCK_SIZE>>>(d_r);
		// Dowload epsilon`s
		cudaMemcpy(tea.h_epsilon, tea.d_epsilon, par.Ntr * par.Ntot * sizeof(float), cudaMemcpyDeviceToHost);
		checkCUDAError("copy epsilon from device");

		//printf("epsilon: [ ");
		for (int t = 0; t < par.Ntr; ++t){
			double epsilon = 0.0;

			for (int i = 0; i < N; ++i){
				//if (i%10==0) printf("[%d] %f\n", i, tea.h_epsilon[t*N + i]);
				epsilon += tea.h_epsilon[t*N + i];
			}

			//epsilon /= 3.*N*(3.*N - 3.); // Averaging, off-diagonal elements only
			epsilon /= 3.*Nnoextra[t]*(3.*Nnoextra[t] - 3.);
			if (epsilon > 1.0){
				if (tea.capricious){
					printf("HI tensor is not diagonal enough for trajectory %d: epsilon = %lf -> 1.0!\n", t, epsilon);
					exit(-1);
				}
				epsilon = 1.0;
			}
			if (epsilon > tea.epsmax){
				printf("HI tensor is not diagonal enough for trajectory %d: epsilon = %lf > %f = tea_epsmax!\n", t, epsilon, tea.epsmax);
				exit(-1);
			}
			//double a = (3.*N-1.)*epsilon*epsilon - (3.*N-2.)*epsilon;
			double a = (3.*Nnoextra[t]-1.)*epsilon*epsilon - (3.*Nnoextra[t]-2.)*epsilon;
			if (fabs(a) < 1e-7){ // To avoid 0/0 division in eq. (26) we explicitly handle small a's
				tea.h_beta_ij[t] = .5f;
				if(tea.capricious && tea.a > 0.0f){
					printf("HI tensor is too diagonal for trajectory %d: a = %lf, beta: %lf -> 0.5!\n", t, a, (1. - sqrt(1. - a)) / a);
					exit(-1);
				}
			} else {
				tea.h_beta_ij[t] = (1. - sqrt(1. - a)) / a; // eq. (26)
			}
			tea.h_epsilon[t] = epsilon; // We slowly overwrite the beginning of h_epsilon with per-trajectory epsilons to later output them
			//printf("%lf ", epsilon);
		}
		//printf("]\n");
		cudaMemcpy(tea.d_beta_ij, tea.h_beta_ij, par.Ntr * sizeof(float), cudaMemcpyHostToDevice);
		checkCUDAError("copy betaij to device");
	}
}


