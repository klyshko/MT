/*
 * bdhitea_kernel.cu
 *
 *  Created on: Mar 4, 2013
 *	  Author: alekseenko
 * 
 * All equations referenced here are from Geyer & Winter, 2009 [doi:10.1063/1.3089668]
 */



#include "bdhitea.cuh"
//#define PRINT_HI_TENSORS // Enable to print tensor values. Should never be used in production. Calling `cudaDeviceSetLimit(cudaLimitPrintfFifoSize, 900000000);` is recommended if you do not want to lose any values!


__global__ void integrateTea_prepare(Coord* d_f, Coord* d_r){
	// Precalculate random forces and apply pulling forces
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_par.Ntot * c_par.Ntr){
		// Random force
		const float var = sqrtf(2.0f * KB * c_par.Temp * c_par.gammaR / c_par.dt);//c_langevin.var;                 /////// sqrt(2 kb T gamma/dt)
		float4 df = rforce(d_i);
		df.x *= var;
		df.y *= var;
		df.z *= var;
		c_tea.rforce[d_i] = df;
		
		float4 f = make_float4(d_f[d_i].x, d_f[d_i].y, d_f[d_i].z, 0.f );//c_gsop.d_forces[d_i];
		// Copy forces and coordinates to auxillary arrays to avoid races during integration phase
		c_tea.mforce[d_i] = f;
		d_f[d_i].x = 0.f;
		d_f[d_i].y = 0.f;
		d_f[d_i].z = 0.f;
		c_tea.coords[d_i] = make_float4(d_r[d_i].x, d_r[d_i].y, d_r[d_i].z, 0.f);
	}
}

__device__ inline float6 integrateTea_RPY(const float4& dr){
	// This functions requires `dr` to be normalized and have its original length in its W component
	// Returns Dij(dr)/Dii, where Dij is Rotne-Prager-Yamakawa tensor submatrix, eq. (3-5)
	float6 ret;
	const float ra = dr.w / c_tea.a;
	float coeffrr, coeffii;
	if (ra > 2.f){
		coeffrr = 0.75f/ra * (1.f - 2.f/ra/ra);
		coeffii = 0.75f/ra * (1.f + 2.f/3.f/ra/ra);
	}else{
		coeffrr = 3.f*ra/32.f;
		coeffii = 1.f - 9.f*ra/32.f;
	}
	ret._XX = dr.x*dr.x*coeffrr + coeffii;
	ret._XY = dr.x*dr.y*coeffrr;
	ret._XZ = dr.x*dr.z*coeffrr;
	ret._YY = dr.y*dr.y*coeffrr + coeffii;
	ret._YZ = dr.y*dr.z*coeffrr;
	ret._ZZ = dr.z*dr.z*coeffrr + coeffii;
	return ret;
}


__device__ inline float4 integrateTea_epsilon_local(const float4& coord1, const int idx2, Coord* d_r){
	// Function calculates various statistics for sub-tensor responsible for interactions between two given particles
	// .x, .y, .z --- sum of squares of normalized tensor components, 3 rows, used for C_i in eq. (19)
	// .w --- sum of normalized tensor components (\sum Dij/Dii), used for \epsilon in eq. (22)
	float4 dr = make_float4(d_r[idx2].x, d_r[idx2].y, d_r[idx2].z, 0.f);
	dr.x -= coord1.x;
	dr.y -= coord1.y;
	dr.z -= coord1.z;
	//printf("%f\n", dr.x);
	dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
	dr.x /= dr.w;
	dr.y /= dr.w;
	dr.z /= dr.w;
	//printf("%d kkk %f %f %f  %f\n",idx2 , dr.x, dr.y, dr.z, dr.w);
	float6 d = integrateTea_RPY(dr);
	dr.w = d._XX + 2*d._XY + 2*d._XZ + d._YY + 2*d._YZ + d._ZZ; // Sum of all 3*3 tensor components
	dr.x = d._XX*d._XX + d._XY*d._XY + d._XZ*d._XZ;
	dr.y = d._YX*d._YX + d._YY*d._YY + d._YZ*d._YZ;
	dr.z = d._ZX*d._ZX + d._ZY*d._ZY + d._ZZ*d._ZZ;
	return dr;
}


__global__ void integrateTea_epsilon_unlisted(Coord* d_r){
	// Like integrateTea_epsilon, but calculate all-vs-all
	const int d_i = blockIdx.x * blockDim.x + threadIdx.x;
	if(d_i < c_par.Ntot * c_par.Ntr){
		const int i0 = ((int)(d_i / c_tea.Ntot))*c_tea.Ntot;
		int i;
		
		float4 coord = make_float4(d_r[d_i].x, d_r[d_i].y, d_r[d_i].z, 0.f);
		float4 sum = make_float4(0.f, 0.f, 0.f, 0.f);
		for(i = i0; i < i0 + c_tea.Ntot; i++){
			if (i != d_i && !c_top.extra[i] && !c_top.extra[d_i]){
				sum += integrateTea_epsilon_local(coord, i, d_r);
			}
		}
		c_tea.d_ci[d_i] = make_float4(sum.x, sum.y, sum.z, 0.f);
		c_tea.d_epsilon[d_i] = sum.w; // Should later be divided by number of non-diagonal degrees-of-freedom in single trajectory
	}
}


__device__ inline float4 integrateTea_force(const float4& coord1, const int idx2, const float3& ci, const int idx1){
	// Calculate the effective force acting on particle with coordinates `coord1` from particle with index `idx2`
	// eq. (13,14,19)

	float4 dr = c_tea.coords[idx2];
	dr.x -= coord1.x;
	dr.y -= coord1.y;
	dr.z -= coord1.z;
	dr.w = sqrtf(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
	dr.x /= dr.w;
	dr.y /= dr.w;
	dr.z /= dr.w;
	float4 f = c_tea.mforce[idx2];
	float4 r = c_tea.rforce[idx2];

	f.x += r.x * ci.x;
	f.y += r.y * ci.y;
	f.z += r.z * ci.z;
	float6 D = integrateTea_RPY(dr);
#ifdef PRINT_HI_TENSORS
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+0, D._XX);
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+1, D._XY);
	printf("d[%d, %d]=%f\n", 3*idx1+0, 3*idx2+2, D._XZ);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+0, D._YX);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+1, D._YY);
	printf("d[%d, %d]=%f\n", 3*idx1+1, 3*idx2+2, D._YZ);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+0, D._ZX);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+1, D._ZY);
	printf("d[%d, %d]=%f\n", 3*idx1+2, 3*idx2+2, D._ZZ);
	printf("b[%d, %d]=%f\n", 3*idx1+0, 3*idx2+0, D._XX*ci.x);
	printf("b[%d, %d]=%f\n", 3*idx1+0, 3*idx2+1, D._XY*ci.x);
	printf("b[%d, %d]=%f\n", 3*idx1+0, 3*idx2+2, D._XZ*ci.x);
	printf("b[%d, %d]=%f\n", 3*idx1+1, 3*idx2+0, D._YX*ci.y);
	printf("b[%d, %d]=%f\n", 3*idx1+1, 3*idx2+1, D._YY*ci.y);
	printf("b[%d, %d]=%f\n", 3*idx1+1, 3*idx2+2, D._YZ*ci.y);
	printf("b[%d, %d]=%f\n", 3*idx1+2, 3*idx2+0, D._ZX*ci.z);
	printf("b[%d, %d]=%f\n", 3*idx1+2, 3*idx2+1, D._ZY*ci.z);
	printf("b[%d, %d]=%f\n", 3*idx1+2, 3*idx2+2, D._ZZ*ci.z);
#endif
	return make_float4( D._XX*f.x + D._XY*f.y + D._XZ*f.z, 
						D._YX*f.x + D._YY*f.y + D._YZ*f.z, 
						D._ZX*f.x + D._ZY*f.y + D._ZZ*f.z, 0.f);
}

__global__ void integrateTea_kernel_unlisted(Coord* d_f, Coord* d_r){
	// Pairist-free version of  integrateTea_kernel
	const int d_i = blockIdx.x*blockDim.x + threadIdx.x;
	if(d_i < c_par.Ntot * c_par.Ntr ){
		int i;
		float4 coord = c_tea.d_ci[d_i]; // Not coord yet!
		float4 f = c_tea.mforce[d_i];
		float4 df = c_tea.rforce[d_i];
		const int tr = d_i / c_tea.Ntot;
		const float beta_ij = c_tea.d_beta_ij[tr];
		float3 ci;
		// Make ci to be actual C_i
		coord.w = beta_ij*beta_ij;
		ci.x = 1.f/sqrtf(1.f + coord.w * coord.x);
		ci.y = 1.f/sqrtf(1.f + coord.w * coord.y);
		ci.z = 1.f/sqrtf(1.f + coord.w * coord.z);
		coord = c_tea.coords[d_i]; // And now it's actually bead coordinates
		f.x += df.x * ci.x;
		f.y += df.y * ci.y;
		f.z += df.z * ci.z;
		ci.x *= beta_ij;
		ci.y *= beta_ij;
		ci.z *= beta_ij;

		// Calculate effective force
		const int i0 = ((int)(d_i / c_tea.Ntot))*c_tea.Ntot;
		for(i = i0; i < i0 + c_tea.Ntot; i++){
			if (i == d_i || c_top.extra[d_i] || c_top.extra[i]) continue;
			df = integrateTea_force(coord, i, ci, d_i);
			f.x += df.x;
			f.y += df.y;
			f.z += df.z;
		}
		
		// Integration step
		// We've replaced all forces with their `effective` counterparts, so this part of integration process stays the same as in simple langevin integrator
		const float mult = c_par.dt / c_par.gammaR ; ///c_langevin.hOverZeta;
		const float3 dr = make_float3(mult*f.x, mult*f.y, mult*f.z);
		coord.x += dr.x;
		coord.y += dr.y;
		coord.z += dr.z;

		Coord ri = d_r[d_i];
		Coord fi = d_f[d_i];

		float4 rf_ang = make_float4(0.0, 0.0, 0.0, 0.0);
		rf_ang = rforce(d_i + c_par.Ntot*c_par.Ntr);

		if(!c_top.fixed[d_i % c_par.Ntot] && !c_top.extra[d_i]){
			ri.x = coord.x;
			ri.y = coord.y;
			ri.z = coord.z;

			ri.fi    += (c_par.dt/(c_par.gammaTheta * c_par.alpha))*fi.fi  + (c_par.varTheta * sqrt(c_par.freeze_temp / c_par.alpha))*rf_ang.x;
			ri.psi   += (c_par.dt/(c_par.gammaTheta * c_par.alpha))*fi.psi + (c_par.varTheta * sqrt(c_par.freeze_temp / c_par.alpha))*rf_ang.y;
			ri.theta += (c_par.dt/c_par.gammaTheta)*fi.theta + c_par.varTheta*rf_ang.z;
		} 
		
		d_r[d_i] = ri;

		// Update energies
		/*coord = c_gsop.d_energies[d_i];
		coord.w += c_langevin.tempNorm*(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);
		c_gsop.d_energies[d_i] = coord;*/
	}
}

#undef _XX
#undef _XY
#undef _XZ
#undef _YX
#undef _YY
#undef _YZ
#undef _ZX
#undef _ZY
#undef _ZZ

