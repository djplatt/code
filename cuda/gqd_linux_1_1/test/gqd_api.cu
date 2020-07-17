#ifndef _GQD_API_CU_
#define _GQD_API_CU_

#include "gqd_api.h"

#include <stdio.h>
#include <stdlib.h>
#include <cutil_inline.h>

#include "cuda_header.cu"
#include "gqd.cu"
#include "gdd.cu"
#include "map.cu"

/**
* c[i] = a[i] + b[i]
*/
template<class T> 
void gpu_add(const T* h_a, const T* h_b, T* h_c, const int len) {
	
	//copy memory to GPU
	T* d_a = NULL;
	GPUMALLOC((void**)&d_a, sizeof(T)*len);
	TOGPU(d_a, h_a, sizeof(T)*len);
	T* d_b = NULL;
	GPUMALLOC((void**)&d_b, sizeof(T)*len);
	TOGPU(d_b, h_b, sizeof(T)*len);
	T* d_c = NULL;
	GPUMALLOC((void**)&d_c, sizeof(T)*len);

	//kernel
	const int numBlock = 128;
	const int numThread = 128;
	unsigned int timer = 0;
	startTimer(timer);
	map_add_kernel<T><<<numBlock, numThread>>>(d_a, d_b, d_c, len);
	cutilCheckMsg("map_add_kernel");
	CUDA_SAFE_CALL(cudaThreadSynchronize());
	endTimer(timer, "GPU add kernel");

	//copy results from GPU
	FROMGPU(h_c, d_c, sizeof(T)*len);

	//clean up
	GPUFREE(d_a);
	GPUFREE(d_b);
	GPUFREE(d_c);
}

void gpu_dd_add(const GPU_dd* h_a, const GPU_dd* h_b, 
		GPU_dd* h_c, const int len) {
	gpu_add<GPU_dd>(h_a, h_b, h_c, len);
}

void gpu_qd_add(const GPU_qd* h_a, const GPU_qd* h_b,
                GPU_qd* h_c, const int len) {
        gpu_add<GPU_qd>(h_a, h_b, h_c, len);
}

template<class T> 
void gpu_sub(const T* h_a, const T* h_b, T* h_c, const int len) {
	
	//copy memory to GPU
	T* d_a = NULL;
	GPUMALLOC((void**)&d_a, sizeof(T)*len);
	TOGPU(d_a, h_a, sizeof(T)*len);
	T* d_b = NULL;
	GPUMALLOC((void**)&d_b, sizeof(T)*len);
	TOGPU(d_b, h_b, sizeof(T)*len);
	T* d_c = NULL;
	GPUMALLOC((void**)&d_c, sizeof(T)*len);

	//kernel
	const int numBlock = 128;
	const int numThread = 128;
	unsigned int timer = 0;
	startTimer(timer);
	map_sub_kernel<T><<<numBlock, numThread>>>(d_a, d_b, d_c, len);
	cutilCheckMsg("map_sub_kernel");
	CUDA_SAFE_CALL(cudaThreadSynchronize());
	endTimer(timer, "GPU sub kernel");

	//copy results from GPU
	FROMGPU(h_c, d_c, sizeof(T)*len);

	//clean up
	GPUFREE(d_a);
	GPUFREE(d_b);
	GPUFREE(d_c);
}

void gpu_dd_sub(const GPU_dd* h_a, const GPU_dd* h_b, 
		GPU_dd* h_c, const int len) {
	gpu_sub<GPU_dd>(h_a, h_b, h_c, len);
}

/**
* c[i] = a[i] * b[i]
*/
template<class T> 
void gpu_mul(const T* h_a, const T* h_b, T* h_c, const int len) {
	
	//copy memory to GPU
	T* d_a = NULL;
	GPUMALLOC((void**)&d_a, sizeof(T)*len);
	TOGPU(d_a, h_a, sizeof(T)*len);
	T* d_b = NULL;
	GPUMALLOC((void**)&d_b, sizeof(T)*len);
	TOGPU(d_b, h_b, sizeof(T)*len);
	T* d_c = NULL;
	GPUMALLOC((void**)&d_c, sizeof(T)*len);

	//kernel
	const int numBlock = 128;
	const int numThread = 128;
	unsigned int timer = 0;
	startTimer(timer);
	map_mul_kernel<T><<<numBlock, numThread>>>(d_a, d_b, d_c, len);
	cutilCheckMsg("map_mul_kernel");
	CUDA_SAFE_CALL(cudaThreadSynchronize());
	endTimer(timer, "GPU mul kernel");

	//copy results from GPU
	FROMGPU(h_c, d_c, sizeof(T)*len);

	//clean up
	GPUFREE(d_a);
	GPUFREE(d_b);
	GPUFREE(d_c);
}

void gpu_dd_mul(const GPU_dd* h_a, const GPU_dd* h_b, 
		GPU_dd* h_c, const int len) {
	gpu_mul<GPU_dd>(h_a, h_b, h_c, len);
}

void gpu_qd_mul(const GPU_qd* h_a, const GPU_qd* h_b,
                GPU_qd* h_c, const int len) {
        gpu_mul<GPU_qd>(h_a, h_b, h_c, len);
}




/**
* c[i] = a[i] / b[i]
*/
template<class T>
void gpu_div(const T* h_a, const T* h_b, T* h_c, const int len) {

        //copy memory to GPU
        T* d_a = NULL;
        GPUMALLOC((void**)&d_a, sizeof(T)*len);
        TOGPU(d_a, h_a, sizeof(T)*len);
        T* d_b = NULL;
        GPUMALLOC((void**)&d_b, sizeof(T)*len);
        TOGPU(d_b, h_b, sizeof(T)*len);
        T* d_c = NULL;
        GPUMALLOC((void**)&d_c, sizeof(T)*len);

        //kernel
        const int numBlock = 128;
        const int numThread = 128;
        unsigned int timer = 0;
        startTimer(timer);
        map_div_kernel<T><<<numBlock, numThread>>>(d_a, d_b, d_c, len);
        cutilCheckMsg("map_div_kernel");
        CUDA_SAFE_CALL(cudaThreadSynchronize());
        endTimer(timer, "GPU div kernel");

        //copy results from GPU
        FROMGPU(h_c, d_c, sizeof(T)*len);

        //clean up
        GPUFREE(d_a);
        GPUFREE(d_b);
        GPUFREE(d_c);
}

void gpu_dd_div(const GPU_dd* h_a, const GPU_dd* h_b,
                GPU_dd* h_c, const int len) {
        gpu_div<GPU_dd>(h_a, h_b, h_c, len);
}

void gpu_qd_div(const GPU_qd* h_a, const GPU_qd* h_b,
                GPU_qd* h_c, const int len) {
        gpu_div<GPU_qd>(h_a, h_b, h_c, len);
}


/**
* c[i] = sqrt(a[i])
*/
template<class T>
void gpu_sqrt(const T* h_a, T* h_c, const int len) {

        //copy memory to GPU
        T* d_a = NULL;
        GPUMALLOC((void**)&d_a, sizeof(T)*len);
        TOGPU(d_a, h_a, sizeof(T)*len);
        T* d_c = NULL;
        GPUMALLOC((void**)&d_c, sizeof(T)*len);

        //kernel
        const int numBlock = 128;
        const int numThread = 128;
        unsigned int timer = 0;
        startTimer(timer);
        map_sqrt_kernel<T><<<numBlock, numThread>>>(d_a, d_c, len);
        cutilCheckMsg("map_sqrt_kernel");
        CUDA_SAFE_CALL(cudaThreadSynchronize());
        endTimer(timer, "GPU sqrt kernel");

        //copy results from GPU
        FROMGPU(h_c, d_c, sizeof(T)*len);

        //clean up
        GPUFREE(d_a);
        GPUFREE(d_c);
}

void gpu_dd_sqrt(const GPU_dd* h_a, GPU_dd* h_c, const int len) {
        gpu_sqrt<GPU_dd>(h_a, h_c, len);
}

void gpu_qd_sqrt(const GPU_qd* h_a, GPU_qd* h_c, const int len) {
        gpu_sqrt<GPU_qd>(h_a, h_c, len);
}


/**
* c[i] = exp(a[i])
*/
template<class T>
void gpu_exp(const T* h_a, T* h_c, const int len) {

        //copy memory to GPU
        T* d_a = NULL;
        GPUMALLOC((void**)&d_a, sizeof(T)*len);
        TOGPU(d_a, h_a, sizeof(T)*len);
        T* d_c = NULL;
        GPUMALLOC((void**)&d_c, sizeof(T)*len);

        //kernel
        const int numBlock = 128;
        const int numThread = 128;
        unsigned int timer = 0;
        startTimer(timer);
        map_exp_kernel<T><<<numBlock, numThread>>>(d_a, d_c, len);
        cutilCheckMsg("map_exp_kernel");
        CUDA_SAFE_CALL(cudaThreadSynchronize());
        endTimer(timer, "GPU exp kernel");

        //copy results from GPU
        FROMGPU(h_c, d_c, sizeof(T)*len);

        //clean up
        GPUFREE(d_a);
        GPUFREE(d_c);
}

void gpu_dd_exp(const GPU_dd* h_a, GPU_dd* h_c, const int len) {
        gpu_exp<GPU_dd>(h_a, h_c, len);
}

void gpu_qd_exp(const GPU_qd* h_a, GPU_qd* h_c, const int len) {
        gpu_exp<GPU_qd>(h_a, h_c, len);
}

/**
* c[i] = sin(a[i])
*/
void gpu_qd_sin(const GPU_qd* h_a, GPU_qd* h_c, const int len) {

#ifdef __DEVICE_EMULATION__
	printf("NOTE: EMULATION mode\n");
#endif

        //copy memory to GPU
        GPU_qd* d_a = NULL;
        GPUMALLOC((void**)&d_a, sizeof(GPU_qd)*len);
        TOGPU(d_a, h_a, sizeof(GPU_qd)*len);
        GPU_qd* d_c = NULL;
        GPUMALLOC((void**)&d_c, sizeof(GPU_qd)*len);

        //kernel
        const int numBlock = 128;
        const int numThread = 128;
        unsigned int timer = 0;
        startTimer(timer);
	if((d_sin_table == NULL) || (d_cos_table == NULL)) {
		printf("!!!sin or cos table is NULL.");
		exit(0);
	}
        map_sin_kernel<GPU_qd><<<numBlock, numThread>>>(d_a, d_c, len, 
			d_sin_table, d_cos_table);
        cutilCheckMsg("map_sin_kernel");
        CUDA_SAFE_CALL(cudaThreadSynchronize());
        endTimer(timer, "GPU sin kernel");

        //copy results from GPU
        FROMGPU(h_c, d_c, sizeof(GPU_qd)*len);

        //clean up
        GPUFREE(d_a);
        GPUFREE(d_c);
}


void gpu_dd_sin(const GPU_dd* h_a, GPU_dd* h_c, const int len) {

	
        //copy memory to GPU
        GPU_dd* d_a = NULL;
        GPUMALLOC((void**)&d_a, sizeof(GPU_dd)*len);
        TOGPU(d_a, h_a, sizeof(GPU_dd)*len);
        GPU_dd* d_c = NULL;
        GPUMALLOC((void**)&d_c, sizeof(GPU_dd)*len);

        //kernel
        const int numBlock = 128;
        const int numThread = 128;
        unsigned int timer = 0;
        startTimer(timer);
        map_sin_kernel<GPU_dd><<<numBlock, numThread>>>(d_a, d_c, len, 
                        d_dd_sin_table, d_dd_cos_table);
        cutilCheckMsg("map_sin_kernel");
        CUDA_SAFE_CALL(cudaThreadSynchronize());
        endTimer(timer, "GPU sin kernel");

        //copy results from GPU
        FROMGPU(h_c, d_c, sizeof(GPU_dd)*len);

        //clean up
        GPUFREE(d_a);
        GPUFREE(d_c);
	
}




#endif // _GQD_API_CU_
