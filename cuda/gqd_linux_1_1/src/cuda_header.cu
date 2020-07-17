
#ifndef _CUDA_HEADER_CU_
#define _CUDA_HEADER_CU_

#include <stdio.h>
#include <stdlib.h>
#include <cutil_inline.h>

//#define CUDA_2

/** for CUDA 2.0 */
#ifdef CUDA_2
    #define cutilCheckMsg CUT_CHECK_ERROR
    #define cutilSafeCall CUDA_SAFE_CALL
#endif


/** macro utility */
#define GPUMALLOC(D_DATA, MEM_SIZE) cutilSafeCall(cudaMalloc(D_DATA, MEM_SIZE))
#define TOGPU(D_DATA, H_DATA, MEM_SIZE) cutilSafeCall(cudaMemcpy(D_DATA, H_DATA, MEM_SIZE, cudaMemcpyHostToDevice))
#define FROMGPU( H_DATA, D_DATA, MEM_SIZE ) cutilSafeCall(cudaMemcpy( H_DATA, D_DATA, MEM_SIZE, cudaMemcpyDeviceToHost))
#define GPUTOGPU( DEST, SRC, MEM_SIZE ) cutilSafeCall(cudaMemcpy( DEST, SRC, MEM_SIZE, cudaMemcpyDeviceToDevice ))
#define GPUFREE( MEM ) cutilSafeCall(cudaFree(MEM));


/** GPU parameters */
#define MAX_SHARED_PER_PRO_IN_BYTES (16000) //16KB
#define MAX_NUM_THREAD_PER_BLOCK (512)
#define NUM_THREADS_WARP (32) //one warp
#define NUM_THREADS_HALF_WARP (NUM_THREADS_WARP/2) //half warp

/** timer utility */
void startTimer( unsigned int & timer ) {
    	timer = 0;
	CUT_SAFE_CALL( cutCreateTimer( &timer));
	CUT_SAFE_CALL( cutStartTimer( timer));
}

void endTimer( unsigned int & timer, char* title ) {
    	CUT_SAFE_CALL( cutStopTimer( timer));
	printf( "%s processing time: %.4f sec\n", title, cutGetTimerValue( timer)/1000.0);
	CUT_SAFE_CALL( cutDeleteTimer( timer));
}

/** utility */
inline void printLine() {
	printf( "------------------\n" );
}

#endif

