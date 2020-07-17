// nvcc -arch=compute_20
//
#include "stdio.h"
// Device code
#define N (10)

typedef struct{
double left;
double right;
} int_double;

__global__ void VecAdd(int_double* A, int_double* B, int_double* C)
{
int i = threadIdx.x;
if (i < N)
C[i].left = __ddiv_rd(A[i].left,B[i].right);
C[i].right = __ddiv_ru(A[i].right,B[i].left);
}

// Host code
int main()
{
int_double h_A[N],h_B[N],h_C[N];
for(int i=0;i<N;i++)
{
h_A[i].left=i;h_A[i].right=i;
h_B[i].left=3.0;h_B[i].right=3.01;
}
// Allocate vectors in device memory
size_t size = N * sizeof(int_double);
int_double* d_A;
cudaMalloc((void**)&d_A, size);
int_double* d_B;
cudaMalloc((void**)&d_B, size);
int_double* d_C;
cudaMalloc((void**)&d_C, size);
// Copy vectors from host memory to device memory
// h_A and h_B are input vectors stored in host memory
cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);
// Invoke kernel
int threadsPerBlock = 256;
int threadsPerGrid = (N + threadsPerBlock -1);
VecAdd<<<threadsPerGrid, threadsPerBlock>>>(d_A, d_B, d_C);
// Copy result from device memory to host memory
// h_C contains the result in host memory
cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);
// Free device memory
cudaFree(d_A);
cudaFree(d_B);
cudaFree(d_C);
for(int i=0;i<N;i++) printf("h_C[%d]=[%20.18e,%20.18e]\n",i,h_C[i].left,h_C[i].right);
}
