// nvcc -arch=compute_20
//
#include "stdio.h"
#include "gqd.cu"
// Device code
#define N (1024)

__global__ void VecSqr(GPU_qd* A, GPU_qd* B)
{
  int i = threadIdx.x;
  if (i < N)
    B[i]=A[i]*A[i];
}


// Host code
int main()
{
  //printf("In main\n");
  GPU_qd h_A[N],h_B[N];
  for(int i=0;i<N;i++)
    {
      h_A[i]=make_qd((double) i,0.0,0.0,0.0);
    }
  //for(int i=0;i<N;i++) printf("%20.18e\n",to_double(h_A[i]));

// Allocate vectors in device memory
  size_t size = N * sizeof(GPU_qd);
  GPU_qd* d_A;
  cudaMalloc((void**)&d_A, size);
  GPU_qd* d_B;
  cudaMalloc((void**)&d_B, size);
  // Copy vectors from host memory to device memory
  // h_A and h_B are input vectors stored in host memory
  cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
  // Invoke kernel
  int threadsPerBlock = 256;
  int threadsPerGrid = (N + threadsPerBlock -1);
  VecSqr<<<threadsPerGrid, threadsPerBlock>>>(d_A, d_B);
  // Copy result from device memory to host memory
  // h_C contains the result in host memory
  cudaMemcpy(h_B, d_B, size, cudaMemcpyDeviceToHost);
  // Free device memory
  cudaFree(d_A);
  cudaFree(d_B);
  for(int i=0;i<20;i++) printf("%20.18e\n",to_double(h_B[i]));
}
