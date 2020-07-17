// nvcc -arch=compute_20
//
#include "stdio.h"
#include "inttypes.h"
#include "time.h"
#include "math.h"
// Device code

#define TPB (256) // number of threads per block
#define MAX_V (200000)

typedef struct{
  double _sqrt;
  double _log;
} table_t;


__global__ void build_table(table_t *d_table)
{
  int thread = blockIdx.x*TPB+threadIdx.x;
  if ((thread==0)||(thread>MAX_V))
    return;
  double n=(double) thread;;
  d_table[thread]._sqrt=1.0/sqrt(n);
  d_table[thread]._log=log(n);
}


#define NUM_T_POWERS (10) // how many powers of t
#define NUM_T2PI_POWERS (10) // how many powers of t/2Pi
__global__ void RS(const double *ts, const long unsigned int *N, double *res, const table_t *table)
{
  int thread = blockIdx.x*TPB+threadIdx.x;
  if (thread >=N[0])
    return;

  double t=ts[thread];

  // set up powers of t
  double t_pows[NUM_T_POWERS];
  t_pows[0]=t;
  for(int j=1;j<NUM_T_POWERS;j++)
    t_pows[j]=t_pows[j-1]*t;

  //set up powers of t/2pi
  double t2pi_pows[NUM_T2PI_POWERS];
  t2pi_pows[0]=t/(2.0*M_PI);
  for(int j=1;j<NUM_T2PI_POWERS;j++)
    t2pi_pows[j]=t2pi_pows[j-1]*t2pi_pows[0];

  // V=floor(sqrt(t/2pi))
  double sqrt_t_2pi=sqrt(t2pi_pows[0]);
  unsigned long int V=floor(sqrt_t_2pi);
  if(V>MAX_V)
    return;
  double theta=t/2.0*(log(t2pi_pows[0])-1)-M_PI/8.0+1.0/(48.0*t);
  theta=theta+7.0/(5760.0*t_pows[2])+31.0/(80640.0*t_pows[4]);
  double res1=0.0;
  for(unsigned long int n=1;n<=V;n++)
    res1=res1+table[n]._sqrt*cos(theta-t*table[n]._log);
  res1=res1+res1;

  double p=sqrt_t_2pi-V;
  double C0=cos(2.0*M_PI*(p*p-p-1.0/16.0))/cos(2.0*M_PI*p);
  if(V&1)
    res1=res1+1.0/sqrt(sqrt_t_2pi)*C0;
  else
    res1=res1-1.0/sqrt(sqrt_t_2pi)*C0;
  res[thread]=res1;
}

#define _N_ (8192) // how many zeros to send to GPU at one time
// Host code


// read our 13 byte structure representing a zero gap
// into a quad
double in_bytes(FILE *infile)
{
  uint64_t a;uint32_t b;uint8_t c;
  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  double r1;
  //printf("Read a=%lX b=%x c=%x\n",a,b,c);
  r1=(double) c/32.0 + (double) b/(32.0*256.0*256.0*256.0*256.0)+(double) (a&0xfff8000000000000)/(32.0*65536.0*65536.0*65536.0*65536.0*65536.0*65536.0); // /2^3, /2^37, /2^101
  //r2=(double) (a&0x0007ffffffffffff)/(32.0*65536.0*65536.0*65536.0*65536.0*65536.0*65536.0); // /2^101
  //printf("r1=%50.48e\nr2=%50.48e\n",r1,r2);

  return(r1);
}

void print_elapsed(clock_t st_time)
{
  clock_t this_time=clock();
  //printf("this_time=%lu start_time %lu\n",this_time,st_time);
  printf("%f seconds elapsed.\n",(double) (this_time-st_time)/(double) CLOCKS_PER_SEC);
}

//
int main(int argc, char **argv)
{
  clock_t st_time=clock();
  printf("In main\n");
  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n");
      exit(0);
    }
  int devcount;
  
  cudaGetDeviceCount(&devcount);
  printf("Found %d CUDA devices.\n",devcount);
  cudaSetDevice(1);
  //cudaDeviceReset();
  
  printf("Allocating device memory for sin_table\n");
  print_elapsed(st_time);
  double h_A[_N_],h_B[_N_];
  
  table_t *d_table;
  cudaMalloc((void **)&d_table,sizeof(table_t)*(MAX_V+1));


  long unsigned int h_N[1];
  h_N[0]=_N_;
  //for(int i=0;i<N;i++) printf("%20.18e\n",to_double(h_A[i]));

// Allocate vectors in device memory

  printf("Allocating other device memory.\n");
  print_elapsed(st_time);
  size_t size = h_N[0] * sizeof(double);
  double* d_A;
  cudaMalloc((void**)&d_A, size);
  double* d_B;
  cudaMalloc((void**)&d_B, size);
  unsigned long int *d_N;
  cudaMalloc((void **)&d_N,sizeof(long unsigned int));
  
  // Copy vectors from host memory to device memory
  // h_A and h_B are input vectors stored in host memory

  int threadsPerBlock = TPB;
  int num_blocks = MAX_V/TPB;
  if(num_blocks*threadsPerBlock<MAX_V)
    num_blocks++;


  build_table<<<num_blocks, threadsPerBlock>>>(d_table);

  num_blocks=h_N[0]/threadsPerBlock;
  if(num_blocks*threadsPerBlock<h_N[0])
    num_blocks++;


  long int num_its;
  fread(&num_its,sizeof(long int),1,infile);
  printf("Doing %ld iterations on file %s.\n",num_its,argv[1]);
  double st[2];
  long int zs[2];
  unsigned long int ptr=0;
  double del_t,t;
  for(long int it=0;it<num_its;it++)
    {
      if(it==100) break;
      
      fread(st,sizeof(double),2,infile);
      fread(&zs[0],sizeof(long int),1,infile);
      if(st[0]==0.0)
	continue;
      t=st[0];
      fread(&zs[1],sizeof(long int),1,infile);
      //printf("Processing zero %ld to %ld=%ld in total.\n",zs[0]+1,zs[1],zs[1]-zs[0]);
      for(long int z=zs[0]+1;z<=zs[1];z++)
	{
	  del_t=in_bytes(infile);
	  if(del_t==0.0)
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  t+=del_t;
	  //printf("Zero at %80.78e\n",to_double(t));
	  h_A[ptr++]=t;
	  if(ptr==h_N[0])
	   {
	     //printf("Invoking GPU with first zero at %f\n",to_double(h_A[0]));
	     //printf("Last zero at %f\n",to_double(h_A[ptr-1]));
	     //print_elapsed(st_time);

	     // copy zeros and N to device memory	     
	     cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
	     cudaMemcpy(d_N, h_N, sizeof(long unsigned int), cudaMemcpyHostToDevice);
	     // Invoke kernel
	     RS<<<num_blocks, threadsPerBlock>>>(d_A, d_N, d_B,d_table);

	     // Copy result from device memory d_B to host memory h_B
	     cudaMemcpy(h_B, d_B, size, cudaMemcpyDeviceToHost);
	     //printf("First result was %f\n",to_double(h_B[0]));
	     //printf("Last result was %f\n",to_double(h_B[ptr-1]));
	     //print_elapsed(st_time);

	     ptr=0;
	   }
       }
   }
 if(ptr!=0) // zeros didn't fit exactly into blocks
   {
     //printf("Invoking GPU with first zero at %f\n",to_double(h_A[0]));
     //printf("Last zero at %f\n",to_double(h_A[ptr-1]));
     
     cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
     cudaMemcpy(d_N, &ptr, sizeof(long unsigned int), cudaMemcpyHostToDevice);
     // Invoke kernel
     
     RS<<<num_blocks, threadsPerBlock>>>(d_A, d_N, d_B,d_table);
     
     // Copy result from device memory to host memory
     // h_B contains the result in host memory
     cudaMemcpy(h_B, d_B, size, cudaMemcpyDeviceToHost);
     //printf("First result was %f\n",to_double(h_B[0]));
     //printf("Last result was %f\n",to_double(h_B[ptr-1]));
   }
 for(long unsigned int i=(ptr>=20?ptr-20:0);i<ptr;i++)
   printf("Zero at %f rs returned %f\n",h_A[i],h_B[i]);
 print_elapsed(st_time);


  // Free device memory
  cudaFree(d_A);
  cudaFree(d_B);
  cudaFree(d_N);
}
