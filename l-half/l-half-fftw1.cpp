#include "stdio.h"
#include "stdlib.h"
#include <fftw3.h>

using namespace std;

int main()
     {
         double *in;
         fftw_plan p;
	 unsigned long int P=9;
	 unsigned long int LOGICAL_N=P-1;
	 unsigned long int PHYSICAL_N=((LOGICAL_N>>1)+1)<<1;
         
	 printf("Trying Q=%ld\n",P);
	 if(!(in = (double *) fftw_malloc(sizeof(double)*PHYSICAL_N)))
	   {
	     printf("Fatal error allocating in.Exiting.\n");
	     exit(0);
	   };
	 

	 p = fftw_plan_r2r_1d(LOGICAL_N,in,in,FFTW_R2HC, FFTW_ESTIMATE);
	 
	 for(int i=0;i<LOGICAL_N;i++)
	   in[i]=i*i;
	 
	 fftw_execute(p); /* repeat as needed */
         
	 for(int i=0;i<PHYSICAL_N;i++)
	   printf("in[%ld]=%f\n",i,in[i]);
	 
	 fftw_destroy_plan(p);
         
         fftw_free(in);
         
	 return(0);
     };
