#include <stdio.h>
#include <fftw3.h>
#include <math.h>
#define pi (3.141592535898)

fftw_complex *gauss_sum(unsigned int G, unsigned int N)
     {
         fftw_complex *in;
         fftw_plan p;
         unsigned int ptr,pow;
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N-1));
	 
	 
	 pow=1;
	 for(ptr=0;ptr<N-1;ptr++)
         {
	   if((pow==1)&ptr)
	     {
	       printf("%d is not a primitive root of %d",G,N);
	       return(NULL);
	     };
             in[ptr][0]=cos(2*pi*pow/N);
             in[ptr][1]=sin(2*pi*pow/N);
	     pow=(pow*G)%N;
         };


	 for(ptr=0;ptr<N-1;ptr=ptr++)
             printf("%d %e + %ei\n",ptr,in[ptr][0],in[ptr][1]);

	 p = fftw_plan_dft_1d(N-1, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
	 fftw_execute(p); /* repeat as needed */
         
	 fftw_destroy_plan(p);
	 printf("\n");
	 for(ptr=0;ptr<N-1;ptr++)
             printf("%d %e + %ei\n",ptr,in[ptr][0],in[ptr][1]);
         return(in);
     }


int main()
{
  gauss_sum(2,5);
}
