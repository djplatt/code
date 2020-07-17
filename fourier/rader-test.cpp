#include <stdio.h>
#include <fftw3.h>
#define N (49919)
#define N_LOOPS 1

int main()
     {
         fftw_complex *in, *out;
         fftw_plan p;
         int ptr,lup;
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	 /*
	   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	 */
	 for(ptr=0;ptr<N;ptr++)
         {
             in[ptr][0]=1.0/(ptr+1);
             in[ptr][1]=0;
         };
         for(lup=0;lup<N_LOOPS;lup++)
         {
            p = fftw_plan_dft_1d(N, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
            fftw_execute(p); /* repeat as needed */
         
	    for(ptr=1;ptr<N;ptr=ptr<<1)
	      printf("%d %e + %ei\n",ptr,in[ptr][0],in[ptr][1]);

            fftw_destroy_plan(p);
         
            p = fftw_plan_dft_1d(N, in, in, FFTW_BACKWARD, FFTW_ESTIMATE);
            fftw_execute(p);
            for(ptr=0;ptr<N;ptr++)
            {
               in[ptr][0]=in[ptr][0]/N;
               in[ptr][1]=in[ptr][1]/N;
            };
            fftw_destroy_plan(p);
         };
	 for(ptr=1;ptr<N;ptr=ptr<<1)
             printf("%d %e + %ei\n",ptr,in[ptr][0],in[ptr][1]);
         fftw_free(in);
         return(0);
     };
