#include <stdio.h>
#include <fftw3.h>
#define N (8)

int main()
     {
         fftw_complex *in, *out;
         fftw_plan p;
         int ptr,lup;
	 int dims[2];
	 dims[0]=2;
	 dims[1]=4;
         in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N*2);
	 /*
	   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	 */
	 for(ptr=0;ptr<N*2;ptr++)
         {
             in[ptr][0]=ptr+1;
             in[ptr][1]=0;
         };
            p = fftw_plan_many_dft(2,dims,2,in,NULL,1,8,
				   in,NULL,1,8,FFTW_FORWARD, FFTW_ESTIMATE);

            fftw_execute(p); /* repeat as needed */
	    
	    for(ptr=0;ptr<N*2;ptr++)
	      printf("%d %e + %ei\n",ptr,in[ptr][0],in[ptr][1]);

            fftw_destroy_plan(p);
	    
         

         fftw_free(in);
         return(0);
     }
