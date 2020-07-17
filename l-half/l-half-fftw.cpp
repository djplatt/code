#include "stdio.h"
#include "stdlib.h"
#include <fftw3.h>
#include "math.h"

using namespace std;

#define H_N (10)
#define K (9)
#define fname ("L_5.dat")

double s_array[K];
double h_bernoulli[K]={1.0/6.0/2.0,-1.0/30.0/24.0,1.0/42.0/720.0,-1.0/30.0/40320.0,5.0/66.0/3628800.0,
		       -691.0/2730.0/479001600.0,7.0/6.0/87178291200.0,-3617.0/510.0/20922789888000.0,
		       43867.0/798.0/6402373705728000.0};

long unsigned int P,LOGICAL_N,PHYSICAL_N,PR;

double hurwitz_half(double alpha)
{
  int i;
  double res=0.0,s1;

  //printf("in hurwitz half with alpha=%20.18e\n",alpha);
  double n_alpha=H_N+alpha;
  double n_alpha_s=sqrt(n_alpha);
  for(s1=0.0;s1<H_N;s1++)
    res=res+1.0/sqrt(s1+alpha);
  //printf("res=%20.18e\n",res);
  res=res-2*n_alpha_s;
  //printf("res=%20.18e\n",res);

  n_alpha_s=1.0/n_alpha_s; // n_alpha_s=(n+alpha)^(-s)
		
  res=res+0.5*n_alpha_s;
  //printf("res=%20.18e\n",res);

  s_array[0]=0.5*n_alpha_s/n_alpha;
  s1=0.5;
  for(i=1;i<K;i++)
    {
      s1=s1+1;
      s_array[i]=s_array[i-1]*s1/n_alpha;
      s1=s1+1;
      s_array[i]=s_array[i]*s1/n_alpha;
    }
  
  for(i=0;i<K;i++)
    {
      res+=s_array[i]*h_bernoulli[i];
    }
  //printf("hurwitx returning %20.18e\n",res);
  return(res);
}

void populate_fft(double *in)
{
  unsigned long int pr,i,pow;
  for(i=0,pow=1;i<LOGICAL_N;i++,pow=(pow*PR)%P)
    {
      in[i]=hurwitz_half((double) pow/P);
      //printf("Hurwitz(0.5,%lu/%lu)=%20.18e\n",pow,P,in[i]);
    }
}

inline void w_cmplx(double *z, FILE *ofile)
{
  fwrite(z,sizeof(double),2,ofile);
  printf("wrote %20.18e +i%20.18e\n",z[0],z[1]);
}

int main(int argc, char ** argv )
     {       

       if(argc!=3)
	 {
	   printf("Usage l-half-fftw <n-code> <filename>\n");
	   exit(0);
	 }

       int ncode=atoi(argv[1]);
       if(ncode==1)
	 {
	   P=11;
	   LOGICAL_N=10;
	   PR=2;
	 }
       else if(ncode==34)
	 {
	   P=(((long unsigned int) 1<<34) + 25); // a large prime
	   LOGICAL_N=P-1;
	   PR=6;
	 }
       else if(ncode==0)
	 {
	   P=3;
	   for(int i=1;i<22;i++)
	     P*=3;
	   LOGICAL_N=2*(P/3);
	   PR=2;
	 }
       else if(ncode==5)
	 {
	   P=100000+3;
	   PR=2;
	   LOGICAL_N=P-1;
	 }
       else if(ncode==6)
	 {
	   P=1000000+3;
	   PR=2;
	   LOGICAL_N=P-1;
	 } 
       else if(ncode==7)
	 {
	   P=10000019;
	   PR=6;
	   LOGICAL_N=P-1;
	 }
       else if(ncode==8)
	 {
	   P=100000000+7;
	   PR=5;
	   LOGICAL_N=P-1;
	 }
       else if(ncode==32)
	 {
	   P=((unsigned long int) 1 <<32) + 15;
	   PR=3;
	   LOGICAL_N=P-1;
	 }
       else if(ncode==9)
	 {
	   P=(long unsigned int) 1000000000+7;
	   PR=5;
	   LOGICAL_N=P-1;
	 }
       else if(ncode==10)
	 {
	   P=(long unsigned int) 10000000000+19;
	   PR=2;
	   LOGICAL_N=P-1;
	 }
       else if(ncode==11)
	 {
	   P=(long unsigned int) 3;
	   for(int i=1;i<19;i++)
	     P*=3;
	   PR=2;
	   LOGICAL_N=P/3*2;
	   printf("Calculating L(1/2) for modulus=%lu with phi=%lu using pr=%lu\n",P,LOGICAL_N,PR);
	 }
       else if(ncode==12)
	 {
	   P=(long unsigned int) 3;
	   for(int i=1;i<18;i++)
	     P*=3;
	   PR=2;
	   LOGICAL_N=P/3*2;
	   printf("Calculating L(1/2) for modulus=%lu with phi=%lu using pr=%lu\n",P,LOGICAL_N,PR);
	 }
       else if(ncode==13)
	 {
	   P=(long unsigned int) 3*333333313;
	   PR=7;
	   LOGICAL_N=2*(P-1);
	   printf("Calculating L(1/2) for modulus=%lu with phi=%lu using pr=%lu\n",P,LOGICAL_N,PR);
	   printf("Not implemented yet.\n");
	   exit(0);
	 }
       else
	 {
	   printf("n-code numst be 0-14,32,34.\n");
	   exit(0);
	 }

       printf("n=%lu pr=%lu phi(n)=%lu\n",P,PR,LOGICAL_N);

       PHYSICAL_N =((LOGICAL_N>>1)+1)<<1;

       FILE *outfile;
       if(!(outfile=fopen(argv[2],"wb")))
	 {
	   printf("Error opening file %s for binary output. Exiting.\n",argv[2]);
	   exit(0);
	 }
         double *in;
         fftw_plan p;
         
	 if(!(in = (double *) fftw_malloc(sizeof(double)*PHYSICAL_N)))
	   {
	     printf("Fatal error allocating in.Exiting.\n");
	     exit(0);
	   };
	 

	 //p = fftw_plan_r2r_1d(LOGICAL_N,in,in,FFTW_R2HC, FFTW_ESTIMATE);
	 p=fftw_plan_dft_r2c_1d(LOGICAL_N,in,(fftw_complex *) in,FFTW_ESTIMATE);
	 /*
	 for(int i=0;i<LOGICAL_N;i++)
	   in[i]=i+1;
	 */
	 populate_fft(in);
	 /*
	 for(int i=0;i<PHYSICAL_N;i++)
	   printf("in[%ld]=%20.18e\n",i,in[i]);
	 */
	 fftw_execute(p);
	 /*         
	 for(int i=0;i<PHYSICAL_N;i++)
	   printf("in[%ld]=%f\n",i,in[i]);
	 */
	 fwrite(&P,sizeof(long unsigned int),1,outfile);
	 /*
	 double z[2];
	 long unsigned int i,j,k;
	 z[0]=in[0];
	 z[1]=0.0;
	 w_cmplx(z,outfile);
	 for(i=1,j=LOGICAL_N-1;i<LOGICAL_N/2-1;i++,j--)
	   {
	     z[0]=in[i];
	     z[1]=in[j];
	     w_cmplx(z,outfile);
	   }
	 z[0]=in[i];
	 z[1]=0.0;
	 w_cmplx(z,outfile);
	 
	 for(long unsigned int i=0;i<PHYSICAL_N;i+=2)
	   printf("%20.18e %20.18e\n",in[i],in[i+1]);
	 */

	 fwrite(in,sizeof(double),PHYSICAL_N,outfile);
	 fftw_destroy_plan(p);
         
         fftw_free(in);
         
	 return(0);
     };
