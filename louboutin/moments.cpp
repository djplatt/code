//
// moments.cpp
//
// Compute S_{2k}^-(p) k=1..9
// p is any odd prime
// p is specified as a power of 2 (>=2, <=40), an offset
// and a primitive root, e.g. ./moments 3 3 2 will compute for p=11
//
// The algorithm works by writing
//
// \theta(1,\chi)=\sum\limits_{n=1}^{p-1} \chi(n) \sum\limits_{m\geq 0} (n+mp)\exp(-\pi (n+mp)^2/p)
//
// which we can efficiently evaluate by a DFT.
// We use FFTW to perform the DFT
//
#include "stdio.h"
#include "stdlib.h"
#include <fftw3.h>
#include "math.h"
#include "inttypes.h"

using namespace std;

// P is our modulus (should have a primitive root)
// LOGICAL_N=phi(P) is length of our DFT
// PHYSICAL_N is the length needed by FFTW
// PR is our primitive root

uint64_t P,LOGICAL_N,PHYSICAL_N,PR;

// populate the DFT with n*exp(-Pi*n*n/P)
// until this is very small
void populate_fft(double *in)
{
  uint64_t i,pow;
  for(i=0;i<PHYSICAL_N;i++) // zero the DFT vector
    in[i]=0.0;
  for(i=0,pow=1;i<P-1;i++,pow=(pow*PR)%P) // i indexes the vector
                                          // pow runs through powers of the
                                          // primitive root
    {
      double x=pow*exp(-M_PI*pow*pow/P);
      if(x<1e-100) // very small
	continue;
      uint64_t n=pow; // first n in theta
      while(x>1e-100)
	{
	  in[i]=in[i]+x;
	  n+=P;                 // next n
	  x=n*exp(-M_PI*n*n/P);
	}
    }
}

// the conjectured dependance on p
// p^(1+3k/2) log^((k-1)^2) (p)
double conjecture(uint64_t p, uint64_t k)
{
  double logp=log(p);
  return(exp((1.0+1.5*k)*logp+(k-1)*(k-1)*log(logp)));
}

// z is already log'd
double term (double z, uint64_t k)
{
  return(exp(k*z));
}

// sum up the odd characters
double sum_odd(fftw_complex *out, uint64_t k)
{
  double res=0;
  uint64_t i;
  for(i=1;i<(P-1)/2;i+=2)
    res+=2*term(out[i][0],k); // 2* because of conjugates
  if((P&3)==3) // P=3(4) so add (P-1)/2 (real) term once
    res+=term(out[i][0],k);
  return res;
}

// replace the real part of all the entries corresponding to odd characters
// with log|z|
void mod_entries(fftw_complex *out)
{
  uint64_t i;
  for(i=1;i<(P-1)/2;i+=2) // the complex results
    out[i][0]=log(out[i][0]*out[i][0]+out[i][1]*out[i][1]);
  if((P&3)==3) // it has a real result at the end
    out[i][0]=log(out[i][0]*out[i][0]); // so log it
}

int main(int argc, char ** argv )
{       
  // check command line
  if(argc!=4)
    {
      printf("Usage:- %s <power of 2> <offset> <pr>\n",argv[0]);
      exit(0);
    }

  int pow=atoi(argv[1]);
  if((pow<2)||(pow>40))
    {
      printf("We expect <pow> to be in [2,40]. Exiting.\n");
      exit(0);
    }

  // compute P as 2^pow+offset
  P=1;
  P<<=pow;
  P+=atoi(argv[2]);
  LOGICAL_N=P-1; // assume P is prime for now
  PR=atol(argv[3]); // primitive root is given, could compute it ourselves

  PHYSICAL_N =((LOGICAL_N>>1)+1)<<1; // see FFTW documentation on r2c

  printf("n=%lu pr=%lu phi(n)=%lu fft_len=%lu\n",P,PR,LOGICAL_N,PHYSICAL_N);

  double *in;

  fftw_plan p;
         
  if(!(in = (double *) fftw_malloc(sizeof(double)*PHYSICAL_N)))
    {
      printf("Fatal error allocating %lu bytes for in. Exiting.\n",sizeof(double)*PHYSICAL_N);
      exit(0);
    };


  fftw_complex *out=(fftw_complex *)in; // same memory, FFT is in place
	 
  // we are doing a 1-dimensional real to complex DFT
  p=fftw_plan_dft_r2c_1d(LOGICAL_N,in,out,FFTW_ESTIMATE);

  // put the n*\exp(-\pi n^2/p) values in the vector in the right places
  populate_fft(in);

  // do the DFT
  fftw_execute(p);

  // convert to log(|z|)
  mod_entries(out);

  // now compute S_{2k}^-
  uint64_t k;
  for(k=1;k<10;k++)
    {
      double sumo=sum_odd(out,k);
      double c=conjecture(P,k);
      printf("p= %lu k= %d c= %30.28e sum= %30.28e sum/c= %30.28e\n",P,k,c,sumo,sumo/c);
    }

  // get rid of our plan
  fftw_destroy_plan(p); 

  // give the memory back 
  fftw_free(in);
         
  // success
  return(0);
}
