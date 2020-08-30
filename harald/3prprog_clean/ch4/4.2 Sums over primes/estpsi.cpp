#include <primesieve.h>
#include <stdio.h>
#include "int_double14.2.h"

/*
Compile with

g++ estpsi.cpp -oestpsi -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lprimesieve

*/

int main(int argc, char *argv[])
{
  uint64_t start;
  uint64_t stop;
  size_t i;
  long n, j, k ,puis , por;
  size_t size;
  int_double psi, pr;
  double mini,maxi;
  long N,M,LB,n0,n1,oldn;
  long *primes;
  long *mangy;
  int loudflag;
  
  _fpu_rndd();
     
  if(argc<3) {
    LB=100000;
    N=1000000;
    M=10000;
    loudflag=0;
  } else if(argc<4) {
    LB=atol(argv[1]);
    N=atol(argv[2]);
    M=2*(-sqrt((int_double) N).right);
    loudflag=0;
  } else {
    LB=atol(argv[1]);
    N=atol(argv[2]);
    M=atol(argv[3]);
    if(argc<5)
      loudflag=0;
    else loudflag = atoi(argv[4]);
  }

  psi=0.0; mini=1.0; maxi=-1.0; oldn=0;
  for(n0=1, n1=M; n0<=N; n0+=M, n1 = min(n1+M,N)) {
    if(loudflag)
      fprintf(stderr,"%ld %g %g\n",n0,mini,maxi);
    mangy = (long *) calloc(M,sizeof(long));

    start=n0; stop=n1;
    primes = (long*)
      primesieve_generate_primes(start, stop, &size, LONG_PRIMES);
    for (i = 0; i < size; i++) 
      mangy[primes[i]-n0] = primes[i];
    primesieve_free(primes);
    
    for(j=2; j<=-(log((int_double) n1)/log((int_double) 2.0)).right; j++) {
      start = pow((int_double) n0, ((int_double) 1.0)/j).left;
      stop = -pow((int_double) n1, ((int_double) 1.0)/j).right;
      primes = (long*)
	primesieve_generate_primes(start, stop, &size, LONG_PRIMES);

      for (i = 0; i < size; i++) {
	puis = primes[i];
	for(k=2; k<=j; k++)
	  puis *= primes[i];

	if(puis>=n0 && puis<=n1)
	  mangy[puis-n0] = primes[i];
      }
      primesieve_free(primes);
    }

    for(n=n0, j=0; n<=n1; n++, j++)
      if(mangy[j]) {
	psi -= (n-oldn);
	if(n>=LB &&psi.left < -(n*((int_double) mini)).right)
	  mini = (psi/n).left;
	psi+=log((int_double) (mangy[j]));
	if(n>=LB && ((-psi.right) > (n*((int_double) maxi)).left)) {
	  maxi = -(psi/n).right;
	}
	oldn = n;
      }
    free(mangy);
  }

  psi -= (N-oldn);
  
  if(psi.left<-(((int_double) N)*mini).right)
    mini = (psi/N).left;

  printf("psi(%ld)-%ld lies in [%.14g %.14g]\n",N,N,psi.left,-psi.right);
  printf("min of (psi(x)-x)/x for %ld<=x<=%ld is at least %.14g\n",LB,N,
         mini);
  printf("max of (psi(x)-x)/x for %ld<=x<=%ld is at most  %.14g\n",LB,N,
         maxi);
}

