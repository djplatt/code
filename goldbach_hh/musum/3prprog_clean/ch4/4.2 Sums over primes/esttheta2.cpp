#include <primesieve.h>
#include <stdio.h>
#include "int_double14.2.h"

/*
Compile with

g++ esttheta2.cpp -oesttheta2 -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lprimesieve

*/

int main(int argc, char *argv[])
{
  uint64_t start; 
  uint64_t stop;
  size_t i;
  size_t size;
  int_double theta, pr, sqt;
  double mini,maxi;
  long N,M,LB;
  long *primes;
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


  theta=0.0; mini=1.0; maxi=-1.0;
  for(start=1, stop=M; start<=N; start+=M,
	stop = (stop+M <= N ? stop+M : N)) {
   if(loudflag)
     fprintf(stderr,"%lu %g %g\n",start,mini,maxi);
    primes = (long*) primesieve_generate_primes(start, stop, &size, LONG_PRIMES);
    for (i = 0; i < size; i++) {
      pr = primes[i]; sqt = sqrt(pr);
      if(primes[i]>=LB && (theta-pr).left < -(mini*sqt).right)
	mini = ((theta-pr)/sqt).left;
      theta+=log(pr);
      if(primes[i]>=LB && -(theta-pr).right > (maxi*sqt).left)
	maxi = ((theta-pr)/sqt).left;
    }
    primesieve_free(primes);
  }
  pr=N; sqt=sqrt(pr);
  if(N>=LB && (theta-pr).left < -(mini*sqt).right)
    mini = ((theta-pr)/sqt).left;  
  if(N>=LB && -(theta-pr).right > (maxi*sqt).left)
    maxi = ((theta-pr)/sqt).left;

  printf("theta(%ld) lies in [%.10g %.10g]\n",N,theta.left,-theta.right);
  printf("min of (theta(x)-x)/sqrt(x) for %ld<=x<=%ld is at least %.10g\n",LB,N,
	 mini);
  printf("max of (theta(x)/x)/sqrt(x) for %ld<=x<=%ld is at most  %.10g\n",LB,N,
	 maxi);
}
