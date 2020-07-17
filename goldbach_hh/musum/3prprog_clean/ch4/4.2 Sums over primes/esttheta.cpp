#include <primesieve.h>
#include <stdio.h>
#include "int_double14.2.h"

/*
Compile with

g++ esttheta.cpp -oesttheta -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm -lprimesieve

*/

int main(int argc, char *argv[])
{
  uint64_t start; 
  uint64_t stop;
  size_t i;
  size_t size;
  int_double theta, pr;
  double mini, maxi;
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
      pr = primes[i];
      if(primes[i]>=LB && theta.left < -(pr*mini).right)
	mini = (theta/pr).left;
      theta+=log(pr);
      if(primes[i]>=LB && -theta.right > (pr*maxi).left)
	maxi = -(theta/pr).right;
    }
    primesieve_free(primes);
  }
  if(theta.left<-(((int_double) N)*mini).right)
    mini = (theta/N).left;
  if(-theta.right>(((int_double) N)*maxi).left)
    maxi = -(theta/N).right;
  
  printf("theta(%ld) lies in [%.8g %.8g]\n",N,theta.left,-theta.right);
  printf("min of theta(x)/x for %ld<=x<=%ld is at least %.8g\n",LB,N,
	 mini);
  printf("max of theta(x)/x for %ld<=x<=%ld is at most  %.8g\n",LB,N,
	 maxi);
}
