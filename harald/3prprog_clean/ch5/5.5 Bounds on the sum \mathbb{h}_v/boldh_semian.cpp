

/* Compile with

g++ boldh_semian.cpp -oboldh_semian -O2 -frounding-math -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

This assumes that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that int_double14.2.h is in the
same directory as this file */

/*
Run with

./boldh_semian

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "int_double14.2.h"
#include "mupsigma2.h"

int_double fod1(long int p)
{
  int_double den;

  den = (p+1)-sqrt((int_double) p);
  return p/(den*den);
}

int_double fod2(long int p)
{
  int_double den, sqt;

  sqt = sqrt((int_double) p);
  den = (p+1)-sqt*sqrt(sqt);
  return (p*sqt)/(den*den);
}

int_double fod1p5(long int p)
{
  int_double den,den2,sqt,ssqt;

  sqt = sqrt((int_double) p);
  ssqt =sqrt(sqt);
  
  den = (p+1)-sqt;
  den2 = (p+1)-sqt*ssqt;
  return p*ssqt/(den*den2);
}

void fillfacblock(int_double *fac, short int *isprime,
		  int_double (*fod)(long int),
		  long int n, long int m)
/* fills fac[0], fac[1],... 
 with f(n), f(n+1), ..., f(n+m-1),
where f is the multiplicative function such that
f(p) = fod(p)  for p an odd prime
and f(n)=0 for n non-square free */
/* assumes isprime is filled and valid up to and including sqrt(n+m) */
{
  long int i;
  long int *tabla;
  long int p, red,psq;
  int_double fp;
  
  
  tabla = (long int *) calloc(m,sizeof(long int));
  p=2; fp=(*fod)(2);
  for(i=0; i<m; i++) {
    if((i+n)%2) {
      fac[i]=1; tabla[i]=1;
    }  else if((i+n)%4) {
      fac[i]=fp; tabla[i] = 2;
    } else {
      fac[i]=0.0; tabla[i]=0;
    }
  }
  
  for(p=3; p*p<=n+m; p+=2) 
    if(isprime[p]) {
      fp = (*fod)(p);
       psq=p*p;

       red = (n%psq);
       for(i= (red? psq-red : 0); i<m; i+=psq) {
	 fac[i] = 0.0; tabla[i]=0;
       }
       
       red = (n%p);
       for(i= (red? p-red : 0); i<m; i+=p)
	 if(tabla[i]) {
	   fac[i] *= fp; tabla[i] *= p;
	 }
    }

  for(i=0; i<m; i++)
    if(tabla[i] && tabla[i]!=n+i) {
      p = ((n+i)/tabla[i]);
      fp = (*fod)(p);
      fac[i] *= fp;
    }
  
  free(tabla);
}

main(int argc, char *argv[])
{
  long int NP,M,N,n,n0,v,maxn;
  short int *isprime, *sqf;
  int_double *fac, sum, L, l;
  double maxr, maxch;
  int opt, loudflag;
  
  _fpu_rndd();

  if(argc<4) {
    N=1000000; v=1; opt = 0; M=10000; loudflag=0; 
  } else {
    N=atol(argv[1]); v=atol(argv[2]); opt = atoi(argv[3]);
    if(argc>=5)
      M=v*(atol(argv[4])/v);
    else M=2*((int) -sqrt((int_double) N).right);
    /* M should be divisible by v */
    if(argc<6)
      loudflag=0;
    else loudflag = atoi(argv[5]);
  }
  NP=sqrt(M*((N-1)/M+1)+v);
  
  isprime = (short int *) calloc(NP+2,sizeof(short int));
  fac = (int_double *) calloc(M+v,sizeof(int_double));
			     
  fillisprime(isprime, NP+2);
  
  maxr=0.0; maxn=0; sum=0.0;
  
  for(n0=1; n0<=N; n0+=M) {
    if(opt==0) {
      fillfacblock(fac,isprime,&fod1,n0,M+v);
      maxch=-(1 + maxr*log((int_double) n0)).right;
    } else if(opt==1) {
      fillfacblock(fac,isprime,&fod1p5,n0,M+v);
      maxch=-(maxr*sqrt(sqrt((int_double) n0))).right;
    } else {
      fillfacblock(fac,isprime,&fod2,n0,M+v);
      maxch=-(maxr*sqrt((int_double) n0)).right;
    }
    
    for(n=n0; n<=N && n<n0+M; n+=v) {
      sum += fac[n-n0];
      if(-sum.right>maxch)
	if(n>1) {
	  if(opt==0) {
	    l = log((int_double) n);
	    if(-sum.right>(1 + maxr*l).left) {
	      maxr = -((sum-1)/l).right; 
	    }
	  } else if(opt==1) {
	    if(-sum.right>(maxr*sqrt(sqrt((int_double) n))).left)
	      maxr = -(sum/sqrt(sqrt((int_double) n))).right;
	  } else {
	    if(-sum.right>(maxr*sqrt((int_double) n)).left)
	      maxr = -(sum/sqrt((int_double) n)).right;
	  }
	}

      if(loudflag)
	if(!((n-1)%1048576))
	  if(n>1)
	    printf("%ld %g\n",n,-sum.right);
    }
  }

  if(opt==0)
    fprintf(stdout,"The bound is 1 + %.10g log y\n",maxr);
  else if(opt==1)
    fprintf(stdout,"The bound is %.10g y^(1/4)\n",maxr);
  else
    fprintf(stdout,"The bound is %.10g sqrt(y)\n",maxr);
}
