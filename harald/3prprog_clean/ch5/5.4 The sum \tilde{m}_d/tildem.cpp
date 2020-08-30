#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "int_double14.2.h"
#include "mupsigma2.h"
#include "optarith.h"

/* Compile with

g++ tildem.cpp -otildem -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

assuming that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that int_double14.2.h is in the
same directory as this file */


/*
Run with

./tildem

*/

main(int argc, char *argv[])
{
  long int M, N, NP, n, v, n0, Mer;
  short int *isprime, *mu;
  long int *sigma;
  int_double mdot, mtilde;
  int i, hund;
  int_double ar0,ar,sar,mt,c1,c1sq,c1qt,step,tw,tq,Q,cent;
  double max0,maxp;
  int loudflag;
  
  _fpu_rndd();

  
  if(argc<3) {
    N=1000000; v=1; M=10000; loudflag = 0;
  } else {
    N=atol(argv[1]); v=atol(argv[2]);
    if(argc<4)
      M=2*((int) -sqrt((int_double) N).right);
    else M = v*(atoi(argv[3])/v);
    if(argc<5)
      loudflag=0;
    else loudflag = atoi(argv[4]);
  }

  NP=sqrt(N+v+1);
    
  isprime = (short int *) calloc(NP+2,sizeof(short int));
  mu = (short int *) calloc(M+v,sizeof(short int));
  sigma = (long int *) calloc(M+v,sizeof(long int));
	
  fillisprime(isprime, NP+2);
  
  mdot=0.0; Mer=0;
  if(v==1) {
    mtilde = (-d_pi*d_pi/6.0); c1 = 1.8;
  } else {
    mtilde = (-d_pi*d_pi/4.0); c1 = 3.9877;
  }
  c1sq = c1*c1; c1qt = c1sq*c1sq; max0=maxp=0;

  for(n0=0; n0<=N; n0+=M) {
    if(loudflag && v==1) {
      if(!((n0/M)%1000))
        printf("Mertens equals %ld at %ld\n",Mer,n0);
    }
    fillmublock(mu,isprime,n0,M+1);
    fillsigmablock(sigma,isprime,n0,M+1);
    for(n=n0+1; n<=N && n<=n0+M; n+=v) {
      if(mu[n-n0]) {
	mdot += ((int_double) mu[n-n0])/sigma[n-n0]; 
	if(loudflag && v==1)
	  Mer += mu[n-n0]; /* We keep track of Mertens only for debugging */
      }
      if(n>1000000) {
        hund = 1;
        cent = 1.0;
      } else {
        hund = 1000000/n;
        cent = ((int_double) 1.0)/hund;
      }

      ar.left=0;
      if(n>100000) 
	ar0 = log1p_small(((int_double) v)/n);
      else
	ar0 = log1p(((int_double) v)/n);

      for(i=1; i<=hund; i++) {
        if(hund>1) {
          step = i*cent*v;
          sar = log1p(step/n);
        } else {
          step = v;
          sar = ar0;
        }
        
        ar.right = sar.right;
	mt = mtilde + ar*mdot;
        /* mt is an interval in which \check{m}(x) lies
           for n<=x<=n+v (or n+(i-1)*cent*v <= x<= n+i*cent*v, if hund>1) */
	
	if(n>1) {
	  if(hund>1)
	    tq = mt*mt*(n+step);
	  else
	    tq = mt*mt*(n+v); /* same thing, but can save us 
                                           some time and precision */
	  if(-tq.right>max0) {
	    max0 = -tq.right;
	    if(loudflag)
	      if(hund>1)
		printf("%g %g\n",
		       -(n+step).right,-sqrt((int_double) max0).right);
	      else
		printf("%ld %g\n",
		       n+v,-sqrt((int_double) max0).right);
	  }

	  if(hund>1)
	    tq = tq*tq*(n+step);
	  else
	    tq = tq*tq*(n+v); 
	  if(-tq.right>c1qt.left) {
	    Q = sqrt(n+step);
	    tw = abs(mt)*Q - c1/sqrt(Q);
	   
            if(-tw.right>maxp) {
              maxp = -tw.right;
              if(loudflag)
                if(hund>1)
                  printf("+ %g %g\n",-(n+step).right,maxp);
                else
                  printf("+ %ld %g\n",n+v,maxp);
            }
          }
        }
        ar.left = -ar.right;
      }

      mtilde += ar0*mdot;
    }
  }

  printf("Let tildem_v(x) = \\sum_{n\\leq x: (n,v)=1} mu(n) log(x/n)/sigma(n)\n"); 
  printf("For x\\leq %ld\n",N);
  printf("\t|tildem_%ld(x) - pi^2/%d|\\leq %.10g/sqrt(x)\n",v,(v==1 ? 6 : 4),
	 -sqrt((int_double) max0).right);
  printf("\t|tildem_%ld(x) - pi^2/%d|\\leq %.10g/x^(3/4) + %.10g/sqrt(x)\n",
	 v, (v==1 ? 6 : 4), -c1.right, maxp);
}
