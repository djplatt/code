/* Compile with

g++ checkmopt.cpp -ocheckmopt -O2 -frounding-math -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

assuming that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that int_double14.2.h is in the
same directory as this file */

/*
Run with

./checkmopt

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "int_double14.2.h"
#include "mupsigma2.h"
#include "optarith.h"
#include "mubitarr.h"

main(int argc, char *argv[])
{
  long int M,N,NP,n,v,n0,ns[3],Mer;
  short int *isprime, *mu;
  int_double kap[3], kai[3], tq, tw;
  int_double ar,ar0;
  int i,j,hund;
  int_double cent, sar, bar[3], step;
  double max0[3], maxp[3];
  int loudflag;

  ulong *imus, *isqf, *gcdind, *lor, Md, Mm, *musbit, *sqfbit;
  unsigned char *indgcd, *pdiff;
  
  _fpu_rndd();

  if(argc<3) {
    N=1000000; v=1;
    ns[0] = 2; ns[1] = 2; ns[2] = 2;
    M=10000; loudflag = 0;
  } else {
    N=atol(argv[1]); v=atol(argv[2]);
    if(argc<6) {
      if(v==1) {
	ns[0] = 2; ns[1] = 2; ns[2] = 2;
      } else {
	ns[0] = 3; ns[1] = 3; ns[2] = 3;
      }
    } else {
      ns[0]=atol(argv[3]); ns[1]=atol(argv[4]); ns[2]=atol(argv[5]);
    }
    if(argc<7)
      M=2*(-sqrt((int_double) N).right);
    else M = atoi(argv[6]);
    if(argc<8)
      loudflag=0;
    else loudflag = atoi(argv[7]);
  }

  NP=sqrt(N+v+1);
  M = ((M-1)/180180+1)*180180;
  /* 180180 = 2*3*2*3*5*7*11*13 */
    
  isprime = (short int *) calloc(NP+2,sizeof(short int));
  mu = (short int *) calloc(M+v+1,sizeof(short int));
  musbit = (ulong *) calloc(M/(8*sizeof(ulong))+1,sizeof(ulong));
  sqfbit = (ulong *) calloc(M/(8*sizeof(ulong))+1,sizeof(ulong));
  
  fillisprime(isprime, NP+2);
  initmubitdat(N,5,17,&imus,&isqf,&indgcd,&gcdind,&lor,&pdiff,&Md,&Mm);
	       
  if(v==1) {
    bar[0] = 0;
    bar[1] = 1;
    bar[2] = 2*d_gamma;
  } else {
    bar[0] = 0;
    bar[1] = 2.9;
    bar[2] = ((int_double) 577)/100;
  }

  for(j=0; j<=2; j++)
    max0[j] = maxp[j] = 0.0;
  Mer = 0;

  kap[0] = 0; kap[1] = -v; kap[2] = 2*v*(d_gamma+log((int_double) v));
  for(n0=0; n0<=N; n0+=M) {
    if(loudflag && v==1) {
      if(!((n0/M)%1000))
	printf("Mertens equals %ld at %ld\n",Mer,n0);
    }

    fillmubitblock(musbit, sqfbit, imus, isqf, indgcd, gcdind, lor,
		   pdiff, n0, M, 5, 17, Md, Mm);
    fillmufrombitarrs(mu,musbit,sqfbit,M);
      
    for(n = n0+1; n<=N && n<=n0+M; n+=v) {
      if(mu[n-n0]) {
	kap[0] += mu[n-n0]/((int_double) n);
	if(loudflag && v==1)
	  Mer+=mu[n-n0]; /* We keep track of Mertens only for debugging */
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
	kai[0] = kap[0];
	kai[1] = kap[1] + ar*kap[0];
	kai[2] = kap[2] + 2*ar*kap[1] + ar*ar*kap[0];
	/* kai[0], kai[1], kai[2] are intervals in which
	   m(x), \check{m}(x), \check{\check{m}}(x) lie
	   for n<=x<=n+v (or n+(i-1)*cent*v <= x<= n+i*cent*v, if hund>1) */

	for(j=0; j<=2; j++)
	  if(n+v>ns[j]) {
	    if(hund>1)
	      tq = kai[j]*kai[j]*(n+step);
	    else
	      tq = kai[j]*kai[j]*(n+v); /* same thing, but can save us 
					   some time and precision */
	    if(-tq.right>max0[j]) {
	      max0[j] = -tq.right;
	      if(loudflag)
		if(hund>1)
		  printf("%d %g %g\n",
			 j,-(n+step).right,-sqrt((int_double) max0[j]).right);
		else
		  printf("%d %ld %g\n",
			 j,n+v,-sqrt((int_double) max0[j]).right);
	    }
	  }
	
	for(j=1; j<=2; j++) {
	  if(hund>1)
	    tq = kai[j]*(n+step);
	  else
	    tq = kai[j]*(n+v);
	  if(-tq.right>bar[j].left) {
	    tw = tq - bar[j]; tw *= tw;
	    if(hund>1)
	      tw /= n+step;
	    else
	      tw /= n+v;
	    if(-tw.right>maxp[j]) {
	      maxp[j] = -tw.right;
	      if(loudflag)
		if(hund>1)
		  printf("+%d %g %g\n",
			 j,-(n+step).right,-sqrt((int_double) maxp[j]).right);
		else
		  printf("+%d %ld %g\n",
			 j,n+v,-sqrt((int_double) maxp[j]).right);
	    }
	  }
	}
	ar.left = -ar.right;
      }

      kap[2] += 2*ar0*kap[1] + ar0*ar0*kap[0];
      kap[1] += ar0*kap[0];
    }
  }

  if(v==1) {
    printf("Let T_0(x) = m_%ld(x),\n",v);
    printf("Let T_1(x) = \\check{m}_%ld(x) - 1,\n",v);
    printf("    T_2(x) = \\check{\\check{m}}_%ld(x) - 2\\log x + 2 \\gamma.\n",v);
  } else {
    printf("Let T_0(x) = m_%ld(x),\n",v);
    printf("Let T_1(x) = \\check{m}_%ld(x) - 2,\n",v);
    printf("    T_2(x) = \\check{\\check{m}}_%ld(x) - 4\\log(x/2) + 4 \\gamma.\n",v);
  }
  
  printf("For x\\leq %ld\n",N);
  for(j=1; j<=2; j++) {
    if(j==2 && v==1)
      printf("\t|T_%d(x)|\\leq 2 gamma/x + %.10g/sqrt(x)\n",j,
	   -sqrt((int_double) maxp[j]).right);
    else
      printf("\t|T_%d(x)|\\leq %.10g/x + %.10g/sqrt(x)\n",j,-bar[j].right,
	   -sqrt((int_double) maxp[j]).right);
  }
  
  for(j=0; j<=2; j++)
    printf("For %ld\\leq x\\leq %ld:\t |T_%d(x)|\\leq %.10g/sqrt(x) \n",
	   ns[j],N,j,
	   -sqrt((int_double) max0[j]).right);
  freemubitdat(imus,isqf,indgcd,gcdind,lor,pdiff);
}
