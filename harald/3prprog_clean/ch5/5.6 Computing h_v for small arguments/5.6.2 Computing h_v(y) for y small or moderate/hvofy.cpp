/* Compile with

g++ hvofy.cpp -ohvofy -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

This assumes that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that 
int_double14.2.h, mupsigma2.h and mdotmtildek4.h
are in the same directory as this file */

/*
Run with

./hvofy

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "int_double14.2.h"
#include "mupsigma2.h"
#include "mdotmtildek4.h"

main(int argc, char *argv[])
{
  int_double hinteg, hdotinteg, errsum, totsum;
  long int N, M, Sqt;
  int v,st;
  long i,j,n1,n0,d;
  short int *mu, *isprime;
  long int *sigma;
  int_double **mdv, **mtv;
  int_double *k0, *k1,*k2, *md2, *mt2;
  int_double k0init,k1init,k2init;
  int_double lar, larint , hint, thtint, hval, r;
  int_double third;
  double hmax, hvalabs;
  int J0=1000000;
  
  _fpu_rndd();

  if(argc<4) {
    n0=100000; N=1000000;
    v=1;   
  } else {
    n0 = atol(argv[1]);
    N = atol(argv[2]); 
    v=atoi(argv[3]);
  }
      
  Sqt = -(sqrt((int_double) N)).right; 
  if(v==1)
    M = max(Sqt,2);
  else
    M = v*(1+Sqt/v);
      
  isprime = (short int *) calloc(Sqt+1,sizeof(short int));
  mu = (short int *) calloc(Sqt+1,sizeof(short int));
  sigma = (long int *) calloc(Sqt+1,sizeof(long int));
  
  fillisprime(isprime,Sqt+1); 
  fillmublock(mu,isprime,0,Sqt+1);
  fillsigmablock(sigma,isprime,0,Sqt+1);

  fillmdmtarray(&mdv,&mtv,v,Sqt,isprime,mu,sigma);
  
  md2 = (int_double *) calloc(Sqt+1,sizeof(int_double));
  mt2 = (int_double *) calloc(Sqt+1,sizeof(int_double));
  k0 = (int_double *) calloc(M+1,sizeof(int_double));
  k1 = (int_double *) calloc(M+1,sizeof(int_double));
  k2 = (int_double *) calloc(M+1,sizeof(int_double));
  third = ((int_double) 1.0)/3.0;
  
  for(i=0, errsum = 0.0, totsum=0.0, hinteg=0.0, hdotinteg=0.0, hmax=0.0,
	k0init = k1init= k2init=0.0; i<N; i+=M) {
    n1 = min(i+M,N-1);
    fillk(k0,k1,k2,&k0init,&k1init,&k2init,
	  v,i,n1,Sqt,md2,mt2,
	  isprime,mu,sigma,mdv,mtv);

    for(j=i+1; j<=n1; j+=v) { 
      if(j+v<=N) {
	st=v;
        lar = log1p(((int_double) v)/j);
      } else {
	st=1;
	lar = log1p(((int_double) 1)/j);
      }
      /* we have just been careful not to overstep N */

      hdotinteg += k2[j-(i+1)]*lar;

      hinteg += k0[j-(i+1)]*lar + k1[j-(i+1)]*0.5*lar*lar
	+ k2[j-(i+1)]*third*lar*lar*lar;
     
      larint=lar; larint.left=0;
      /* larint is the interval [0,log(min(j+v,N)/j)] */

      if(j+(v-1)>=n0) {
	hval = k0[j-(i+1)] + larint*k1[j-(i+1)] + larint*larint*k2[j-(i+1)];
	hvalabs = -(abs(hval)*(j+v)).right;
      /* hvalabs gives an upper bound on h_v(t) t in [j,min(j+v,N)] */
	if(hvalabs>hmax)
	  hmax = hvalabs;
      }
      
      thtint = (k0[j-(i+1)]+k1[j-(i+1)]) 
        + larint*(k1[j-(i+1)]+2*k2[j-(i+1)])
        + larint*larint*k2[j-(i+1)];
      /* thtint is an interval containing all values of (t h(t))'
	 for t in the interval [j,min(j+v,N)] */

      if(thtint.left<=0 && -thtint.right>=0)
	/* if (t h(t))' might change signs in the interval */ {
	r = ((int_double) st)/j;
	totsum+=
	  abs(k0[j-(i+1)] + k1[j-(i+1)])
	  + abs(k1[j-(i+1)] + 2*k2[j-(i+1)])*
	  (j > J0 ? 0.5*st*r : j*((r+1)*lar-r))
	  + abs(k2[j-(i+1)])*
	  (j > J0 ? third*st*r*r : j*((r+1)*lar*lar-2*(r+1)*lar+2*r));
      } else
	totsum+=
	  abs(v*k0[j-(i+1)] + (j+v)*
	      (k1[j-(i+1)]*lar+k2[j-(i+1)]*lar*lar));
      
      if(j+v>=N)  /* at the last step... */
	hval = k0[j-(i+1)] + lar*k1[j-(i+1)]
          + lar*lar*k2[j-(i+1)];
      /* store the final value of h(t), namely, h(N) */
    }
  }

    printf("The integral of h_%d(t)/t from %ld to %ld lies in [%.12g %.12g]\n",v, (long) 1, N, hinteg.left,-hinteg.right);

    printf("The integral of \\dot{h}_%d(t)/t from %ld to %ld lies in [%.12g %.12g]\n",v, (long) 1, N, hdotinteg.left,-hdotinteg.right);
     printf("The maximum of h_%d(t)/t from %ld to %ld is at most %.10g\n",v, n0, N, hmax);     
  printf("The integral of |(t h_%d)'| from %ld to %ld is at most: [%.10g %.10g]\n",v, (long) 1, N, totsum.left,-totsum.right);
  /* It would be enough to output -totsum.right, but we output an interval
     to get a sense of the loss of accuracy */
  
  printf("If we add |%ld h_%d(%ld)| to that, we obtain at most: %.10g\n",N,v,N,-(totsum+abs(hval)*N).right);
}
