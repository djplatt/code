/* Computes hbold(y,y) brutishly */

/* Compile with

g++ hboldvofy.cpp -ohboldvofy -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

This assumes that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that int_double12.0.h is in the
same directory as this file */

/*
Run with

./hboldvofy

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "int_double14.2.h"
#include "mupsigma2.h"
#include "mdotmtildek4.h"

int_double zeta2, zeta4;

void compkcircb(int_double *kc0b, int_double *kc1b, int_double *kc2b,
	    int v, long n, short *mu, long *sigma,
	    int_double **mdv, int_double **mtv)
/* determines g_v(y,y) for n<=y<=n+v:
   g_v(y,y) = *kc0 + (*kc1) log(y/n) + (*kc2) (log(y/n))^2 */
{
  long d;
  int_double hf, lf, isds, sd;

  *kc0b = *kc1b = *kc2b = 0.0;
  for(d=1; d<=n; d+=v)
    if(mu[d]) {
      sd = sigma[d];
      isds = ((int_double) 1.0)/(sd*sd);
      hf = mdv[d][n/d];
      lf = mtv[d][n/d]+(zeta2*sigma[d]*sigma[v])/(d*v)
	+mdv[d][n/d]*(log((int_double) n)-log((int_double) (d*(n/d))));
      *kc0b += lf*lf*isds;
      *kc1b += 2*hf*lf*isds;
      *kc2b += hf*hf*isds;
    }
}

main(int argc, char *argv[])
{
  long int N, Sqt;
  int v;
  long int i, j, n;
  int_double kc0b, kc1b, kc2b, k0b, k1b, k2b, wabs, x;
  double y, cmax, clmax;
  short int *mu, *isprime;
  long int *sigma, *phi;
  int_double **mdv, **mtv;
  int_double lar, ctm, ccf;
  int_double tot0, tot1, c1;
  int iter, ittil=10000, ittil2=1000000;
  int_double step, stx, lst, lint;
  int flag;
  
  _fpu_rndd();

  if(argc<3) {
    N=100000; v=1; flag = 0;
  } else {
    N=atol(argv[1]); v=atoi(argv[2]);
    if(argc<4)
      flag=1;
    else
      flag = atoi(argv[3]);
  }

  zeta2 = d_pi*d_pi/6.0;
  zeta4 = d_pi*d_pi*d_pi*d_pi/90.0;

  if(v==1) {
    ctm = zeta2*zeta2*zeta2/zeta4;
    ccf = 2.0*zeta2;
  } else {
    ctm = (((int_double) 9.0)/5.0)*zeta2*zeta2*zeta2/zeta4;
    ccf = 3.0*zeta2;
  }
  /* ctm = (sigma(v)/v)^2 * zeta(2)^3/zeta(4) 
     ccf =  2 (sigma(v)/v) zeta(2) */

  if(v==1) 
    c1 = ((int_double) 858)/1000.0;
  else 
    c1 = ((int_double) 175)/100.0;

  Sqt = sqrt(N);
  isprime = (short int *) calloc(Sqt+1,sizeof(short int));
  mu = (short int *) calloc(N+1,sizeof(short int));
  sigma = (long int *) calloc(N+1,sizeof(long int));
  phi = (long int *) calloc(N+1,sizeof(long int));
    
  fillisprime(isprime,Sqt+1); 
  fillmublock(mu,isprime,0,N+1);
  fillsigmablock(sigma,isprime,0,N+1);
  fillphiblock(phi,isprime,0,N+1);
  
  mdv=(int_double **) calloc(N+1,sizeof(int_double *));
  mtv=(int_double **) calloc(N+1,sizeof(int_double *));
  for(i=1; i<=N; i++) {
    mdv[i] = (int_double *) calloc(N/i+1,sizeof(int_double));
    mtv[i] = (int_double *) calloc(N/i+1,sizeof(int_double));
    fillmdot(mdv[i],isprime,v*i,N/i,mu,sigma);
    fillmtilde(mtv[i],mdv[i],v*i,N/i,mu,sigma,sigma[i]*sigma[v]);
  }

  cmax = clmax = 0;
  for(n=1, tot0=0.0, tot1=0.0; n<=N; n+=v) {
    compkcircb(&kc0b, &kc1b,&kc2b,v,n,mu,sigma,mdv,mtv);
    if(n==1)
      tot0 = 0.0;
    else
      tot0 += tot1 * lar;
      
    if(mu[n])
      tot1 += mu[n]*phi[n]/(((int_double) n) * sigma[n]);
    
    k2b = kc2b;
    k1b = kc1b - ccf* tot1;
    k0b = kc0b + ctm - ccf*tot0;

    /* for n<=y<=n+v,
       h_v(y,y) = k0b + k1b log(y/n) + k2b log(y/n)^2  */
    
    lar = log1p(((int_double) v)/n);
    
    if(n>ittil) {
      lint.left=0; lint.right = lar.right; 
      wabs = abs(k0b+lint*k1b+lint*lint*k2b);

      x = n;
      y = -((x+v)*wabs).right;
      if(y>cmax)
	cmax = y;
      y = -(y-c1*log(x)).right;
      if(y>clmax)
	clmax = y;

      if(flag)
	printf("%ld %g\n",n,-((n+v)*wabs).right);
    } else {
      if((n==7 || n==8 ||n ==1) && v==1)
	iter = ittil2/n+1;
      else
	iter = ittil/n+1;
      step = ((int_double) v)/iter;
      lint.left=0;
      
      for(j=0, stx=step; j<iter; j++, stx+=step) {
	lst = log1p(stx/n);
	lint.right=lst.right;
	wabs = abs(k0b+lint*k1b+lint*lint*k2b);

	x = n + stx - step;
	y = -((n+stx)*wabs).right;
	if(y>cmax)
	  cmax = y;
	y = -(y-c1*log(x)).right;
	if(y>clmax)
	  clmax = y;
	
	if(flag)
	  printf("%g %g\n",-x.right,y);

	lint.left=lst.left;
      }
    }
  }

  printf("Let boldh_v(y,y) = sum_{d: (d,v)=1} mu^2(d)/sigma(d)^2\n");
  printf("                               (tildem_{dv}(y/d) - zeta(2) sigma(d v)/d v)^2 \n");
  printf("For all y<=%ld,\n",N);
  printf("boldh_v(y,y) <= %.10g/y\n",cmax);
  printf("boldh_v(y,y) <= (%.10g + %.10g log y)/y\n",clmax,-c1.right);
}

    
