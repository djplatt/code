/* Compile with

g++ Hv.cpp -oHv -O2 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

This assumes that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that int_double14.2.h is in the
same directory as this file */

/*
Run with

./Hv

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
  int_double sum[3], sumh, L;
  long int N, Sqt, y, s, cent;
  int v,k;
  long i,j;
  int_double *k0, *k1,*k2, *md2, *mt2, sdphz,svz, ym;
  int_double k0init,k1init,k2init,S,NS,sint,sintm,c0,c1;
  int_double *mch, *m; 
  short int *mu, *isprime;
  long int *sigma;
  int_double **mdv, **mtv;
  int flag, checkfl;
  double maxerr;
  
  _fpu_rndd();

  if(argc<3) {
    N=1000000; v=1; flag = 0;  
  } else {
    N=atol(argv[1]); 
    v=atoi(argv[2]);
    if(argc<4)
      flag = 0;
    else
      flag = atoi(argv[3]);
  }

  Sqt=N;
  isprime = (short int *) calloc(Sqt+1,sizeof(short int));
  mu = (short int *) calloc(Sqt+1,sizeof(short int));
  sigma = (long int *) calloc(Sqt+1,sizeof(long int));
  fillisprime(isprime,Sqt+1); 
  fillmublock(mu,isprime,0,Sqt+1);
  fillsigmablock(sigma,isprime,0,Sqt+1);
  
  mdv=(int_double **) calloc(Sqt+1,sizeof(int_double *));
  mtv=(int_double **) calloc(Sqt+1,sizeof(int_double *));
  for(i=1; i<=Sqt; i++) {
    mdv[i] = (int_double *) calloc(Sqt/i+1,sizeof(int_double));
    mtv[i] = (int_double *) calloc(Sqt/i+1,sizeof(int_double));
    fillmdot(mdv[i],isprime,v*i,Sqt/i,mu,sigma);
    fillmtilde(mtv[i],mdv[i],v*i,Sqt/i,mu,sigma,sigma[i]*sigma[v]);
  }

  md2 = (int_double *) calloc(Sqt+1,sizeof(int_double));
  mt2 = (int_double *) calloc(Sqt+1,sizeof(int_double));
  k0 = (int_double *) calloc(Sqt+1,sizeof(int_double));
  k1 = (int_double *) calloc(Sqt+1,sizeof(int_double));
  k2 = (int_double *) calloc(Sqt+1,sizeof(int_double));
  
  for(i=1; i<=Sqt; i+=v)     
    mt2[i] = -(d_pi*d_pi/6.0)*((int_double) (sigma[i]*sigma[v]))/(i*v);

  if(v==1) {
    sdphz = d_pi*d_pi/6.0;
    svz = d_pi*d_pi/6.0;
  } else {
    sdphz = d_pi*d_pi/2.0;
    svz = d_pi*d_pi/4.0;
  }
    
  fillk(k0,k1,k2,&k0init,&k1init,&k2init,
	  v,0,N,Sqt,md2,mt2,
	isprime,mu,sigma,mdv,mtv);

  m = (int_double *) calloc(N+1,sizeof(int_double));
  mch = (int_double *) calloc(N+1,sizeof(int_double));
  
  fillm(m, v, N, mu);
  fillmch(mch, m, v, N, mu);

  if(v==1) {
    c0 = -0.0495; c1 = -0.998982;
  } else {
    c0 = 1.31742;
    c1 = - 1.817075;
  }

  checkfl=1; maxerr = 0.0; sumh = 0.0;
  for(y=1; y<N; y+=v) {
    sum[0] = sum[1] = sum[2] = 0; 
    for(s=1; s<=y; s+=v) {
      L = log1p(((int_double) (y%s))/(s*(y/s)));
      sum[0] += 
	(k0[y/s - 1] +
	 k1[y/s - 1]*L +
	 k2[y/s - 1]*L*L)/s;
      sum[1] += 
	(k1[y/s - 1] +
	 2*k2[y/s - 1]*L)/s;
      sum[2] += 
	k2[y/s - 1]/s;
    }
    sumh += ((int_double) 1.0)/y;

    if(y<100000)
      cent = 1000000/y;
    else
      cent = 5;

    S=0;
    for(j=0; j<cent; j++) {
      ym = v*((int_double) (j+1))/cent;
      NS = log1p(ym/y); S.right = NS.right;
      sint = sum[0] + sum[1]*S + sum[2]*S*S;
      sintm = sint - sumh*sdphz;
      sintm += svz*log(y+ym);
      if(flag && !j)
	printf("%ld %.10g %.10g\n",
	       y,-sint.right, -sintm.right);
      if(y>=2) {
	if(-(1-sintm/c1).right>
	   pow(y+ym,-(1.0 +((int_double) 1.0)/(v+1))).left) {
	  checkfl=0;
	  if(flag)
	    printf("The second bound may not hold around %g\n",(y+ym).left);
	}
	if(-((sint-c0)*(y+ym)).right > maxerr)
	  maxerr = -((sint-c0)*(y+ym)).right;
       }
      S.left = NS.left;
    }
  }

  printf("For 2<=Y<=%ld, v=%d,\n",N,v);
  printf("\t sum_{1<=s<=Y: (s,v)=1} h_v(Y/s)/s <= %.10g + %.10g/Y\n",
	 -c0.right,maxerr);
  if(checkfl) {
    printf("and sum_{1<=s<=Y: (s,v)=1} (h_v(Y/s)-(zeta(2) sigma(v)/phi(v)))/s is at most\n");
    if(v==1)
      printf("\t  - (zeta(2) sigma(v)/v) log y - %.7g (1 - y^(-3/2))\n",-(-c1).right);
    else
      printf("\t  - (zeta(2) sigma(v)/v) log y - %.7g (1 - y^(-4/3))\n",-(-c1).right);
  }
}
