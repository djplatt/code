/* Compile with

g++ boundhvxy.cpp -oboundhvxy -O2 -frounding-math -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

This assumes that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that int_double14.0.h is in the
same directory as this file */

/*
Run with

./betterhmin4

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "int_double14.2.h"
#include "mupsigma2.h"
#include "mdotmtildek4.h"

int_double *log1par;

void kap11vstart(long N, int_double *kap,
		 short *mu, long *sigma, int v)
/* initializes kap[n2] to
   kap_{1,1,v}(1,n2) = \sum_{r_2\leq n_2: (r,v)=1} mu(r_2)/sigma(r_2) */
{
  long n2;

  kap[1] = 1;

  for(n2=1+v; n2<=N; n2+=v) {
    if(mu[n2])
      kap[n2] = kap[n2-v] + ((int_double) mu[n2])/sigma[n2];
    else
      kap[n2] = kap[n2-v];
  }
}

void kap01vstart(long N, int_double *kap, int_double *m, int v, int_double mc)
/* initializes kap[n2] to

   -mc m(n2) = - mc \sum_{r\leq n_2: (r,v)=1} mu(r)/d(r) */
{
  long n2;
  
  for(n2=1; n2<=N; n2+=v)
    kap[n2] = -mc*m[n2];
}

void kapj0vbuild(long N, int_double *kapj0, int_double *kapj1, int v,
		 int_double firstkap)
/* initialize kapj0[n2] = kap_{j,0}(n_1,n_2) for given n_1 using
   kapj0[n2+v] = kapj1[n2] * log((n_2+v)/n_2) */
{
  long n2;

  kapj0[1] = firstkap;

  for(n2=1; n2<=N-v; n2+=v)
    kapj0[n2+v] = kapj0[n2] + kapj1[n2]*log1par[n2];
}

void kap0jvbuild(long n1, long N, int_double *kap0j, int_double *kap1j,
		 int v)
/* builds kap0j[n2] = kap_{0,j}(n_1+v, n_2) for given n1,
   given an old value kap0j[n2] = kap_{0,j}(n_1, n_2) 
   and kap1j[n2] = kap_{1,j,v}(n_1,n_2) */
{
  long n2;

  for(n2=1; n2<=N; n2+=v)
    kap0j[n2] += kap1j[n2]*log1par[n1];
}

void kap11vbuild(long n1, long N, int_double *kap,
	     short *mu, long *sigma, short *cop, int v)
/* for 1<=n2<=N, (n2,v)=1,
lets kap[n2] +=
       mu[n1+v]*mu[n2]/(sigma[n1+v]*sigma[n2]) if cop[n2]
                                                   (meaning (n1+v,n2)=1)
     kapnew[n2] = kapold[n2] otherwise */
/* v = 1 or v=2*/
{
  long n2;
  int_double tot;

  tot=0;
  for(n2=1; n2<=N; n2+=v) {
    if(cop[n2] && mu[n1+v] && mu[n2])
      tot+=	((int_double) (mu[n1+v]*mu[n2]))/(sigma[n1+v]*sigma[n2]);
    kap[n2] +=tot;
    
  }
}


#define conform(n,v) ((v)==1 ? (n) : ((n)%2 ? (n) : (n)-1))

main(int argc, char *argv[])
{
  long N, *sigma, Sqt, i, n2start, minthere;
  int v, R, S, B, S1, r, s;
  short *isprime, *mu, *cop;
  int_double *m;
  long n1,n2,rat;
  int_double hres, mc, hrat;
  int_double *kap00, *kap01, *kap10, *kap11;
  int_double sqrtip, *sqrtar, **sqrtor, **log1por;
  double minh,lowh, minrat, minratthere;
  int_double integ, neginteg, mininteg, step, igrand;
  long firstn;
  int loudflag;
  
  _fpu_rndd();

  if(argc<4) {
    N=100000; rat=1; v=1; loudflag = 0;
  } else {
    N=atol(argv[1]); rat=atol(argv[2]); v=atol(argv[3]);
    if(argc<5)
      loudflag = 0;
    else
      loudflag = atoi(argv[4]);
  }
    R=10000; B=5000; N=max(N,R); Sqt = sqrt(N+v);


  isprime = (short int *) calloc(Sqt+2,sizeof(short int));
  mu = (short int *) calloc(N+1,sizeof(short int));
  cop = (short int *) calloc(N+1,sizeof(short int));
  sigma = (long int *) calloc(N+1,sizeof(long int));

  fillisprime(isprime, Sqt+2);
  fillmublock(mu,isprime,0,N+1);
  fillsigmablock(sigma,isprime,0,N+1);

  m = (int_double *) calloc(N+1,sizeof(int_double));
  sqrtar = (int_double *) calloc(N+v+1,sizeof(int_double));
  sqrtor = (int_double **) calloc(N+1,sizeof(int_double));
  log1par = (int_double *) calloc(N+1,sizeof(int_double));
  log1por = (int_double **) calloc(N+1,sizeof(int_double));
  kap00 = (int_double *) calloc(N+1,sizeof(int_double));
  kap01 = (int_double *) calloc(N+1,sizeof(int_double));
  kap10 = (int_double *) calloc(N+1,sizeof(int_double));
  kap11 = (int_double *) calloc(N+1,sizeof(int_double));
  fillm(m, v, N, mu);
  mc = ( v==1 ? d_pi*d_pi/6 : d_pi*d_pi/4);

  /* time to initialize some variables for later usage */

  for(i=1; i<=N; i+=v) 
    log1par[i] = log1p(v/((int_double) i));
    
  for(i=1; i<=N; i+=v) {
    S = (i<=R ? R/i : 1);
    log1por[i] = (int_double *) calloc(S,sizeof(int_double));
    sqrtor[i] = (int_double *) calloc(S,sizeof(int_double));
    log1por[i][0] = log1p(v/(S*((int_double) i)));
    (log1por[i][0]).left = 0;
    sqrtor[i][0] = sqrt(i + ((int_double) v)/S);
    (sqrtor[i][0]).left = (sqrt((int_double) i)).left;
    for(r=1; r<S; r++) {
      log1por[i][r] = log1p(((r+1)*v)/(S*((int_double) i)));
      (log1por[i][r]).left = -(log1por[i][r-1]).right;
      sqrtor[i][r] = sqrt(i + (r+1)*((int_double) v)/S);
      (sqrtor[i][r]).left = -(sqrtor[i][r-1]).right;
    } /* we are computing log(1+x) and sqrt(x) on subintervals of [i,i+v] */
  }

  lowh=0;
  /* compute kapij[1][n] */

  neginteg = mininteg = integ = 0; minrat = 0; firstn=1;

  for(n1=1; n1<=N; n1+=v) {
    if(n1==1) {
      kap11vstart(N, kap11, mu, sigma, v);
      kap01vstart(N, kap01, m, v, mc);    
      kapj0vbuild(N, kap10, kap11, v, -mc);
      kapj0vbuild(N, kap00, kap01, v, (v==1 ? d_pi*d_pi/6 : d_pi*d_pi/2));
    } else {
      fillcopblock(cop,isprime,0,N+1,n1);
      kap0jvbuild(n1-v, N, kap01, kap11, v);
      kap0jvbuild(n1-v, N, kap00, kap10, v);
      kap11vbuild(n1-v, N, kap11, mu, sigma, cop, v);
      kapj0vbuild(N, kap10, kap11, v, -mc*m[n1]);
    }
    if(firstn<=N) {
      S1 = (n1<=B ? R/n1 : 1);   step = ((int_double) v)/S1;
      for(s=0; s<S1; s++) {
	minh = 1000000; minthere=firstn;

	sqrtip = sqrtor[n1][s];

	/*  We must start at rat*(1+s/S1)), or rather at the largest
	    integer \equiv 1 mod v that is no greater than rat*(1+s/S1)  */
	firstn = 1+v*((rat*(n1*S1+s)-S1)/(v*S1));

	for(n2=firstn; n2<=N; n2+=v) {
	  S = (n2<=R ? R/n2 : 1);
	  if(n2!=firstn)
	    r = 0;
	  else
	    r = (S*((rat*(n1*S1+s)-S1)%(v*S1)))/(v*S1);
	  for(r; r<S; r++) {
	    hres = kap00[n2]+kap10[n2]*log1por[n1][s]+kap01[n2]*log1por[n2][r]
	      + kap11[n2]*log1por[n1][s]*log1por[n2][r];
	    	/* hres is an interval in which h_v(x,y) lies
	   for x in [1+s*v/S1,1+(s+1)*v/S1]
	   and y in [n2+r*v/S1,n +(r+1)* v/S1] */
	    hres *= sqrtip*sqrtor[n2][r];
	    	/* we multiply hres by maximum of sqrt(x y) for x,y
	   in the intervals above */
	    if(hres.left<minh) {
	      minh = hres.left;
	      minthere = n2;
	    }
	  }
	}

	if(firstn<=N) {
	  if(loudflag)
	    printf("%.10g %.10g\n",(n1+(s*step)).left,minh);
	  if(minh<lowh)
	    lowh=minh;
	  igrand = step/((n1+(s*step))*(n1+((s+1)*step)));
      /* (x_1-x_0)/x_0 x_1 = 1/x_0 - 1/x_1 is the
	 integral of 1/x^2 from x_0 to x_1 */ 	  
	  integ += minh*igrand;
	  if(minh<0)
	    neginteg += minh*igrand;
	  mininteg += lowh*igrand;

	  hrat = minh/(n1+s*step);
	  if(hrat.left<minrat) {
	    minrat = hrat.left;
	    minratthere = (n1+s*step).left;
	  }
	}
      }
    }
  }

  printf("The following bounds concern f_{%d,%d}(x), defined as\n",r,v);
  printf("\t inf_{r x <= y <= %ld} sqrt(x y) h_%d(x,y) for r=%d.\n",N,v,r);
  printf("Then, for 1<=x<=%ld/r, where r=%d,\n",N,r);
  printf("f_{%d,%d}(x) >= %.10g\n",r,v,lowh);
  printf("and f_{%d,%d}(x)/x >= %.10g; that value may be reached near x=%.10g\n",r,v,minrat,minratthere);
  printf("The integral of f_{%d,%d}/x^2 from 0 to infty is at least %.10g\n",r,v,integ.left);
  printf("The integral of min(0,f_{%d,%d})/x^2 from 0 to infty is at least %.10g\n",r,v,integ.left);
  printf("Define F_{%d,%d}(x) = max_{1<=t<=x} (-f_{%d,%d}(t)).\n",r,v,r,v);
  printf("Then the integral of F_{%d,%d}(x)/x^2 from 0 to infty is at most  %.10g\n",r,v,-mininteg.left);
}
