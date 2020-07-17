/* Checks whether

|\sum_{n\leq x} \mu(n)/n|\leq sqrt(2/x)

using interval arithmetic.
*/

/* Compile with

g++ musumv4.cpp -omusumv2 -O1 -frounding-math -fomit-frame-pointer -finline-functions -mfpmath=387 -I$CRDIR -L$CRDIR -lcrlibm

This assumes that crlibm is installed (in directory $CRDIR), that the
underlying processor is Intel, and that int_double12.0.h is in the
same directory as this file */

/* Test with

./musum 1000000000

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../includes/int_double12.0.h"
#include <omp.h>

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define min(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

typedef unsigned long ulong;

void fillisprime(short int *isprime, long int N)
/*sets isprime[i]=1 if i is a prime, isprime[i]=0 otherwise,
  for 2<=i<N */
{
  long int i,j,k,p;
  long int M, M2;
  long int *plist;
  int th_id;
  double t0, t1;

    
  isprime[2] = isprime[3] = isprime[5] = 1;
  
#pragma omp parallel
  {
  #pragma omp for schedule(dynamic, 100)
    for(i=6; i<N-5; i+=6) 
      isprime[i+1]=isprime[i+5]=1;
  }
  if(i<N-1)
    isprime[i+1]=1;

 M = sqrt(N); M2 = sqrt(M);
 
 for(i=2; i<M2; i++)
    if(isprime[i])
      for(j=i; j*i<M; j++)
        isprime[j*i]=0;

   k=0;
   for(i=2; i<M; i++)
     k+=isprime[i];

 plist = (long int *) calloc(k,sizeof(long int));
 for(i=2,j=0; i<M; i++)
   if(isprime[i])
     plist[j++] = i;

     t0 = omp_get_wtime();
     
#pragma omp parallel private(j,p)
 {
#pragma omp for schedule(dynamic, 1)
   for(i=2; i<k; i++) {
     p = plist[i];
     if((p%6)==1) {
       for(j=p*p; j<N-4*p; j+=6*p)
	 isprime[j]=isprime[j+4*p]=0;
       if(j<N) isprime[j]=0;
     }  else {
       for(j=p*p; j<N-2*p; j+=6*p) 
	 isprime[j]=isprime[j+2*p]=0;
       if(j<N) isprime[j]=0;
     }
     th_id = omp_get_thread_num();
     /*     printf("%ld %ld %d\n",i,p,th_id); */
   }
 }
     t1 = omp_get_wtime();

     printf("Central time: %g\n",t1-t0);
 free(plist);
}

void initmudiv(short int *initmu, long *initpr, short int *isprime,
	       long int minsq, long int minpr, long int M)
/*sets initmu[0],...,initmu[M-1] and initpr[0],...,initpro[M-1] where 
   where M = \prod_{p<minsq} p  *   \prod_{p<minpr} p */
/* assumes initmu and initpr have
  been allocated; assumes minpr>=minsq, minpr>=3*/
/* assumes isprime is filled correctly for m<minpr */
/*initmu[k] = 0 if p^2|k for some p<minsq
  initmu[k] = \prod_{p|k, p<minpr} (-1) otherwise*/
/* initpro[k] = \prod_{p|k, p<minpr} p if p^2 does not divide k for any p<minsq 
 */
{
  long int i, j, p;

  if(minsq>2)
    for(i=0; i<M; i+=4) {
      initmu[i] = 0; initmu[i+1] = initmu[i+3] = 1; initmu[i+2] = -1;
      initpr[i+1] = initpr[i+3] = 1; initpr[i+2] = 2;
    }
  else
    for(i=0; i<M; i+=2) {
      initmu[i] = -1; initmu[i+1] = 1;
      initpr[i] = 2; initpr[i+1] = 1;
    }

  

  for(p=3; p<minpr; p+=2)
    if(isprime[p])
      for(j=p; j<M; j+=p) {
	initmu[j] = - initmu[j];
	initpr[j]*=p;
      }

  for(p=3; p<minsq; p+=2)
    if(isprime[p])
      for(j=p*p; j<M; j+=p*p)
	initmu[j] = 0;
}

void copyblock(ulong *mus, ulong *sqf, short int *initmu,
	       long *tabla, long *initpr, long m, long n, long M)
{
  long int i;
  int a, b;
  ulong u;

  /* can optimize... */
  memset(mus,0,m/8+1);
  memset(sqf,255,m/8+1);

  for(i=0,a=0,b=0,u=1; i<m; i++) {
    switch(initmu[i%M]) {
    case 1 : mus[a] |= u;
      tabla[i] = (n+i)/initpr[i%M]; break;
    case -1 : mus[a] &= ~u;
      tabla[i] = (n+i)/initpr[i%M]; break;  
    case 0: sqf[a] &= ~u;
    }
    b++; u<<=1; if(b==(8*sizeof(ulong))) {
      b=0; a++; u=1;
    }
  }
}

void fillmublock(ulong *mus, ulong *sqf,
		 short int *isprime,
		 long int n, long int m,
		 short int *initmu, long *initpr,
		 long int minsq, long int minpr, long int M)
/* sets mus[0].. mus[m/sizeof(long)],
        sqf[0].. sqf[m/sizeof(long)] in such a way that
	the bth bit of mus[j] stores (1+mu(n+sizeof(long)*j+b))/2
            if (n+sizeof(long)*j+b) is square-free, and
        the bth bit of sqf[j] stores zhether
               n+sizeof(long)*j+b is square-free */

/* assumes isprime is filled and valid up to and including sqrt(n+m-1) */
/* assumes initmu is filled and valid from 0 up to M-1,
   where M = \prod_{p<minsq} p  *   \prod_{p<minpr} p */
/* assumes n and m are divisible by M */
/* assumes minsq <= minpr */
/* convention: mu(0) = 0 */
/* assumes minpr>2*/
{
  long int i,j;
  long int *tabla;
  long int maxp;
  long int p,p2;

  tabla = (long int *)calloc(m,sizeof(long int));
  /* printf("sizes: m %ld tabla %ld mun %ld\n",m,m*sizeof(long int),
     m*sizeof(short int));*/

  copyblock(mus, sqf, initmu, tabla, initpr, m, n, M);

  maxp = (long int) sqrt(n+m);

    for(p=minsq; p<minpr; p++)
      if(isprime[p]) {
	for(j=((n+p*p-1)/(p*p))*p*p-n; j<m; j+=p*p)
	  sqf[j/(8*sizeof(ulong))]  &= ~(((ulong) 1)<<(j%(8*sizeof(ulong))));;
      }
    

  
  for(p=minpr; p<=maxp; p+=2) 
    if(isprime[p]) {
      for(j=((n+p-1)/p)*p-n; j<m; j+=p) {
	mus[j/(8*sizeof(ulong))]  ^= (((ulong) 1)<<(j%(8*sizeof(ulong))));
	    tabla[j] /= p; 
	  }
      p2 = p*p;
      for(j=((n+p2-1)/p2)*p2-n; j<m; j+=p2)
	sqf[j/(8*sizeof(ulong))]  &= ~(((ulong) 1)<<(j%(8*sizeof(ulong))));;
    }


  for(i=0; i<m; i++)
     if(tabla[i]!=1)
 	mus[i/(8*sizeof(ulong))]  ^= (((ulong) 1)<<(i%(8*sizeof(ulong))));

  free(tabla);
}


	  int_double summu(long int N, double *maxr, int nthreads, long minsq, long minpr,
int chunk)
/* returns \sum_{n\leq N} \mu(n)/n */
/* stores the maximum of (\sum_{n\leq x} \mu(n)/n)*sqrt(x) for x\leq N+1 */
{
  long maxp, M,i, j,k,k2;
  short int *isprime, *initmu;
  ulong **mus, **sqf;
  long *initpr;
  long  m, p;
  int th_id;
  double t0, t1;
  int_double psum, sum, *psumar;
  int i0;
  double maxrat, pmax, *pmaxar;
  int_double ms;
  
  sum = 0.0;
  
  maxp = (long int) sqrt(N);

  isprime = (short int *) calloc(maxp+1,sizeof(short int));
  fillisprime(isprime,maxp+1); 
  
  for(M=1, p=2; p<minsq; p++)
    if(isprime[p])
      M *= p;
  
  for(p=2; p<minpr; p++)
    if(isprime[p])
      M *= p;

  printf("%ld\n",M);
  initmu = (short int *) calloc(M,sizeof(short int));
  initpr = (long int *) calloc(M,sizeof(long int));
  initmudiv(initmu, initpr, isprime, minsq, minpr, M);

  printf("%ld\n",M);
  m = ((chunk*((long int) sqrt(N))+M-1)/M)*M;

  mus = (ulong **) calloc(nthreads,sizeof(ulong *));
  sqf = (ulong **) calloc(nthreads,sizeof(ulong *));
  psumar = (int_double *) calloc(nthreads,sizeof(int_double));
  pmaxar = (double *) calloc(nthreads,sizeof(double));

  t0 = omp_get_wtime();

  maxrat=0.0;
  for(i=0; i<=N; i+=m*nthreads) {
    if(!((i/(m*nthreads))%100))
      printf("pecar %ld\n",i);
#pragma omp parallel private(k, k2, psum, pmax, i0, th_id)
    {
#pragma omp for schedule(dynamic, 1)
      for(j=i; j<i+m*nthreads; j+=m) {
	pmax = 0.0;
	i0 = ((j-i)/m);
	if(j<=N) {
	  mus[i0] = (ulong *) calloc(m/(8*sizeof(ulong))+1,sizeof(ulong));
	  sqf[i0] = (ulong *) calloc(m/(8*sizeof(ulong))+1,sizeof(ulong));
	  fillmublock(mus[i0],sqf[i0],isprime,j,m,initmu, initpr, minsq, minpr,M);
	  th_id = omp_get_thread_num();
	  k2 = min(j+m,N+1);
	  psum = 0.0;
	  /*	for(k=max(j,1); k<max(j,1)+100; k++)
		printf("%ld %d\t",k,mun[i0][k-j]);*/
	  for(k=max(j,1); k<k2; k++) {
	    if(sqf[i0][(k-j)/(8*sizeof(ulong))]&(((ulong) 1)<<((k-j)%(8*sizeof(ulong))))) {
	      if(mus[i0][(k-j)/(8*sizeof(ulong))]&(((ulong) 1)<<((k-j)%(8*sizeof(ulong)))))
		psum += 1.0/((int_double) k);
	      else 
		psum -= 1.0/((int_double) k);
	      if(max(-psum.left,-psum.right)>pmax)
		pmax = max(-psum.left,-psum.right);
	    }
	  }
	  psumar[i0] = psum; pmaxar[i0] = pmax;
	  free(mus[i0]); free(sqf[i0]);
	} else {
	  psumar[i0] = 0.0; pmaxar[i0] = 0.0;
	}
      }
	/*	printf("%ld %d: [%g,%g] total: [%g,%g]\n",j,th_id,
		psum.left,-psum.right,sum.left,-sum.right);*/
    }
 
    for(i0=0; i0<nthreads; i0++) {
      j = i + m*i0;
      if(-((((int_double) max(-sum.right,-sum.left)) + pmaxar[i0]) *
	   sqrt((int_double) j+m)).right> maxrat && j<=N) {
	  printf("slow %ld %ld %g %g %g %g\n",j,j+m,
		 max(-sum.right,-sum.left),
		 max(-(sum+pmaxar[i0]).right,-(sum+pmaxar[i0]).left),
		 -(max(-(sum+pmaxar[i0]).right,-(sum+pmaxar[i0]).left)*
		   sqrt((int_double) j+m)).right,maxrat);
	  ms = sum + psumar[i0];
	  mus[0] = (ulong *) calloc(m/(8*sizeof(ulong))+1,sizeof(ulong));
	  sqf[0] = (ulong *) calloc(m/(8*sizeof(ulong))+1,sizeof(ulong));
	  fillmublock(mus[0],sqf[0],isprime,j,m,initmu, initpr, minsq, minpr,M);
	  k2 = min(j+m,N+1);
	  for(k=max(j,1); k<k2; k++) { 
	    if(sqf[0][(k-j)/(8*sizeof(ulong))]&
	       (((ulong) 1)<<((k-j)%(8*sizeof(ulong))))) {
	      if(mus[0][(k-j)/(8*sizeof(ulong))]
		 &(((ulong) 1)<<((k-j)%(8*sizeof(ulong)))))
		sum += 1.0/((int_double) k);
	      else 
		sum -= 1.0/((int_double) k);  
	    }
	    if(-(max(-sum.left,-sum.right)*sqrt((int_double) k+1)).right
	       >maxrat && k>=3) {
	      maxrat = 
		-(max(-sum.left,-sum.right)*sqrt((int_double) k+1)).right;
	      printf("newmax %ld %g\n",k,maxrat); 
	    }
	  }
	  free(mus[0]); free(sqf[0]);
      }
      else sum += psumar[i0];
    }
  }
  
  t1 = omp_get_wtime();
  printf("Central time: %g\n",t1-t0);

  printf("Final result: %ld %.15g %.15g\n",i,sum.left,sum.right);
  free(mus);
  free(sqf);  
  free(psumar); free(pmaxar); 
  return maxrat;
}

/*
  
  long int i, j, m;
  int_double sq,lev;
  int *Mu;
  int flag=0;

  Mu = (int *)calloc(M,sizeof(int));
  for(i=N0; i<=N; i+=M) {
    if(i + M <= N)
      m = M;
    else m = N-i;
    if(!((i/M)%10))
      printf("Checking at %ld: [%20.18e,%20.18e] \n",i,sum.left,-sum.right);
    fillmublock(Mu,i,m);
    for(j=i; j<i+m; j++) 
      if(Mu[j-i]) {
	sum += int_double(Mu[j-i])/j;
	if(!flag) {
	  sq = 1.0/(2*sqrt(int_double(j+1)));
          //print_int_double_str("sq=",sq);exit(0);
	  if((-sum.right)>sq.left || (-sum.left)>sq.left) {
	    if(j>=3) {
	      flag=1;
	      printf("The first minor crime - at %ld: [%20.18e,%20.18e] [%20.18e,%20.18e]\n",
		     j,sum.left,-sum.right,sq.left,-sq.right);
	    }
	  }
	} else {
	  lev = sqrt(int_double(2.0)/(j+1));
	  if(-sum.right>lev.left || (-sum.left)>lev.left)
	    if(j>1)
	      printf("FELONY! at %ld: [%20.18e,%20.18e] [%20.18e,%20.18e]\n",j,sum.left,-sum.right,lev.left,-lev.right);
	} 
      }
  }
  free(Mu);
  
  return sum;
}
*/

main(int argc, char *argv[])
{
  int nthreads, th_id;
  /*  _fpu_rndd();*/
  long int N;
  long minsq, minpr;
  double maxr;
  int_double s;
  int chunk;

   _fpu_rndd();
  
  if(argc>1)
    N = atol(argv[1]);
  else 
    N = 1000000000000;
  
  if(argc>2) {
    minsq = atol(argv[2]);
    minpr = atol(argv[3]);
    if(minpr<minsq)
      minpr = minsq;
  } else {
    minsq = 5;
    minpr = 19;
  }


  if(argc>4)
    chunk = atoi(argv[4]);
  else
    chunk = 5;

  if(argc>5)
    omp_set_num_threads(atoi(argv[5]));

#pragma omp parallel private(th_id)
  {
    th_id = omp_get_thread_num();
    if ( th_id == 0 ) {
      nthreads = omp_get_num_threads();
      printf("There are %d threads\n",nthreads);
    }
  }

    
  printf("%ld %ld\n",minsq,minpr);
  s=summu(N, &maxr, nthreads, minsq, minpr, chunk);
  printf("[%.15g,%.15g]\n",s.left,-s.right);
}
