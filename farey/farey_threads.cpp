// compute Farey discrepancy sigma (p/q-l/m)^2
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>
#include "A005728_table.h" // from https://oeis.org/A005728/b005728.txt
#include "pthread.h"
#include "slint.h"
#include "inttypes.h"

// globals
long N,n_threads,M;
double DM,*s1,*s2; // sum p/q
long *ls;

#define CPU_TIME (getrusage(RUSAGE_SELF,&ruse), ruse.ru_utime.tv_sec + \
  ruse.ru_stime.tv_sec + 1e-6 * \
  (ruse.ru_utime.tv_usec + ruse.ru_stime.tv_usec))

typedef unsigned long ulong;
struct rusage ruse;
double t0,t1;

ulong A005728(ulong n) {
  // A005728 Number of fractions in Farey series of order n.
  // a(5)=11 because the fractions are 0/1, 1/5, 1/4, 1/3, 2/5, 1/2, 3/5, 2/3, 3/4, 4/5, 1/1.
  // a(n) = n(n+3)/2 - Sum(k = 2 to n, a([n/k])) = 1 + Sum_{i=1..n} phi(i).
  // a[0]=1, a[1]=2, 3, 5, 7, 11, 13, 19, 23, 29, 33, 43, 47, 59, 65, 73, 81,
  ulong k,sum;
  // replace below by global table from A005728_table.h
  //unsigned int A005728_table[]={1,2,3,5,7,11,13,19,23,29,33,43,47,59,65,73,81,97,103,121,129,141,151,173,181,201,213,231,243,271,279,309,325,345,361,385,397,433,451,475,491,531,543,585,605,629,651,697,713,755,775,807,831,883,901,941,965};
  if (n<A005728_table_size) return A005728_table[n];
  sum=0;
  for (k=2; k<=n; k++) sum+=A005728(n/k);
  return (n*(n+3))/2-sum;
}


// compute the next farey fraction after a/b in sequence F_n
void compute_next(long &c, long &d, long p, long q, long n)
{
  if(p==0)
    {
      c=1;d=n;
      //printf("Next after %ld/%ld in F_%ld is %ld/%ld.\n",p,q,n,c,d);
    }
  else
    {
      long r=q-InvMod(p,q);
      long m=(n-r)/q;
      d=m*q+r;
      c=(d*p+1)/q;
      //printf("Next after %ld/%ld in F_%ld is %ld/%ld.\n",p,q,n,c,d);
    }
}

// compute the farey fractions from thread_no/(2*no_threads)
// to (thread_no+1)/(2*no_threads)
// we dont know where l should start so we guess and adjust afterwards
void *farey(void *arg)
{
  long *thread_no_ptr=(long *)arg;
  long thread_no=thread_no_ptr[0];
  long a,b,c,d,m,cc,dd,k,l=thread_no*M/(2*n_threads); // our guess
  b=n_threads*2;
  a=thread_no;
  long a_stop=thread_no+1,b_stop=n_threads*2;
  long g=GCD(a,b);
  a/=g;b/=g;
  g=GCD(a_stop,b_stop);
  a_stop/=g;b_stop/=g;
  compute_next(c,d,a,b,N);
  double dist,tmp;
  long double ss1=0.0,ss2=0.0;  // to avoid referencing s1 etc. in loop
  while((a!=a_stop)||(b!=b_stop))
    {
      dist=a;dist/=b;
      ss1+=dist; // sum p/q needed to correct afterwards
      tmp=dist-(double)l/DM;
      ss2+=tmp*tmp; // our first guess at the discrepancy
      l++;
      cc=c; dd=d;
      k=(N+b)/d;
      d=k*d-b;
      c=k*c-a;
      a=cc; b=dd;
    }
  s1[thread_no]=ss1;
  s2[thread_no]=ss2;
  ls[thread_no]=l-thread_no*M/(2*n_threads);
  
}

// this is the adjustment needed if our guessed l wasn't right
double dl_term(long th, long l0, long l1, double ss1, double ss2)
{

  long l=th*M/(2*n_threads); // we actually started with l
  long d=l-l0;
  double res=d/DM*((double)(l1-l0)*(double)(l1+l0-1)/DM+d*(double)(l1-l0)/DM-2*ss1);
  return(ss2-res);
}
  
  
int main(int argc, char* argv[]) 
{
  if(argc!=3)
    {
      printf("Usage:- %s <n> <no threads>.\n",argv[0]);
      return 0;
    }
  printf("Command line:- %s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  N=atol(argv[1]);
  n_threads=atol(argv[2]);
  M=A005728(N)-1;
  printf("m=%ld\n",M);
  DM=M;
  s1=(double *)malloc(sizeof(double)*n_threads);
  s2=(double *)malloc(sizeof(double)*n_threads);
  ls=(long *)malloc(sizeof(long)*n_threads);
  for(int i=0;i<n_threads;i++)
    {
      s1[i]=0.0;s2[i]=0.0;ls[i]=0;
    }
  // a vector to hold the thread id's
  pthread_t *thread_id;
  thread_id=(pthread_t *)malloc(sizeof(pthread_t)*n_threads);
  
  // a vector containing 0..n_threads-1 to pass to the sieve routine
  // so it knows which a's to sieve
  long *thread_indices;
  thread_indices=(long *)malloc(sizeof(long)*n_threads);
  
  // number the threads from 0 to n_threads-1
  for(int n=0;n<n_threads;n++)
    thread_indices[n]=n;
  
  // fire off the threads
  for(int thread=0,err;thread<n_threads;thread++)
    {
      err=pthread_create(thread_id+thread,NULL,farey,thread_indices+thread);
      if(err!=0)
	{
	  printf("Error creating thread %d. Exiting.\n",thread);
	  return 0;
	}
    }

  // wait for all threads to complete
  for(int thread=0;thread<n_threads;thread++)
    pthread_join(thread_id[thread],NULL);

  double tot=s2[0]; // the zero'th result needs no adjustment
  long l0=ls[0];
  for(int i=1;i<n_threads;i++)
    {
      tot+=dl_term(i,l0,l0+ls[i],s1[i],s2[i]);
      l0+=ls[i];
    }
  tot*=2.0; // we've only done to fraction 1/2-epsilon
  printf("Total = %8.6e %8.6e\n",tot,tot*N);

  return 0;
}
