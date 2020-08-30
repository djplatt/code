#include "stdlib.h"
#include "stdio.h"
#include "inttypes.h"
#include "pthread.h"

/*
#define MUTEX_MULT_SHIFT (12)
#define MUTEX_MULT (1LL<<MUTEX_MULT_SHIFT) 

pthread_mutex_t *mutexes;
*/

uint64_t N,n_threads;
uint8_t *ns;

void crossout(uint64_t n)
{
  /*
  printf("Crossing out %lu\n",n);
  asm goto (
	    "xor %%rdx,%%rdx\n\t"
	    "movq %0,%%rax\n\t"
	    "divq 18\n\t" /* rdx has remainder, rax the quotient *//*
	    "testq %%rdx,%%rdx\n\t"
	    "je %l1\n\t"
	    "movq ns(%%rip),%%rdx\n\t"
	    "xor %%bl,%%bl\n\t"
	    "movb %%bl,(%%rax,%%rdx)"
	    : /* no outputs *//*
	    : "r" (n)
	    : "rax","rdx","bl"
	    : no_conj);
 no_conj:
  return;
}
*/

  uint64_t res,quo;
  res=n%18;
  quo=n/18;
  if(res==0)
    {
      //uint64_t mutex_ptr=quo>>mutex_shift;
      //printf("Crossing out 18*%lu so grabbing mutex %lu\n",quo,mutex_ptr);
      //pthread_mutex_lock(mutexes+mutex_ptr);
      ns[quo]=0;
      //pthread_mutex_unlock(mutexes+mutex_ptr);
    }
}

void *thread_function(void *arg)
{
  uint64_t *thread_no_ptr=(uint64_t *)arg;
  uint64_t thread_no=thread_no_ptr[0];
  printf("In thread function with thread=%lu\n",thread_no);
  for(uint64_t a=thread_no+1;;a+=n_threads)
    {
      //printf("Thread %lu doing a=%lu\n",thread_no,a);
      uint64_t a2=a*a;
      if(a2>N)
	break;
      for(uint64_t b=a;;b++)
	{
	  //printf("Thread %lu doing b=%lu\n",thread_no,b);
	  uint64_t a2b2=a2+b*b;
	  if(a2b2>N) break;
	  crossout(a2b2);
	  for(uint64_t alpha=1;;alpha<<=1)
	    {
	      //printf("Thread %lu doing alpha=%lu\n",thread_no,alpha);
	      uint64_t a2b2a=a2b2+alpha;
	      if(a2b2a>N) break;
	      crossout(a2b2a);
	      for(uint64_t beta=1;;beta<<=1)
		{
		  //printf("Thread %lu doing beta=%lu\n",thread_no,beta);
		  uint64_t a2b2ab=a2b2a+beta;
		  if(a2b2ab>N) break;
		  crossout(a2b2ab);
		}
	    }
	}
    }
}

int main(int argc, char**argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <M> <T>\nWhere N=18*2^M and n_threads=2^T.\n",argv[0]);
      return(0);
    }
  printf("Looking for solutions  congruent to 0 mod 18\n");
  N=18LL<<atol(argv[1]);
  printf("Running with N=%lu\n",N);
  n_threads=1LL<<atol(argv[2]);
  printf("Using %lu threads.\n",n_threads);
  ns=(uint8_t *) malloc(sizeof(uint8_t)*(N/18+1));
  for(uint64_t n=0;n<=N/18;n++)
    ns[n]=1;

  pthread_t *thread_id;
  thread_id=(pthread_t *)malloc(sizeof(pthread_t)*n_threads);
  
  uint64_t *thread_indices;
  thread_indices=(uint64_t *)malloc(sizeof(uint64_t)*n_threads);
  
  // number the threads from 0 to n_threads-1
  for(uint64_t n=0;n<n_threads;n++)
    thread_indices[n]=n;
  
  // fire off the threads
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_create(thread_id+thread,NULL,thread_function,thread_indices+thread);
  
  // wait for all threads to return
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_join(thread_id[thread],NULL);

  // now find if any 1's remain
  uint64_t count=0;
  for(uint64_t i=1,n=18;n<=N;i++,n+=18,count++)
    if(ns[i]==1)
      printf("Found one with n=%lu\n",n);
  printf("We sieved %lu cases.\n",count);

  return(0);
}
