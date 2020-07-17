#include "stdlib.h"
#include "stdio.h"
#include "inttypes.h"
#include "pthread.h"

uint64_t N; // largest number we will sieve
uint64_t n_threads; // how many threads
uint8_t *ns; // byte array, one byte per n<=N with n\equiv 0 (18)

// we can express n as two squares and at most two powers of two
// if n \equiv 0 (18) then cross out (set to zero) the corresponding entry
// in the byte array ns
// Note that we don't bother ensuring exclusive use when writing to ns.
// Since all threads write the same value to ns (0) it doesn't matter
// what order the threads are serviced 
void crossout(uint64_t n)
{
  uint64_t res,quo;
  res=n%18;
  quo=n/18;
  if(res==0)
    ns[quo]=0;
}

// this does all the work
// it runs through a\equiv thread_no+1 (n_threads)
// it crosses out a^2+b^2 (a2b2)
//                a^2+b^2+2^alpha (a2b2a)
// and            a^2+b^2+2^alpha+2^beta (a2b2ab)
void *sieve(void *arg)
{
  uint64_t *thread_no_ptr=(uint64_t *)arg;
  uint64_t thread_no=thread_no_ptr[0];
  printf("In thread function with thread=%lu\n",thread_no);
  for(uint64_t a=thread_no+1;;a+=n_threads)
    {
      uint64_t a2=a*a;
      if(a2>N)
	break;
      for(uint64_t b=a;;b++)
	{
	  uint64_t a2b2=a2+b*b;
	  if(a2b2>N) break;
	  crossout(a2b2);
	  for(uint64_t alpha=1;;alpha<<=1)
	    {
	      uint64_t a2b2a=a2b2+alpha;
	      if(a2b2a>N) break;
	      crossout(a2b2a);
	      for(uint64_t beta=1;;beta<<=1)
		{
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

  // set all entries in ns to 1.
  // should multi thread this too, but it takes neg. time
  for(uint64_t n=0;n<=N/18;n++)
    ns[n]=1;

  // a vector to hold the thread id's
  pthread_t *thread_id;
  thread_id=(pthread_t *)malloc(sizeof(pthread_t)*n_threads);
  
  // a vector containing 0..n_threads-1 to pass to the sieve routine
  // so it knows which a's to sieve
  uint64_t *thread_indices;
  thread_indices=(uint64_t *)malloc(sizeof(uint64_t)*n_threads);
  
  // number the threads from 0 to n_threads-1
  for(uint64_t n=0;n<n_threads;n++)
    thread_indices[n]=n;
  
  // fire off the threads
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_create(thread_id+thread,NULL,sieve,thread_indices+thread);
  
  // wait for all threads to complete
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
