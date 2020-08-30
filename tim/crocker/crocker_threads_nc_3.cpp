// way more efficient implementation due to P Zimmerman
#include "stdlib.h"
#include "stdio.h"
#include "inttypes.h"
#include "pthread.h"

uint64_t N; // largest number we will sieve
uint64_t n_threads; // how many threads
uint8_t *ns; // byte array, one byte per n<=N
uint64_t count=0;
uint64_t *counts;
//pthread_mutex_t print_mutex;

// we can express n as two squares and at most three powers of two
// Note that we don't bother ensuring exclusive use when writing to ns.
// Since all threads write the same value to ns (0) it doesn't matter
// what order the threads are serviced 
inline void crossout(uint64_t n)
{
    ns[n]=0;
}

void *set(void *arg)
{
  uint64_t *thread_no_ptr=(uint64_t *)arg;
  uint64_t thread_no=thread_no_ptr[0];
  printf("In thread set function with thread=%lu\n",thread_no);
  for(uint64_t a=thread_no+1;a<=N;a+=n_threads)
    ns[a]=1;
}

// this does all the work
// it runs through a\equiv thread_no+1 (n_threads)
// it crosses out a^2+b^2 (a2b2)
void *sieve(void *arg)
{
  uint64_t *thread_no_ptr=(uint64_t *)arg;
  uint64_t thread_no=thread_no_ptr[0];
  printf("In thread sieve function with thread=%lu\n",thread_no);
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
	}
    }
}


bool check_n(uint64_t n)
{
  if(ns[n]==0) return(false); // its sum of two squares
  uint64_t a=1,b;
  int64_t na=n-1,nab;
  while(na>0)  
    {
      if(ns[na]==0) return(false); // = 2 squares + a
      b=a+a;
      nab=na-b;
      while(nab>0)
	{
	  if(ns[nab]==0) return(false); // =2 squares + a + b
	  nab-=b;
	  b+=b;
	}
      na-=a;
      a+=a;
    }
  return(true);
}

void *check(void *arg)
{
  uint64_t *thread_no_ptr=(uint64_t *)arg;
  uint64_t thread_no=thread_no_ptr[0];
  printf("In thread check function with thread=%lu\n",thread_no);
  for(uint64_t n=thread_no+2;n<=N;n+=n_threads)
    if(check_n(n))
      {
	//pthread_mutex_lock(&print_mutex);
	counts[thread_no]++;
	//printf("Can't write %lu as sum of two squares and zero, one or two powers of 2.\n",n);
	//pthread_mutex_unlock(&print_mutex);
      }
}

int main(int argc, char**argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <M> <T>\nWhere N=2^M and n_threads=2^T.\n",argv[0]);
      return(0);
    }
  N=1LL<<atol(argv[1]);
  printf("Running with N=%lu\n",N);
  n_threads=1LL<<atol(argv[2]);
  printf("Using %lu threads.\n",n_threads);
  ns=(uint8_t *) malloc(sizeof(uint8_t)*(N+1));

  // set all entries in ns to 1.
  // should multi thread this too, but it takes neg. time
  for(uint64_t n=0;n<=N;n++)
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
  
  // fire off the threads to set ns[n]
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_create(thread_id+thread,NULL,set,thread_indices+thread);

  // wait for all threads to complete
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_join(thread_id[thread],NULL);

  // fire off the threads to sieve
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_create(thread_id+thread,NULL,sieve,thread_indices+thread);
  
  // wait for all threads to complete
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_join(thread_id[thread],NULL);

  /*
  if(pthread_mutex_init(&print_mutex,NULL)!=0)
    {
      printf("Failed to create mutex for printing. Exiting.\n");
      return(0);
    }
  */

  counts=(uint64_t *)malloc(sizeof(uint64_t)*n_threads);
  for(uint64_t n=0;n<n_threads;n++)
    counts[n]=0;
  // fire off the threads to check
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_create(thread_id+thread,NULL,check,thread_indices+thread);
  
  // wait for all threads to complete
  for(uint64_t thread=0;thread<n_threads;thread++)
    pthread_join(thread_id[thread],NULL);
  count=counts[0];
  for(uint64_t n=1;n<n_threads;n++)
    count+=counts[n];
  printf("We found %lu exceptions.\n",count);
  return(0);
}
