#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "math.h"
#include "primesieve.h"

#define SIEVE_LEN ((uint64_t) 1<<31) // 1Gbyte appears to be the sweet spot
#define LAST_P2 (1849) // the largest p^2 we will try before giving up

bool nsfree[SIEVE_LEN+LAST_P2]; // nsfree[0]=>st1-LAST_P2

uint64_t st,en,st1; // st,en are start and end of square free vector
                    // st1 

void cross_out(uint64_t p)
{
  //printf("p=%lu\n",p);
  uint64_t p2=p*p;
  uint64_t ptr=p2-(st%p2);
  if(ptr==p2) ptr=0;
  while(ptr<SIEVE_LEN+LAST_P2)
    {
      //printf("Crossing out %lu\n",ptr+st);
      nsfree[ptr]=true;
      ptr+=p2;
    }
}

int main(int argc, char **argv)
{

  if(argc!=4)
    {
      printf("Usage:- %s <sieve size> <start n> <num its>.\n",argv[0]);
      printf("<sieve size> = L1 cache size per core in KBytes.\n");
      exit(0);
    }

  //PrimeSieve ps;

  primesieve_set_sieve_size(atol(argv[1]));
  st1=atol(argv[2]);
  uint64_t num_its=atol(argv[3]);
  if(st1%4!=0)
    {
      printf("Expect start to be 0 mod 4. Exiting.\n");
      exit(0);
    }

  if(st1<LAST_P2)
    {
      printf("Need start >=%lu. Exiting.\n",LAST_P2);
      exit(0);
    }

  st=st1-LAST_P2;

  en=st1+SIEVE_LEN; // end of sieve range +1

  printf("Starting sieve at n=%lu\n",st1);

  for(uint64_t it=0;it<num_its;it++)
    {
      //printf("sqfree runs from [%lu,%lu]\n",st,en-1);
      //printf("checking in [%lu,%lu]\n",st1,en-1);

      for(uint64_t i=0;i<SIEVE_LEN+LAST_P2;i++)
	nsfree[i]=false;

      primesieve_callback_primes(2,(uint64_t) sqrt(en-1),cross_out);

#ifdef DEBUG
      uint64_t count=0;
      for(uint64_t ptr=0;ptr<SIEVE_LEN+LAST_P2;ptr++)
	if(nsfree[ptr]) count++;
      printf("There were %lu out of %lu not square free %5.2e%%\n",count,SIEVE_LEN+LAST_P2,(double)(100*count)/(double) (SIEVE_LEN+LAST_P2));
#endif

      for(uint64_t n=st1,ptr=LAST_P2;n<en;)
	{ // n=0 mod 4, so don't check n-4!
	  if(nsfree[ptr-9]&&nsfree[ptr-25]&&nsfree[ptr-49]&&
	     nsfree[ptr-121]&&nsfree[ptr-169]&&nsfree[ptr-289]&&nsfree[ptr-361]&&
	     nsfree[ptr-529]&&nsfree[ptr-841]&&nsfree[ptr-961]&&nsfree[ptr-1369]&&
	     nsfree[ptr-1681]&&nsfree[ptr-1849])
	    printf("Failed with n=%lu\n",n);
	  n+=2;ptr+=2;
	  if(nsfree[ptr-4]&&nsfree[ptr-9]&&nsfree[ptr-25]&&nsfree[ptr-49]&&
	     nsfree[ptr-121]&&nsfree[ptr-169]&&nsfree[ptr-289]&&nsfree[ptr-361]&&
	     nsfree[ptr-529]&&nsfree[ptr-841]&&nsfree[ptr-961]&&nsfree[ptr-1369]&&
	     nsfree[ptr-1681]&&nsfree[ptr-1849])
	    printf("Failed with n=%lu\n",n);
	  n++;ptr++;
	  if(nsfree[ptr-4]&&nsfree[ptr-9]&&nsfree[ptr-25]&&nsfree[ptr-49]&&
	     nsfree[ptr-121]&&nsfree[ptr-169]&&nsfree[ptr-289]&&nsfree[ptr-361]&&
	     nsfree[ptr-529]&&nsfree[ptr-841]&&nsfree[ptr-961]&&nsfree[ptr-1369]&&
	     nsfree[ptr-1681]&&nsfree[ptr-1849])
	    printf("Failed with n=%lu\n",n);
	  n++;ptr++;
	}
      st+=SIEVE_LEN;
      en+=SIEVE_LEN;
      st1+=SIEVE_LEN;
    }  

  printf("sieve ended at n=%lu\n",en-SIEVE_LEN-1);

  return(0);
}
