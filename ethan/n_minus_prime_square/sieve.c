// identify square-free numbers
// then printout n in [start,end] where n-p^2 is square-free
// for all prime p in [2,sqrt(n)]

#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "math.h"
#include "primesieve.h"

// set_bit sets the ptr'th bit of sieve to 1
// test_bit returns > 0 iff ptr'th bit of sieve is set
#define set_bit(ptr,sieve,bits) sieve[ptr>>6]|=bits[ptr&63]
#define test_bit(ptr,sieve,bits) ((uint64_t) sieve[ptr>>6]&bits[ptr&63])

int main (int argc, char** argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <start> <end>. Exiting.\n",argv[0]);
      exit(0);
    }

  uint64_t start=atol(argv[1]);
  uint64_t end=atol(argv[2]);

  if(end<start)
    {
      printf("Must have end (%lu) >= start (%lu). Exiting.\n",end,start);
      exit(0);
    }

  uint64_t sieve_len=(end>>6)+1;
  
  uint64_t *sieve=(uint64_t *)malloc(sizeof(uint64_t)*sieve_len);
  if(!sieve)
    {
      printf("Failed to allocate memory for sieve. Exiting.\n");
      exit(0);
    }

  size_t n_primes;
  
  int* primes = (int*) primesieve_generate_primes(2, sqrt((double) end)+1, &n_primes, INT_PRIMES);
  
  for(uint64_t n=0;n<sieve_len;n++)
    sieve[n]=0; // 0 = square free, 1 = not

  uint64_t bits[64];
  bits[0]=1;
  for(uint64_t n=1;n<64;n++)
    bits[n]=bits[n-1]<<1;

  for(uint64_t n=0;n<n_primes;n++)
    {
      uint64_t prime=primes[n];
      uint64_t p2=prime*prime;
      uint64_t ptr=p2;
      while(ptr<=end)
	{
	  set_bit(ptr,sieve,bits);
	  ptr+=p2;
	}
    }

  // only need to check n = {2,3} mod 4
  // but {0,1} get eliminated by p=2,3 resp. anyway
  
  for(uint64_t n=start;n<=end;n++)
    for(uint64_t nthp=0;;nthp++)
      {
	uint64_t p=primes[nthp];
	uint64_t p2=p*p;
	if(p2>n)
	  {
	    printf("%lu passed.\n",n);
	    break;
	  }
	if(test_bit(n-p2,sieve,bits))
	  break;
	else
	  continue;
      }
  


  return 0;

}
