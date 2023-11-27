#include "stdio.h"
#include "stdlib.h"
#include "stdbool.h"
#include "primesieve.h"

void set_bit(uint64_t ptr, uint64_t *sieve, uint64_t *bits)
{
  sieve[ptr>>6]|=bits[ptr&63];
}

int test_bit(uint64_t ptr, uint64_t *sieve, uint64_t *bits)
{
  return sieve[ptr>>6]&bits[ptr&63];
}

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

  for(uint64_t n=0;n<sieve_len;n++)
    sieve[n]=0; // 0 = square free, 1 = not

  uint64_t bits[64];
  bits[0]=1;
  for(uint64_t n=1;n<64;n++)
    bits[n]=bits[n-1]<<1;

  primesieve_iterator it;
  primesieve_init(&it);

  uint64_t prime,p2;

  while (true)
    {
      prime = primesieve_next_prime(&it);
      p2=prime*prime;
      if(p2>end)
	break;
      uint64_t ptr=p2;
      while(ptr<=end)
	{
	  set_bit(ptr,sieve,bits);
	  ptr+=p2;
	}
    }

  for(uint64_t n=start;n<=end;n++)
    if(test_bit(n,sieve,bits)==0)
      printf("%lu is square free.\n",n);

  for(uint64_t n=0;n<sieve_len;n++)
    printf("sieve[%lu]=%lX\n",n,sieve[n]);
  for(uint64_t n=0;n<64;n++)
    printf("bits[%lu]=%lX\n",n,bits[n]);
  
  primesieve_free_iterator(&it);

  return 0;

}
