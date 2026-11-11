#include "stdbool.h"
#include "inttypes.h"
#include "stdlib.h"
#include "stdio.h"
#include "primesieve.h"

int main(int argc, char **argv)
{

  uint64_t lim=1000;
  if(argc>1)
    lim=atol(argv[1]);
  bool *sieve;
  sieve=(bool *)malloc(sizeof(bool)*lim);
  for(uint64_t i=0;i<lim;i++)
    sieve[i]=true;
  primesieve_iterator it;
  primesieve_init(&it);
  uint64_t p=primesieve_next_prime(&it);
  while(true)
    {
      p=primesieve_next_prime(&it); // start at 3
      uint64_t p2=p*p;
      if(p2>lim)
	break;
      uint64_t ptr=p2-1;
      while(ptr<lim)
	{
	  sieve[ptr]=false;
	  ptr+=p2;
	}
    }
  primesieve_free_iterator(&it);

  uint64_t count=0;
  for(uint64_t ptr=2;ptr<lim;ptr+=16) // 3 mod 16
    if(sieve[ptr])
      count++;
  for(uint64_t ptr=3;ptr<lim;ptr+=16) // 4 mod 16
    if(sieve[ptr])
      count++;
  for(uint64_t ptr=6;ptr<lim;ptr+=16) // 7 mod 16
    if(sieve[ptr])
      count++;
  for(uint64_t ptr=7;ptr<lim;ptr+=16) // 8 mod 16
    if(sieve[ptr])
      count++;
  for(uint64_t ptr=10;ptr<lim;ptr+=16) // 11 mod 16
    if(sieve[ptr])
      count++;
  for(uint64_t ptr=14;ptr<lim;ptr+=16) // 15 mod 16
    if(sieve[ptr])
      count++;

  printf("There are %lu fundamental discriminants in [-%lu,-1].\n",count,lim);
  

  return 0;
}
