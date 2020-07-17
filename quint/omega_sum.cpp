/// @example count_primes.cpp
/// This example shows how to count primes.

#include "stdlib.h"
#include "stdio.h"
#include <primesieve.hpp>
#include <stdint.h>
#include <iostream>

#define MAX_PRIMES (11)
#define MAX_PRIMORIAL (200560490130LL)
#define BOUND (6980000000LL)
uint8_t *counts;
uint64_t start,end;

void callback(uint64_t prime)
{
  //printf("Doing prime %lu\n",prime);
  for(uint64_t ptr=prime;ptr<=BOUND;ptr+=prime)
    counts[ptr]++;
  
}

int main(int argc, char ** argv)
{
  if(argc!=1)
    {
      printf("Usage:- %s\n",argv[0]);
      exit(0);
    }


  counts=(uint8_t *)malloc(sizeof(uint8_t)*(BOUND+1));
  if(!counts)
    {
      printf("Failed to allocate memory for counts. Exiting.\n");
      exit(0);
    }
  for(uint64_t i=0;i<BOUND;i++)
    counts[i++]=0;

  primesieve::callback_primes(2,BOUND,callback);

  //for(uint64_t n=1;n<=100;n++)printf("Omega(%lu)=%u\n",n,counts[n]);

  uint64_t fours[MAX_PRIMES+1];
  fours[0]=1;
  for(uint64_t ptr=1;ptr<=MAX_PRIMES;ptr++)
    fours[ptr]=fours[ptr-1]<<2;

  uint64_t sum=0;
  for(uint64_t ptr=2;ptr<=BOUND;ptr++)
    sum+=fours[counts[ptr]];

  printf("Sum 4^omega{n} n=[2,%lu]=%lu\n",BOUND,sum);
  return 0;
}
