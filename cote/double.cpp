#include "primesieve.hpp"
#include <stdio.h>
#include <stdlib.h>

uint64_t prime_count=0;
double theta=0.0;

void callback(uint64_t p)
{
  prime_count++;
  theta+=log(p);
}

int main(int argc, char **argv)
{
  if(argc!=2)
    return 0;
  uint64_t maxn=atol(argv[1]);
  primesieve::callback_primes(2,maxn,callback);
  printf("pi(%lu)=%lu\n",maxn,prime_count);
  printf("theta(%lu) approx %20.18e\n",maxn,theta);
  return 0;
}
