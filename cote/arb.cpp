#include "arb.h"
#include "primesieve.hpp"

uint64_t prime_count=0;
arb_t theta,tmp;
int64_t prec;

void callback(uint64_t p)
{
  prime_count++;
  arb_log_ui(tmp,p,prec);
  arb_add(theta,theta,tmp,prec);
}

int main(int argc, char **argv)
{
  if(argc!=3)
    return 0;
  uint64_t maxn=atol(argv[1]);
  prec=atol(argv[2]);
  primesieve::callback_primes(2,maxn,callback);
  printf("pi(%lu)=%lu\n",maxn,prime_count);
  printf("theta(%lu) in ",maxn);
  arb_printd(theta,40);
  printf("\n");
  return 0;
}
