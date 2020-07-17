#include "arb.h"
#include <primesieve.hpp>
#include "inttypes.h"


arb_t tmp1,tmp2,tmp3,sum;
int64_t prec;

uint64_t last_prime,last_2prime,count;
void callback(uint64_t p)
{
  if(p==last_prime+2)
    {
      count++;
      last_2prime=last_prime;
      arb_set_ui(tmp1,p);
      arb_inv(tmp2,tmp1,prec);
      arb_add(sum,sum,tmp2,prec);
      arb_set_ui(tmp1,p-2);
      arb_inv(tmp2,tmp1,prec);
      arb_add(sum,sum,tmp2,prec);
    }
  last_prime=p;
}

int main(int argc, char **argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <P> <prec>\n",argv[0]);
      return 0;
    }
  uint64_t P=atol(argv[1]);
  prec=atol(argv[2]);
  arb_init(sum);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);

  last_prime=3;
  count=0;
  primesieve::callback_primes(5,P+2,callback);

  printf("B(%lu)=",P);arb_printd(sum,20);
  printf("\n%lu twin primes found.\nlast twin prime found = %lu\n",count,last_2prime);

  return 0;
}
