#include "arb.h"
#include <primesieve.hpp>
#include "inttypes.h"

arb_t sum,tmp1,tmp2;
int64_t prec;

void callback(uint64_t p)
{
  arb_set_ui(tmp1,p);
  arb_inv(tmp2,tmp1,prec);
  arb_neg(tmp2,tmp2);
  arb_log1p(tmp1,tmp2,prec);
  arb_add(sum,sum,tmp1,prec);
}

int main(int argc, char **argv)
{
  printf("Command line:- %s ",argv[0]);
  for(uint64_t a=1;a<argc,a++)
    printf("%s ",argv[a]);
  printf("\n");
  if(argc!=4)
    {
      printf("Usage:- %s <start> <end> <prec>.\n",argv[0]);
      return 0;
    }

  arb_init(sum);arb_zero(sum);
  arb_init(tmp1);arb_init(tmp2);
  prec=atol(argv[3]);
  primesieve::callback_primes(atol(argv[1]),atol(argv[2]),callback);

  printf("Sum returned ");arb_printd(sum,40);printf("\n");
  arb_clear(sum);
  arb_clear(tmp1);
  arb_clear(tmp2);
  return 0;
}
