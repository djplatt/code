#include "../includes/arb_sum.h"

int main(int argc, char ** argv)
{
  int64_t prec=200; // default
  if(argc>=2)
    prec=atol(argv[1]);

  arb_t res;
  arb_init(res);
  uint64_t count=arb_sum(res,prec);
  printf("%lu records from Stdin summed to ",count);
  arb_printd(res,30);
  printf("\n");
  arb_clear(res);
  return(0);
}
