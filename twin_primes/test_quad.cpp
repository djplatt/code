#include "inttypes.h"
#include "gmp.h"
#include "acb.h"
#include "quad.h"

// use DJP's implementation of Molin's double exponential quadrature
// to integrate x^2 from 0 to 1 rigorously.

// need real version of test function for integration
void mysqr(arb_t res, const arb_t x, int64_t prec)
{
  arb_pow_ui(res,x,2,prec);
}
// ned complex version to compute maxd
void mysqr(acb_t res, const acb_t x, int64_t prec)
{
  acb_pow_ui(res,x,2,prec);
}
int main(int argc, char **argv)
{
  if(argc!=4)
    {
      printf("Usage:- %s <n> <m> <prec>.\n",argv[0]);
      printf("n is how many steps to use in quadrature.\n");
      printf("m is how many segments to divide circle |z|=2 into to find maxd.\n");
      exit(0);
    }
  int64_t n=atol(argv[1]);
  int64_t m=atol(argv[2]);
  int64_t prec=atol(argv[3]);
  arb_t res,maxd,lo,hi;
  arb_init(res);arb_init(maxd);arb_init(lo);arb_init(hi);
  arb_set_ui(lo,0);
  arb_set_ui(hi,1);
  arb_maxd(maxd,mysqr,lo,hi,prec,m);
  printf("Maximum aound circle was ");arb_printd(maxd,prec);printf("\n"); 
  molin_int(res,n,mysqr,maxd,lo,hi,prec);
  printf("int x^2 dx, x=0..1 in ");arb_printd(res,prec);printf("\n");
  return(0);
}
