#include "quad.h"
#include "acb.h"
#include <stdbool.h>
#include "inttypes.h"
#define PREC (100)

arb_t b_pi,h2;

void acb_f1(acb_t res, const acb_t z, int64_t prec)
{
  acb_sinc(res,z,prec);
  acb_mul_arb(res,res,b_pi,prec);
}

void acb_f2(acb_t res, const acb_t z, int64_t prec)
{
  static bool init=false;
  static acb_t tmp;
  if(!init)
    {
      init=true;
      acb_init(tmp);
    }
  acb_mul_arb(tmp,z,h2,prec);
  acb_cosh(res,tmp,prec);
}


int main()
{
  arb_t tmp,maxd,hi,lo,err;
  arb_init(tmp);arb_init(maxd);arb_init(hi);arb_init(lo);arb_init(err);
  arb_init(b_pi);
  arb_init(h2);
  arb_const_pi(tmp,PREC);
  arb_set_ui(h2,32); // b
  arb_div(b_pi,h2,tmp,PREC); // b/pi
  arb_set_ui(h2,2); //h/2
  arb_set_ui(lo,0);
  arb_set_ui(hi,100);

  arb_maxd2(maxd,acb_f1,acb_f2,lo,hi,PREC,100);
  comp_error(err,maxd,600,PREC);
  arb_printd(err,30);printf("\n");

  arb_clear(h2);
  arb_clear(b_pi);
  arb_clear(tmp);
  arb_clear(lo);
  arb_clear(hi);
  arb_clear(maxd);
  arb_clear(err);
  return 0;
}
