#include "acb.h"
#include "stdio.h"

acb_t tmp,tmp1;

void gamr(acb_t res, acb_t z, int64_t prec)
{
  acb_mul_2exp_si(tmp,z,-1);
  acb_gamma(res,tmp,prec);
}

void lam(acb_t res, acb_t z, int64_t prec)
{
  acb_zeta(tmp1,z,10);
  gamr(res,z,prec);
  acb_mul(res,res,tmp1,prec);
}

int main()
{
  int64_t prec=200;
  acb_init(tmp);
  acb_init(tmp1);
  arb_t rho;
  arb_init(rho);
  arb_set_ui(rho,2599761596147898551LL);
  arb_mul_2exp_si(rho,rho,-64);
  arb_add_ui(rho,rho,1301780508LL,prec);

  arb_t del;
  arb_init(del);
  arb_set_ui(del,1);
  arb_mul_2exp_si(del,del,-64);

  acb_t z;
  acb_init(z);
  acb_t tmp2;
  acb_init(tmp2);
  arb_set_d(acb_realref(z),0.5);
  arb_sub(acb_imagref(z),rho,del,prec);
  lam(tmp2,z,prec);
  acb_printd(tmp2,100);printf("\n");fflush(stdout);
  arb_add(acb_imagref(z),rho,del,prec);
  lam(tmp2,z,prec);
  acb_printd(tmp2,100);printf("\n");fflush(stdout);

  return 0;
}
