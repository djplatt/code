#include "arb.h"
#include "acb_hypgeom.h"
#include "inttypes.h"

void arb_hypgeom_ei(arb_t res, arb_t z, int64_t prec)
{
  static bool init=false;
  static acb_t rr,zz;
  if(!init)
    {
      init=true;
      acb_init(rr);
      acb_init(zz);
    }
  arb_set(acb_realref(zz),z);
  acb_hypgeom_ei(rr,zz,prec);
  arb_set(res,acb_realref(rr));
  //printf("li returning ");arb_printd(res,20);printf("\n");
}

void pi_error(arb_t res, arb_t lx, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3,c1,c2,c3;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(c1);
      arb_init(c2);
      arb_init(c3);
      arb_set_ui(tmp1,62);
      arb_div_ui(c1,tmp1,10,prec);
      arb_set_ui(tmp1,1009);
      arb_div_ui(c2,tmp1,1000,prec);
      arb_set_ui(tmp1,83742);
      arb_div_ui(c3,tmp1,10000,prec);
      arb_neg(c3,c3);
    }
  arb_pow(tmp1,lx,c2,prec);
  arb_mul(tmp2,tmp1,c1,prec);
  arb_sqrt(tmp1,lx,prec);
  arb_mul(tmp3,tmp1,c3,prec);
  arb_exp(tmp1,tmp3,prec);
  arb_mul(tmp3,tmp2,tmp1,prec);
  arb_exp(tmp1,lx,prec);
  arb_mul(res,tmp3,tmp1,prec);
  //printf("pi_error returning ");arb_printd(res,20);printf("\n");
} 

void pip(arb_t res, arb_t lx, int64_t prec)
{
  static bool init=false;
  static arb_t err,li;
  if(!init)
    {
      init=true;
      arb_init(err);
      arb_init(li);
    }
  pi_error(err,lx,prec);
  printf("pip_error = ");arb_printd(err,20);printf("\n");
  arb_hypgeom_ei(li,lx,prec);
  printf("Ei(lx) = ");arb_printd(li,20);printf("\n");
  arb_add(res,li,err,prec);
}

void pim(arb_t res, arb_t lx, int64_t prec)
{
  static bool init=false;
  static arb_t err,li;
  if(!init)
    {
      init=true;
      arb_init(err);
      arb_init(li);
    }
  pi_error(err,lx,prec);
  printf("pim_error = ");arb_printd(err,20);printf("\n");
  arb_hypgeom_ei(li,lx,prec);
  printf("Ei(lx-1) = ");arb_printd(li,20);printf("\n");
  arb_sub(res,li,err,prec);
}

void g(arb_t res, arb_t lx, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3,tmp4;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp4);
    }
  pip(tmp1,lx,prec);
  printf("pip(lx) = ");arb_printd(tmp1,20);printf("\n");
  arb_mul(tmp2,tmp1,tmp1,prec); // pi^2(x)
  printf("pip^2(lx) = ");arb_printd(tmp2,20);printf("\n");
  arb_sub_ui(tmp1,lx,1,prec);
  pim(tmp3,tmp1,prec);
  printf("pim(lx-1) = ");arb_printd(tmp3,20);printf("\n");
  arb_div(tmp1,tmp3,lx,prec);
  arb_add_ui(tmp3,lx,1,prec);
  arb_exp(tmp4,tmp3,prec);
  arb_mul(tmp3,tmp4,tmp1,prec);
  printf("ex/lx pim(lx-1) = ");arb_printd(tmp3,20);printf("\n");
  arb_sub(res,tmp2,tmp3,prec);
}

int main()
{
  arb_t lx,err,res;
  arb_init(lx);
  arb_init(err);
  arb_init(res);
  arb_set_ui(lx,4000);
  arb_set_ui(err,1);
  arb_mul_2exp_si(err,err,-30);
  arb_add_error(lx,err);
  g(res,lx,1000);
  arb_printd(res,30);
  printf("\n");
  return 0;
}

