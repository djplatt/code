/*
   compute 1/2 int K(u-omega)li(u/2)exp(-u/2)u du for
   u from omega-eta to omega+eta

   where li(x)=Ei(log(x))
   K(y)=sqrt(alpha/(2Pi)) exp(-alpha y^2/2)

*/
#include "acb.h"
#include "../trace/quad.h"

// Li(exp(u))/exp(u)
void li_re(arb_t res, const arb_t u, uint64_t n, int64_t prec)
{
  static bool init=false;
  static arb_t rhou,rhouk,ctmp1,ctmp2;
  if(!init)
    {
      init=true;
      arb_init(rhou);
      arb_init(rhouk);
      arb_init(ctmp1);
      arb_init(ctmp2);
    }

  uint64_t kfac=1;

  arb_inv(rhou,u,prec); // rhou=1/u
  arb_set(rhouk,rhou); // 1/u^1
  arb_set(ctmp1,rhou);
  for(uint64_t k=2;k<=n;k++)
    {
      kfac*=k-1;
      arb_mul(rhouk,rhouk,rhou,prec); // 1/(u)^k
      arb_mul_ui(ctmp2,rhouk,kfac,prec); // (k-1)!/(u)^k
      arb_add(ctmp1,ctmp1,ctmp2,prec);
    }

  arb_set(res,ctmp1);
}


// Re (Li(exp(u(1/2 + i gam)))+Li(exp(u(1/2-i gam))))/exp(u/2)
void li_star(arb_t res, const arb_t u, const arb_t gam, uint64_t n, int64_t prec)
{
  static bool init=false;
  static arb_t gamu,tmp1,tmp2;
  static acb_t rho,ctmp1,ctmp2,ctmp3,rhou,rhouk;
  if(!init)
    {
      init=true;
      arb_init(gamu);
      arb_init(tmp1);
      arb_init(tmp2);
      acb_init(rho);
      acb_init(ctmp1);
      acb_init(ctmp2);
      acb_init(ctmp3);
      acb_init(rhou);
      acb_init(rhouk);
    }
  uint64_t kfac=1;

  arb_set_d(acb_realref(rho),0.5);
  arb_set(acb_imagref(rho),gam); // rho =1/2+i gam

  acb_mul_arb(ctmp1,rho,u,prec);
  acb_inv(rhou,ctmp1,prec); // rhou=1/rho u
  acb_set(rhouk,rhou); // 1/(rho u)^1
  acb_set(ctmp1,rhou);
  for(uint64_t k=2;k<=n;k++)
    {
      kfac*=k-1;
      acb_mul(rhouk,rhouk,rhou,prec); // 1/(rho u)^k
      acb_mul_ui(ctmp2,rhouk,kfac,prec); // (k-1)!/(rho u)^k
      acb_add(ctmp1,ctmp1,ctmp2,prec);
    }

  arb_zero(acb_realref(ctmp2));
  arb_set(acb_imagref(ctmp2),gamu);
  acb_exp(ctmp3,ctmp2,prec);

  acb_mul(ctmp2,ctmp1,ctmp3,prec);

  arb_mul(gamu,gam,u,prec);
  arb_pow_ui(tmp1,gamu,n+1,prec);
  arb_inv(tmp2,tmp1,prec);
  arb_mul_ui(tmp1,tmp2,kfac*n,prec);
  arb_add_error(acb_realref(ctmp2),tmp1);
  arb_mul_2exp_si(res,acb_realref(ctmp2),1);
}

double alpha;
arb_t arb_omega,min_al_by_2,sqrt_al_by_2pi;

// Gaussian Kernel K(y-omega) as per S. T. D. 2.3
void K(arb_t res, const arb_t y, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
    }
  //printf("K called with u = ");arb_printd(y,20);printf("\n");
  arb_sub(tmp2,y,arb_omega,prec);
  //printf("Doing exp with y = ");arb_printd(tmp2,20);printf("\n");
  arb_mul(tmp1,tmp2,tmp2,prec); // y^2
  arb_mul(tmp2,tmp1,min_al_by_2,prec); // -y^2 alpha/2
  arb_exp(tmp1,tmp2,prec);
  arb_mul(res,tmp1,sqrt_al_by_2pi,prec);
  //printf("K returning ");arb_printd(res,20);printf("\n");
}

//u*Li(exp(u/2))*exp(-u/2)
void li_bit_re(arb_t res, const arb_t u, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_mul_2exp_si(tmp,u,-1);
  li_re(res,tmp,10,prec);
  arb_mul(res,res,u,prec);
}

int main()
{
  int64_t prec=200,n=1<<16;
  alpha=1.0e16;

  arb_init(arb_omega);
  arb_set_ui(arb_omega,727951335426);
  arb_div_ui(arb_omega,arb_omega,1000000000,prec);

  arb_t arb_eta;
  arb_init(arb_eta);
  arb_set_ui(arb_eta,1);
  arb_div_ui(arb_eta,arb_eta,1000000,prec);

  arb_t res,maxd,lo,hi,tmp1,tmp2;
  arb_init(res);
  arb_init(maxd);
  arb_init(lo);
  arb_init(hi);
  arb_init(tmp1);
  arb_init(tmp2);

  arb_sub(lo,arb_omega,arb_eta,prec); // left pt of integration
  arb_add(hi,arb_omega,arb_eta,prec); // right pt

  arb_init(min_al_by_2); // constants used by K
  arb_set_d(min_al_by_2,-alpha/2.0);
  arb_init(sqrt_al_by_2pi);
  arb_const_pi(tmp1,prec);
  arb_div(tmp2,min_al_by_2,tmp1,prec); // -alpha/2Pi
  //printf("-al/2 = ");arb_printd(min_al_by_2,20);printf("\n");
  arb_neg(tmp2,tmp2);
  arb_sqrt(sqrt_al_by_2pi,tmp2,prec);
  //printf("sqrt(al/2Pi) = ");arb_printd(sqrt_al_by_2pi,20);printf("\n");


  // compute max K((u-omega))Li(exp(u/2))exp(-u/2)*u
  // on omega+2*eta*exp(2*Pi*I*t), t in [0,1)
  // < 4*K(2i eta)
  arb_set_d(maxd,alpha);
  arb_mul(tmp1,maxd,arb_eta,prec);
  arb_mul(maxd,tmp1,arb_eta,prec);
  arb_mul_2exp_si(maxd,maxd,1);
  arb_exp(tmp1,maxd,prec);
  arb_mul(maxd,tmp1,sqrt_al_by_2pi,prec); // K(2 i eta)

  arb_mul_2exp_si(maxd,maxd,2); // max Ei(u/2)exp(-u/2)u is at theta=Pi and is less than 4

  molin_int(res,n,K,li_bit_re,maxd,lo,hi,prec);
  arb_mul_2exp_si(res,res,-1);
  arb_printd(res,20);printf("\n");

  return(0);
}
