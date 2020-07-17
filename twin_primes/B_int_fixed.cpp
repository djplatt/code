#include "inttypes.h"
#include "gmp.h"
#include "acb.h"
#include "quad.h"

// use DJP's implementation of Molin's double exponential quadrature
// to integrate 16C/tlog(t)(log(t)+F(x)) as per Platt & Trudgian
// Improved bounds on Brun's constant.


arb_t alpha;



void F(arb_t res, const arb_t lx, int64_t prec)
{
  static bool init=false;
  static arb_t A6,A7,A8,A9,tmp1,tmp2,tmp3;
  if(!init)
    {
      init=true;
      arb_init(A6);
      arb_init(A7);
      arb_init(A8);
      arb_init(A9);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
#include "consts.h"
      printf("A6=");arb_printd(A6,20);printf("\n");
      printf("A7=");arb_printd(A7,20);printf("\n");
      printf("A8=");arb_printd(A8,20);printf("\n");
      printf("A9=");arb_printd(A9,20);printf("\n");
    }
  arb_div(tmp1,A7,lx,prec);
  arb_add(tmp2,tmp1,A6,prec);
  //arb_div_ui(tmp1,lx,5,prec); // s=2/5, s*log(x)/2
  arb_mul(tmp1,lx,alpha,prec);
  arb_mul_2exp_si(tmp1,tmp1,-1);
  arb_exp(tmp3,tmp1,prec);
  arb_mul(tmp1,tmp3,lx,prec);
  arb_div(tmp3,A8,tmp1,prec);
  arb_sub(tmp1,tmp2,tmp3,prec);
  arb_mul_2exp_si(tmp2,lx,-1);
  arb_exp(tmp3,tmp2,prec);
  arb_mul(tmp2,tmp3,lx,prec);
  arb_div(tmp3,A9,tmp2,prec);
  arb_sub(res,tmp1,tmp3,prec);
}


// need real version of test function for integration
void f(arb_t res, const arb_t x, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,C;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(C);
      arb_set_ui(tmp1,1320324);
      arb_div_ui(C,tmp1,1000000,prec);
      arb_mul_2exp_si(C,C,4); // 16C
    }
  F(tmp2,x,prec);
  arb_add(tmp1,tmp2,x,prec);
  arb_mul(tmp2,tmp1,x,prec);
  arb_div(res,C,tmp2,prec); // 16C/(x(x+F(exp(x))))
}


// need complex version to compute maxd
void mysqr(acb_t res, const acb_t x, int64_t prec)
{
  acb_pow_ui(res,x,2,prec);
}


void f_maxd(arb_t res, arb_t lo, arb_t hi, int64_t prec)
{
  static bool init=false;
  static arb_t A6,A7,A8,A9,tmp1,tmp2,tmp3,C;
  if(!init)
    {
      init=true;
      arb_init(A6);
      arb_init(A7);
      arb_init(A8);
      arb_init(A9);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(C);
      arb_set_ui(tmp1,1320324);
      arb_div_ui(C,tmp1,1000000,prec);
      arb_mul_2exp_si(C,C,4); // 16C

#include "consts.h"

    }
  arb_mul_ui(tmp1,lo,3,prec);
  arb_mul_2exp_si(tmp1,tmp1,-1);
  arb_sub(tmp3,tmp1,hi,prec);
  //printf("t' = ");arb_printd(tmp3,20);printf("\n");
  arb_mul_2exp_si(tmp1,tmp3,-1);
  arb_exp(tmp2,tmp1,prec);
  arb_div(res,A9,tmp2,prec);
  //printf("A9 term = ");arb_printd(res,20);printf("\n");
  arb_mul(tmp2,tmp3,alpha,prec);
  arb_mul_2exp_si(tmp2,tmp2,-1);
  arb_exp(tmp1,tmp2,prec);
  arb_div(tmp2,A8,tmp1,prec);
  //printf("A8 term = ");arb_printd(tmp2,20);printf("\n");
  arb_add(res,res,A7,prec);
  arb_div(res,res,tmp3,prec);
  arb_neg(res,res);
  //printf("neg term = ");arb_printd(res,20);printf("\n");
  arb_add(res,res,A6,prec);
  arb_add(res,res,tmp3,prec);
  arb_mul(res,res,tmp3,prec);
  arb_div(res,C,res,prec);
  if(!arb_is_positive(res))
    {
      printf("f_maxd returning non-positive interval. Exiting.\n");
      exit(0);
    }
  /*
  else
    {
      printf("f_max_d returning ");
      arb_printd(res,20);
      printf("\n");
    }
  */
}

#define MAX_H ((uint64_t) 20000)
int main(int argc, char **argv)
{
  if(argc!=5)
    {
      printf("Usage:- %s <n> <alpha n> <alpha d> <prec>.\n",argv[0]);
      printf("n is how many steps to use in quadrature.\n");
      exit(0);
    }
  int64_t n=atol(argv[1]);
  int64_t prec=atol(argv[4]);
  arb_t res,maxd,lo,hi,tmp;
  arb_init(res);arb_init(maxd);arb_init(lo);arb_init(hi);arb_init(tmp);
  arb_init(alpha);
  arb_set_ui(alpha,atol(argv[2]));
  arb_div_ui(alpha,alpha,atol(argv[3]),prec);
  printf("alpha set to ");arb_printd(alpha,20);printf("\n");
  //arb_set_ui(lo,100);
  //F(hi,lo,prec);
  //printf("F(exp(100)=");arb_printd(hi,30);printf("\n");

  
  arb_set_ui(lo,4000000000);
  arb_mul_ui(res,lo,1000000000,prec); // 4e18
  arb_log(lo,res,prec);
  arb_set_ui(hi,44);
  f_maxd(maxd,lo,hi,prec);
  //arb_maxd(maxd,f,lo,hi,prec,m);
  //printf("Maximum aound circle was ");arb_printd(maxd,prec);printf("\n"); 
  molin_int(res,n,f,maxd,lo,hi,prec);
  uint64_t h=46;
  while(h<MAX_H)
    {
      arb_set(lo,hi);
      arb_set_ui(hi,h);
      f_maxd(maxd,lo,hi,prec);
      //arb_maxd(maxd,f,lo,hi,prec,m);
      //printf("Maximum aound circle was ");arb_printd(maxd,prec);printf("\n"); 
      molin_int(tmp,n,f,maxd,lo,hi,prec);
      arb_add(res,res,tmp,prec);
      h=h*46;
      h=h/44;
    }
  arb_set(lo,hi);
  arb_set_ui(hi,MAX_H);
  f_maxd(maxd,lo,hi,prec);
  //arb_maxd(maxd,f,lo,hi,prec,m);
  //printf("Maximum aound circle was ");arb_printd(maxd,prec);printf("\n"); 
  molin_int(tmp,n,f,maxd,lo,hi,prec);
  arb_add(res,res,tmp,prec);
  

#ifndef BEST_POSS
  // now add 8/sqrt(4e18) = 2 int x_0 infty 2sqrt(x)/x^2 dx 
  arb_set_d(tmp,4e-9);
  arb_add(res,res,tmp,prec);
#endif

  // now add 16C/MAX_H
  arb_set_ui(tmp,1320324*16);
  arb_div_ui(tmp,tmp,MAX_H*1000000ul,prec);
  arb_add(res,res,tmp,prec);

  // now add B(x0)
  arb_set_ui(tmp,1840518);
  arb_div_ui(tmp,tmp,1000000,prec);
  arb_add(res,res,tmp,prec);

  // now substract 2pi(x0)/x0
  arb_set_ui(tmp,3023463123235320); // pi(x0)
  arb_div_ui(tmp,tmp,2000000000,prec);
  arb_div_ui(tmp,tmp,1000000000,prec);
  arb_sub(res,res,tmp,prec);

  // now add 8/sqrt(x0)
  arb_set_ui(tmp,2000000000); // sqrt(x0)
  arb_inv(tmp,tmp,prec);
  arb_mul_2exp_si(tmp,tmp,3);
  arb_add(res,res,tmp,prec); 

  printf("B in ");arb_printd(res,prec);printf("\n");
  return(0);
}
