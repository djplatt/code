#include "inttypes.h"
#include "gmp.h"
#include "acb.h"
#include "quad.h"

// use DJP's implementation of Molin's double exponential quadrature
// to integrate 16C/tlog(t)(log(t)+F(x)) as per Platt & Trudgian
// Improved bounds on Brun's constant.

//#define BEST_POSS

void F(acb_t res, const acb_t lx, int64_t prec)
{
  static bool init=false;
  static acb_t A6,A7,A8,A9,tmp1,tmp2,tmp3;
  if(!init)
    {
      init=true;
      acb_init(A6);
      acb_init(A7);
      acb_init(A8);
      acb_init(A9);
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(tmp3);
#ifdef BEST_POSS
      acb_set_ui(tmp1,927436);
      acb_div_ui(A6,tmp1,100000,prec);
#else
      acb_set_ui(tmp1,872606);
      acb_div_ui(A6,tmp1,100000,prec);
#endif
      acb_set_si(tmp1,-813199);
      acb_div_ui(A7,tmp1,100000,prec);
      //acb_set_ui(tmp1,2226754749);
      //acb_div_ui(A8,tmp1,10000,prec);
      acb_set_ui(tmp1,910879209);
      acb_div_ui(A8,tmp1,100000,prec);
      acb_set_ui(tmp1,276335978);
      acb_div_ui(A9,tmp1,10000000,prec);
#ifdef BEST_POSS
      acb_zero(A7);
      acb_zero(A8);
      acb_zero(A9);
#endif
    }
  acb_div(tmp1,A7,lx,prec);
  acb_add(tmp2,tmp1,A6,prec);
  acb_div_ui(tmp1,lx,5,prec); // s=2/5, s*log(x)/2
  acb_exp(tmp3,tmp1,prec);
  acb_mul(tmp1,tmp3,lx,prec);
  acb_div(tmp3,A8,tmp1,prec);
  acb_sub(tmp1,tmp2,tmp3,prec);
  acb_mul_2exp_si(tmp2,lx,-1);
  acb_exp(tmp3,tmp2,prec);
  acb_mul(tmp2,tmp3,lx,prec);
  acb_div(tmp3,A9,tmp2,prec);
  acb_sub(res,tmp1,tmp3,prec);
}


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
#ifdef BEST_POSS
      arb_set_ui(tmp1,927436);
      arb_div_ui(A6,tmp1,100000,prec);
#else
      arb_set_ui(tmp1,872606);
      arb_div_ui(A6,tmp1,100000,prec);
#endif
      arb_set_si(tmp1,-813199);
      arb_div_ui(A7,tmp1,100000,prec);
      //arb_set_ui(tmp1,2226754);
      //arb_div_ui(A8,tmp1,100,prec);
      arb_set_ui(tmp1,910879209);
      arb_div_ui(A8,tmp1,100000,prec);
      arb_set_ui(tmp1,2763359);
      arb_div_ui(A9,tmp1,100000,prec);
#ifdef BEST_POSS
      arb_zero(A7);
      arb_zero(A8);
      arb_zero(A9);
#endif
      printf("A6=");arb_printd(A6,20);printf("\n");
      printf("A7=");arb_printd(A7,20);printf("\n");
      printf("A8=");arb_printd(A8,20);printf("\n");
      printf("A9=");arb_printd(A9,20);printf("\n");
    }
  arb_div(tmp1,A7,lx,prec);
  arb_add(tmp2,tmp1,A6,prec);
  arb_div_ui(tmp1,lx,5,prec); // s=2/5, s*log(x)/2
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

void f(acb_t res, const acb_t x, int64_t prec)
{
  static bool init=false;
  static acb_t tmp1,tmp2,C;
  if(!init)
    {
      init=true;
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(C);
      acb_set_ui(tmp1,1320324);
      acb_div_ui(C,tmp1,1000000,prec);
      acb_mul_2exp_si(C,C,4); // 16C
    }
  F(tmp2,x,prec);
  acb_add(tmp1,tmp2,x,prec);
  acb_mul(tmp2,tmp1,x,prec);
  acb_div(res,C,tmp2,prec); // 16C/(x(x+F(exp(x))))
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

#define MAX_H ((uint64_t) 20000)
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
  arb_t res,maxd,lo,hi,tmp;
  arb_init(res);arb_init(maxd);arb_init(lo);arb_init(hi);arb_init(tmp);

  //arb_set_ui(lo,100);
  //F(hi,lo,prec);
  //printf("F(exp(100)=");arb_printd(hi,30);printf("\n");


  arb_set_ui(lo,4000000000);
  arb_mul_ui(res,lo,1000000000,prec); // 4e18
  arb_log(lo,res,prec);
  arb_set_ui(hi,50);
  arb_maxd(maxd,mysqr,lo,hi,prec,m);
  //printf("Maximum aound circle was ");arb_printd(maxd,prec);printf("\n"); 
  molin_int(res,n,f,maxd,lo,hi,prec);
  uint64_t h=75; // increase by 1.5 each time
  while(h<MAX_H)
    {
      arb_set(lo,hi);
      arb_set_ui(hi,h);
      arb_maxd(maxd,mysqr,lo,hi,prec,m);
      //printf("Maximum aound circle was ");arb_printd(maxd,prec);printf("\n"); 
      molin_int(tmp,n,f,maxd,lo,hi,prec);
      arb_add(res,res,tmp,prec);
      h=h*3;
      h=h/2;
    }
  arb_set(lo,hi);
  arb_set_ui(hi,MAX_H);
  arb_maxd(maxd,mysqr,lo,hi,prec,m);
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


  printf("B in ");arb_printd(res,prec);printf("\n");
  return(0);
}
