/*

File: skewes_error.c

Created: 24 January 2012

Version: 1.0

Last Modified: 

Dialect: C

Requires: GMP
          MPFR
          MPFI

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 1.0 Initial implementation


Build instructions: gcc -oskewes_error skewes_error -O2 -lmpfi

By: DJ Platt
    Bristol University

Copyright 2012.

The author is funded by the UK ESPRC. */

#include "stdio.h"
#include "stdlib.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"

#define mpfi_print(str,x) {printf(str);mpfi_out_str(stdout,10,0,x);printf("\n");}
mpfi_t tmp1;

void K(mpfi_t res, mpfi_t alpha, mpfi_t y)
{
  mpfi_t(tmp);mpfi_init(tmp);
  mpfi_const_pi(res);
  mpfi_mul_2ui(res,res,1);
  mpfi_div(res,alpha,res);
  mpfi_sqrt(res,res);
  mpfi_sqr(tmp,y);
  mpfi_mul(tmp,tmp,alpha);
  mpfi_div_d(tmp,tmp,-2.0);
  mpfi_exp(tmp,tmp);
  mpfi_mul(res,res,tmp);
  mpfi_clear(tmp);
}

// S1=2/(om-et)+9.336/(om-et)^2 +(om+et)*(log 2*exp(-(om-et)/2)+2/log 2 * exp(-(om-et)/6))
void do_S1(mpfi_t res, mpfi_t eta, mpfi_t omega)
{
  mpfi_t tmp1,tmp2,tmp3,tmp4;
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  mpfi_init(tmp3);
  mpfi_init(tmp4);
  mpfi_sub(tmp1,omega,eta); // om-et
  mpfi_set_ui(tmp2,9336);
  mpfi_div_ui(res,tmp2,1000);
  mpfi_div(tmp2,res,tmp1);
  mpfi_add_ui(res,tmp2,2);
  mpfi_div(tmp2,res,tmp1); // 2/(om-et)+9.336/(om-et)^2
  //mpfi_print_str("2/()+9.336/()=",tmp2);
  mpfi_set_ui(res,2);
  mpfi_log(tmp3,res);
  mpfi_div_d(res,tmp1,-2.0);
  mpfi_exp(tmp4,res);
  mpfi_mul(res,tmp4,tmp3); // log2*exp(-(om-et)/2)
  mpfi_div_d(tmp4,tmp1,-6.0);
  mpfi_exp(tmp1,tmp4);
  mpfi_mul_2ui(tmp1,tmp1,1);
  mpfi_div(tmp4,tmp1,tmp3); // 2/log2*exp(-(om-et)/6)
  mpfi_add(tmp1,tmp4,res);
  mpfi_add(tmp3,omega,eta);
  mpfi_mul(tmp4,tmp1,tmp3);
  mpfi_add(res,tmp4,tmp2);
  mpfi_clear(tmp1);
  mpfi_clear(tmp2);
  mpfi_clear(tmp3);
  mpfi_clear(tmp4);
}

// not needed!
void do_S2(mpfi_t res, mpfi_t alpha, mpfi_t eta)
{
  mpfi_t tmp1,tmp2;
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  mpfi_sqr(res,eta); // eta^2
  mpfi_mul(tmp1,alpha,res); // alpha*eta^2
  mpfi_div_d(res,tmp1,-2.0); // -alpha*eta^2/2
  mpfi_exp(tmp1,res); // exp(-alpha*eta^2/2)
  //mpfi_print("alpha=",alpha);
  //mpfi_print("eta=",eta);
  //mpfi_print("exp(-alpha*eta^2/2)=",tmp1);
  mpfi_mul_2ui(tmp1,tmp1,1); // 2exp(-alpha*eta^2/2)
  mpfi_const_pi(res); // pi
  mpfi_mul_2ui(res,res,1); // 2 pi
  mpfi_mul(tmp2,res,alpha); // 2 pi alpha
  mpfi_sqrt(res,tmp2); // sqrt(2 pi alpha)
  mpfi_mul(tmp2,res,eta); // sqrt(2 pi alpha) eta
  mpfi_div(res,tmp1,tmp2);
  mpfi_clear(tmp1);
  mpfi_clear(tmp2);
}  

// K(eta)*(0.2+0.016/(om-et))
void do_S3(mpfi_t res, mpfi_t alpha, mpfi_t eta, mpfi_t omega)
{
  mpfi_t tmp1,tmp2;
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  K(tmp1,alpha,eta); // K(eta)
  mpfi_set_ui(tmp2,16);
  mpfi_div_ui(res,tmp2,1000); // 0.016
  mpfi_sub(tmp2,omega,eta);
  mpfi_div(res,res,tmp2); // 0.016/(omega-eta)
  mpfi_set_ui(tmp2,2);
  mpfi_div_ui(tmp2,tmp2,10); // 0.2
  mpfi_add(tmp2,tmp2,res); // 0.2+0.016/(omega-eta)
  mpfi_mul(res,tmp2,tmp1);
  mpfi_clear(tmp1);
  mpfi_clear(tmp2);
}

// exp(-T^2/2alpha)[alpha/PiT^2*log(T/2Pi)+8logT/T+4alpha/T^3][1+1/T(om-et)]
void do_S4(mpfi_t res, mpfi_t alpha, double T, mpfi_t et, mpfi_t om)
{
  mpfi_t T_,T2,pi_,tmp1,tmp2,tmp3;
  mpfi_init(T_);
  mpfi_init(T2);
  mpfi_init(pi_);
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  mpfi_init(tmp3);
  mpfi_set_d(T_,T);
  mpfi_sqr(T2,T_);
  mpfi_const_pi(pi_);
  mpfi_mul(res,pi_,T2); // pi T^2
  mpfi_div(tmp1,alpha,res); // alpha/(pi T^2)
  mpfi_div(res,T_,pi_);
  mpfi_mul_d(tmp2,res,0.5);
  mpfi_log(res,tmp2); // log(T/2pi)
  mpfi_mul(tmp3,tmp1,res);
  mpfi_log(tmp1,T_);
  mpfi_div_d(tmp2,tmp1,T);
  mpfi_mul_2ui(tmp2,tmp2,3); // 8log(T)/T
  mpfi_add(tmp1,tmp2,tmp3);
  mpfi_mul_2ui(res,alpha,2);
  mpfi_div_d(tmp2,res,T);
  mpfi_div(res,tmp2,T2);
  mpfi_add(tmp2,res,tmp1);
  mpfi_div(tmp1,T2,alpha);
  mpfi_mul_d(tmp3,tmp1,-0.5);
  mpfi_exp(tmp1,tmp3);
  mpfi_mul(res,tmp1,tmp2);
  mpfi_sub(tmp1,om,et);
  mpfi_mul_d(tmp2,tmp1,T);
  mpfi_set_ui(tmp1,1);
  mpfi_div(tmp3,tmp1,tmp2); // 1/T(om-et)
  mpfi_add_ui(tmp3,tmp3,1);
  mpfi_mul(res,res,tmp3);
  mpfi_clear(T_);
  mpfi_clear(T2);
  mpfi_clear(pi_);
  mpfi_clear(tmp1);
  mpfi_clear(tmp2);
  mpfi_clear(tmp3);
}

// S5 = 0.0032/(om-et)^2
void do_S5(mpfi_t res, mpfi_t eta, mpfi_t omega)
{
  mpfi_t tmp1,tmp2;
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  mpfi_sub(tmp2,omega,eta);
  mpfi_sqr(tmp2,tmp2);
  mpfi_set_ui(res,1);
  mpfi_div(tmp1,res,tmp2);
  mpfi_mul_2ui(res,tmp1,5); // 32/(om-et)^2
  mpfi_div_ui(res,res,1000);
  mpfi_clear(tmp1);
  mpfi_clear(tmp2);
}

// S6=AlogA exp(-A^2/2alpha +(om+et)/2)(3.2 al^-0.5+14.4et)
void do_S6(mpfi_t res, mpfi_t alpha, mpfi_t eta, mpfi_t omega, double A)
{
  mpfi_t A_,A2,tmp1,tmp2,tmp3;
  mpfi_init(A_);
  mpfi_init(A2);
  mpfi_init(tmp1);
  mpfi_init(tmp2);
  mpfi_init(tmp3);

  mpfi_set_d(A_,A);
  mpfi_sqr(A2,A_);
  mpfi_sqrt(tmp1,alpha);
  mpfi_set_ui(res,1);
  mpfi_div(tmp3,res,tmp1);
  mpfi_mul_ui(tmp2,tmp3,32); // 32 alpha^-1/2
  mpfi_div_ui(tmp2,tmp2,10); // 3.2
  mpfi_mul_ui(tmp3,eta,144);
  mpfi_div_ui(tmp3,tmp3,10);
  mpfi_add(tmp1,tmp2,tmp3); // 3.2alpha^-1/2+14.4eta
  mpfi_add(tmp2,omega,eta);
  mpfi_mul_2si(tmp2,tmp2,-1); // (om+et)/2
  mpfi_div(tmp3,A2,alpha);
  mpfi_mul_2si(tmp3,tmp3,-1); // A^2/2alpha
  mpfi_sub(res,tmp2,tmp3); // -A^2/2alpha+(om+et)/2
  //mpfi_print("-A^2/2 alpha+(om+et)/2=",res);
  mpfi_exp(tmp2,res);
  mpfi_mul(tmp3,tmp1,tmp2);
  mpfi_log(res,A_);
  mpfi_mul(tmp1,res,tmp3);
  mpfi_mul_d(res,tmp1,A);
  mpfi_clear(A_);
  mpfi_clear(A2);
  mpfi_clear(tmp1);
  mpfi_clear(tmp2);
  mpfi_clear(tmp3);
}

int main(int argc, char **argv)
{
  int i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=8)
    {
      printf("Usage:- skewes_error <prec> <om den> <om num> <eta> <alpha> <A> <T>\n");
      exit(0);
    }
  mpfr_set_default_prec(atoi(argv[1]));
  mpfi_t gamma,omega,eta,res,S1,S2,S3,S4,S5,S6,alpha;
  mpfi_init(gamma);mpfi_init(omega);mpfi_init(res);mpfi_init(eta);
  mpfi_init(alpha);
  mpfi_init(S1);
  mpfi_init(S2);
  mpfi_init(S3);
  mpfi_init(S4);
  mpfi_init(S5);
  mpfi_init(S6);
  mpfi_set_ui(omega,atol(argv[2]));
  mpfi_div_ui(omega,omega,atol(argv[3]));
  printf("omega=");mpfi_out_str(stdout,10,0,omega);printf("\n");
  double A=atof(argv[6]);
  double T=atof(argv[7]);
  printf("A=%20.18e\nT=%20.18e\n",A,T);
  mpfi_set_d(eta,atof(argv[4]));
  mpfi_print_str("eta=",eta);
  mpfi_set_ui(alpha,atol(argv[5]));
  mpfi_print_str("alpha=",alpha);

  do_S1(S1,eta,omega);
  //do_S2(S2,alpha,eta);
  mpfi_set_ui(S2,0);
  do_S3(S3,alpha,eta,omega);
  do_S4(S4,alpha,T,eta,omega);
  do_S5(S5,eta,omega);
  do_S6(S6,alpha,eta,omega,A);

  printf("S1=");mpfi_out_str(stdout,10,0,S1);printf("\n");
  printf("S2=");mpfi_out_str(stdout,10,0,S2);printf("\n");
  mpfi_print("S3=",S3);
  mpfi_print("S4=",S4);
  mpfi_print("S5=",S5);
  mpfi_print("S6=",S6);
  mpfi_add(S1,S1,S2);
  mpfi_add(S3,S3,S4);
  mpfi_add(S5,S5,S6);
  mpfi_add(S1,S1,S3);
  mpfi_add(S1,S1,S5);
  mpfi_print("Total S=",S1);
  return(0);
}

  
