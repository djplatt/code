#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"

int main(int argc, char ** argv)
{
  if(argc!=3)
    {
      printf("Usage:- euler <prec> <t>\n");
      exit(0);
    }

  double t=atof(argv[2]);
  printf("t=%f\n",t);
  mpfi_c_setup(atoi(argv[1]));
  
  mpfi_c_t res,lng;
  mpfi_t one;

  mpfi_c_init(res);
  mpfi_c_init(lng);
  mpfi_init(one);
  mpfi_set_ui(one,1);
  mpfi_c_hurwitz(res,0.5,t,one);
  mpfi_c_print_str("Zeta(1/2+It) = ",res);
  mpfi_c_lngamma_hi(lng,0.25,t/2.0);
  mpfi_const_pi(one);
  mpfi_log(one,one);
  mpfi_mul_d(one,one,-t/2);
  mpfi_add(one,lng->im,one);
  mpfi_cos(lng->re,one);
  mpfi_sin(lng->im,one);
  mpfi_c_mul(res,lng,res);
  mpfi_c_print_str("Z(t) = ",res);

  return(1);
}
