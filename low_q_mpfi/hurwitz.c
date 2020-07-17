#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "../includes/mpfi.h"
#include "../includes/mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "./f_defs.h"

int main()
{
  mpfi_c_setup(300);

  mpfi_t a,b;
  mpfi_c_t res,res1,q_s;
  mpfi_init(a);mpfi_init(b);
  mpfi_c_init(res);mpfi_c_init(res1);mpfi_c_init(q_s);
  mpfi_set_ui(a,1);mpfi_set_ui(b,2);
  mpfi_div_ui(a,a,3);mpfi_div_ui(b,b,3);

  double t;
  double del=5.0/64.0;
  for(t=100000.0;t<100100.0;t+=del)
    {
      mpfi_c_pow_double_to_doubles(q_s,3.0,-0.5,-t);
      mpfi_c_hurwitz(res,0.5,t,a);
      mpfi_c_hurwitz(res1,0.5,t,b);
      mpfi_c_sub(res,res,res1);
      mpfi_c_mul(res1,res,q_s);
      printf("L(1/2+%eit)=",t);mpfi_c_print(res1);printf("\n");
    }
  return(0);
}
