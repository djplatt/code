#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"

#include "mpfi.h"
#include "mpfi_io.h"
#include "../G/mpfi_c.c"
#define PREC 53
#define LIM 1000000

int main()
{
  int i;
  mpfi_c_t z1,z2,*z3;
  mpfi_c_setup(PREC);

  mpfi_c_init(z1);
  mpfi_c_init(z2);

  mpfi_set_d(z1->re,10.0);
  mpfi_set_d(z1->im,6.2);
  mpfi_set_d(z2->re,4.0);
  mpfi_set_d(z2->im,-1.5);

  z3=malloc(sizeof(mpfi_c_t)*LIM);

  for(i=0;i<LIM;i++)
    {
      mpfi_c_init(z3[i]);
      mpfi_c_add(z3[i],z1,z2);
      mpfi_c_mul(z3[i],z3[i],z1);
      mpfi_c_div(z3[i],z3[i],z1);
      mpfi_c_exp(z3[i],z3[i]);
    };


  return(0);
};
