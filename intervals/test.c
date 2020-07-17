
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"





int main()
{
  int i;
  mpfi_t foo;

  mpfr_set_default_prec(53);

  mpfi_init(foo);
  mpfi_set_d(foo,0.999999999999999);

  for(i=0;i<1000000000;i++)
    mpfi_mul(foo,foo,foo);

  printf("[%20.18e,%20.18e]\n",foo->left,foo->right);
    return(0);
}
