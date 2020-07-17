
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "math.h"
#include "float.h"


int main()
{
  mpfi_t x;
  int i;

  mpfr_set_default_prec(53);
  mpfi_init(x);

  mpfi_interv_d(x,1.0,nextafter(1.0,DBL_MAX));

  for(i=0;i<50;i++)
    mpfi_mul(x,x,x);

  mpfi_out_str(stdout,10,0,x);
  return(0);
};
