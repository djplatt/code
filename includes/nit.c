#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"

mpfi_t nit_n,nit_logn,nit_tlogn,nit_res;
mpfr_t nit_x;

void nit_init()
{
  //printf("In init_nit.\n");
  mpfr_set_default_prec(100);
  mpfi_init(nit_n);
  mpfi_init(nit_logn);
  mpfi_init(nit_tlogn);
  mpfi_init(nit_res);
  mpfr_init(nit_x);
}

unsigned long int last_n=0;
void mpfi_nit(unsigned long int n, double t, double *rleft, double *rright, double *ileft, double *iright)
{
  mpfi_set_ui(nit_n,n);
  if(n!=last_n)
    {
      mpfi_log(nit_logn,nit_n);
      last_n=n;
    }
  mpfi_mul_d(nit_tlogn,nit_logn,t);
  mpfi_cos(nit_res,nit_tlogn);
  mpfi_get_left(nit_x,nit_res);
  rleft[0]=mpfr_get_d(nit_x,GMP_RNDD);
  mpfi_get_right(nit_x,nit_res);
  rright[0]=-mpfr_get_d(nit_x,GMP_RNDU);
  mpfi_sin(nit_res,nit_tlogn);
  mpfi_get_left(nit_x,nit_res);
  ileft[0]=mpfr_get_d(nit_x,GMP_RNDD);
  mpfi_get_right(nit_x,nit_res);
  iright[0]=-mpfr_get_d(nit_x,GMP_RNDU);
}
