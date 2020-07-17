#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"

mpfi_t nit_n,nit_logn,nit_tlogn,nit_res;
mpfr_t nit_x;

void init_nit()
{
  mpfr_set_default_prec(100);
  mpfi_init(nit_n);
  mpfi_init(nit_logn);
  mpfi_init(nit_tlogn);
  mpfi_init(nit_res);
  mpfr_init(nit_x);
}

void mpfi_nit(unsigned long int n, double t, double *rleft, double *rright, double *ileft, double *iright)
{
  printf("In mpfi_nit with n=%lu and t=%g\n",n,t);
  mpfi_set_ui(nit_n,n);
  mpfi_log(nit_logn,nit_n);
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
