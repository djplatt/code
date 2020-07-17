//
// mpfi_hurwitz.c
//
// Uses arb routines to compute zeta(s,alpha)
//
#include "inttypes.h"
#include "fmpr.h"
#include "fmprb.h"
#include "fmpcb.h"

// res[4] for left,right of real/imag
// z[4] for z
// alpha[4] for alpha
fmpr_t mhl,mhr;
fmpcb_t mhcs,mhca,mhcres;
fmpz_t mha,mhb,mhexp;
mpz_t mhz;
int mpfi_c_hurwitz_first=1;
void mpfi_c_hurwitz1(mpfr_t *res, mpfr_t *z, mpfr_t *alpha, uint64_t PREC)
{
  if(mpfi_c_hurwitz_first)
    {
      fmpz_init(mha);fmpz_init(mhb);fmpz_init(mhexp);
      fmpr_init(mhl);fmpr_init(mhr);
      fmpcb_init(mhcs);
      fmpcb_init(mhca);
      fmpcb_init(mhcres);
      mpz_init(mhz);
      mpfi_c_hurwitz_first=0;
    }
  fmpr_set_mpfr(mhl,z[0]);
  fmpr_set_mpfr(mhr,z[1]);
  fmprb_set_interval_fmpr(fmpcb_realref(mhcs),mhl,mhr,PREC);
  fmpr_set_mpfr(mhl,z[2]);
  fmpr_set_mpfr(mhr,z[3]);
  fmprb_set_interval_fmpr(fmpcb_imagref(mhcs),mhl,mhr,PREC);
  fmpr_set_mpfr(mhl,alpha[0]);
  fmpr_set_mpfr(mhr,alpha[1]);
  fmprb_set_interval_fmpr(fmpcb_realref(mhca),mhl,mhr,PREC);
  fmpr_set_mpfr(mhl,alpha[2]);
  fmpr_set_mpfr(mhr,alpha[3]);
  fmprb_set_interval_fmpr(fmpcb_imagref(mhca),mhl,mhr,PREC);
  fmpcb_hurwitz_zeta(mhcres,mhcs,mhca,PREC);
  //printf("fmpcb_hur=");fmpcb_printd(mhcres,100);printf("\n");
  fmprb_get_interval_fmpz_2exp(mha,mhb,mhexp,fmpcb_realref(mhcres));
  int64_t e=fmpz_get_si(mhexp);
  fmpz_get_mpz(mhz,mha);
  mpfr_set_z_2exp(res[0],mhz,e,GMP_RNDD);
  fmpz_get_mpz(mhz,mhb);
  mpfr_set_z_2exp(res[1],mhz,e,GMP_RNDU);
  fmprb_get_interval_fmpz_2exp(mha,mhb,mhexp,fmpcb_imagref(mhcres));
  e=fmpz_get_si(mhexp);
  fmpz_get_mpz(mhz,mha);
  mpfr_set_z_2exp(res[2],mhz,e,GMP_RNDD);
  fmpz_get_mpz(mhz,mhb);
  mpfr_set_z_2exp(res[3],mhz,e,GMP_RNDU);
}

void mpfi_c_lngamma1(mpfr_t *res, mpfr_t *z, uint64_t PREC)
{
  if(mpfi_c_hurwitz_first)
    {
      fmpz_init(mha);fmpz_init(mhb);fmpz_init(mhexp);
      fmpr_init(mhl);fmpr_init(mhr);
      fmpcb_init(mhcs);
      fmpcb_init(mhca);
      fmpcb_init(mhcres);
      mpz_init(mhz);
      mpfi_c_hurwitz_first=0;
    }
  fmpr_set_mpfr(mhl,z[0]);
  fmpr_set_mpfr(mhr,z[1]);
  fmprb_set_interval_fmpr(fmpcb_realref(mhcs),mhl,mhr,PREC);
  fmpr_set_mpfr(mhl,z[2]);
  fmpr_set_mpfr(mhr,z[3]);
  fmprb_set_interval_fmpr(fmpcb_imagref(mhcs),mhl,mhr,PREC);
  fmpcb_lgamma(mhcres,mhcs,PREC);
  fmprb_get_interval_fmpz_2exp(mha,mhb,mhexp,fmpcb_realref(mhcres));
  int64_t e=fmpz_get_si(mhexp);
  fmpz_get_mpz(mhz,mha);
  mpfr_set_z_2exp(res[0],mhz,e,GMP_RNDD);
  fmpz_get_mpz(mhz,mhb);
  mpfr_set_z_2exp(res[1],mhz,e,GMP_RNDU);
  fmprb_get_interval_fmpz_2exp(mha,mhb,mhexp,fmpcb_imagref(mhcres));
  e=fmpz_get_si(mhexp);
  fmpz_get_mpz(mhz,mha);
  mpfr_set_z_2exp(res[2],mhz,e,GMP_RNDD);
  fmpz_get_mpz(mhz,mhb);
  mpfr_set_z_2exp(res[3],mhz,e,GMP_RNDU);
}
