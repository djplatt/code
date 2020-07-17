//
// test_hurwitz.c
//
// Uses arb routines to compute zeta(s,alpha)
//
#include "inttypes.h"
#include "fmpr.h"
#include "fmprb.h"
#include "fmpcb.h"

#define PREC (300)
int main()
{
  fmpcb_t z,alpha,res;
  fmpcb_init(z);fmpcb_init(alpha);fmpcb_init(res);
  fmprb_set_ui(fmpcb_realref(z),1);
  fmprb_div_ui(fmpcb_realref(z),fmpcb_realref(z),2,PREC);
  fmpz_t mz;
  fmpz_init(mz);
  fmpz_set_str(mz,"32994123041857409634601541066702486179099146220380442395407908848237177323881041957065463066",10);
  fmprb_set_fmpz(fmpcb_imagref(z),mz);
  fmpz_set_str(mz,"10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000",10);
  fmprb_div_fmpz(fmpcb_imagref(z),fmpcb_imagref(z),mz,PREC);
  fmprb_set_ui(fmpcb_realref(alpha),4);
  fmprb_div_ui(fmpcb_realref(alpha),fmpcb_realref(alpha),100,PREC);
  fmprb_set_ui(fmpcb_imagref(alpha),0);
  fmpcb_hurwitz_zeta(res,z,alpha,PREC);
  fmpcb_printd(res,100);
  return(0);
}
