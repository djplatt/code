#ifndef MPFI_HURWITZ
#define MPFI_HURWITZ
#include "mpfr.h"
#include "inttypes.h"

#ifdef __cplusplus
extern "C" {
#endif

  void mpfi_c_hurwitz1(mpfr_t *, mpfr_t *, mpfr_t *, uint64_t);

  void mpfi_c_lngamma1(mpfr_t *, mpfr_t *, uint64_t);

#ifdef __cplusplus
}
#endif
#endif
