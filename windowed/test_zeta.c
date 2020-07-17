#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "time.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "../includes/mpfi_fft.h"
#include "../includes/fft_defs.h"


#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define debug printf("Reached line number %d\n",__LINE__)
#define bool int

int main()
{
  mpfi_c_setup(500);

  mpfi_c_t zeta;
  mpfi_c_init(zeta);
  mpfi_t one;
  mpfi_init(one);
  mpfi_set_ui(one,1);
  double i=100000.0;

      mpfi_c_hurwitz(zeta,0.5,i,one);
      printf("zeta(1/2+%gi)=",i);mpfi_c_print(zeta);
      printf("With a relative/absolute accuracy of %d %d/%d %d\n",mpfi_rel_error(zeta->re),mpfi_rel_error(zeta->im),mpfi_abs_error(zeta->re),mpfi_abs_error(zeta->im));

  return(0);
}
