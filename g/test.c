/*

File: G1.6.c

Created: 15 March 2011

Version: 1.6

Last Modified: 

Dialect: C

Requires: GMP
          MPFR
          MPFI

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 1.0 Initial implementation
V 1.1 Control over maximum size of h added
      Calculation of G(1) added
      mpfi_c routines moved to mpfi_c.c
V 1.2 Revised Taylor Expansion used and taylor error added
V 1.3 Now takes zeros file in standard form
v 1.4 Taylor Error Fixed. Lambda now entered as 2^(-n)
v 1.5 rewritten for uber-accurate zeros from windowed zeta
      now uses mpfi_c.h
v 1.6 Improved Taylor error


Build instructions: gcc -oG1.5 G1.5.c -O2 -lmpfi

By: DJ Platt
    Bristol University

Copyright 2008,2009,2010,2011.

The author is funded by the UK ESPRC. */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"
#include "../includes/pi_x.h"
#include "../windowed/win_zeta.h"

int main()
{
  mpfi_c_setup(200);
  long int z=(long int) 1 << 33;
  mpfi_c_t x;
  mpfi_c_init(x);
  mpfi_set_d(x->re,4.5);
  mpfi_set_d(x->im,-1.3);
  mpfi_c_mul_ui(x,x,z);
  mpfi_c_print_str("2^33*(4.5,-1.3)=",x);
  return(0);
}
