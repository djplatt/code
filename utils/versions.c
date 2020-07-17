/*

File: versions.c

Created: 1 March 2011

Version: 1.0

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

Build instructions: gcc -O2 -lmpfi -lmpfr -lgmp

By: DJ Platt
    Bristol University

Copyright 2011.

This work is funded by the UK ESPRC. */


#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"

int main()
{
  printf("gmp  version = %s\n",gmp_version);
  printf("mpfr version = %s\n",mpfr_get_version());
  printf("mpfi version = %s\n",mpfi_get_version());
  return(0);
}
