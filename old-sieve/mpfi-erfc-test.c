/*

File: mpfi-erfc-test.c

Created: 26 March 2008

Version: 1.0

Last Modified: 26 March 2008

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Build instructions: for use as an include file
                    requires -lmpfi -lmpfr -lgmp (in that order!)

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */



#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "mpfi-erfc.c"

#ifndef SUCCESS
#define SUCCESS 0
#define FAILURE 1
#define QUIT_CODE 0
#endif

int main()
{
    mpfi_t ans,z;
    int i;
    mpfr_set_default_prec(53);
    if (erfc_init()!=SUCCESS)
        return(QUIT_CODE);

    mpfi_init(ans);
    mpfi_init(z);

    for(i=0;i<=40;i++)
    {
        mpfi_set_d(z,(i-20)*0.5);
        mpfi_erfc(ans,z,100);
    printf("\nErfc Called with ");mpfi_out_str(stdout,10,0,z);
    printf("\nand returned ");mpfi_out_str(stdout,10,0,ans);
    };

    erfc_clear();
    mpfi_clear(z);
    mpfi_clear(ans);
    return(QUIT_CODE);
};