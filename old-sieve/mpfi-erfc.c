/*

File: mpfi-erfc.c

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


/*

Implements an interval arithmetic version of the erfc function.

*/

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"

#ifndef SUCCESS
#define SUCCESS 0
#define FAILURE 1
#define QUIT_CODE 0
#endif

#define ERFC_MAX_K 200

mpfi_t *erfc_denom,ksum,z_2k_plus_1,zsqr,ktmp;

int erfc_init()
{

// called once to initialise global vector erfc_denom
// and other mpfi variables needed by erfc
// assumes precision of mpfr library already set

    int i;
    mpz_t kbang;
    mpfi_t half_sqrt_pi;

    erfc_denom=malloc((ERFC_MAX_K+2)*sizeof(mpfi_t));
    if(erfc_denom==NULL)
    {
        printf("Error allocationg memory for erfc denominator.\n");
        return(FAILURE);
    };


    mpfi_init(ksum);
    mpfi_init(z_2k_plus_1);
    mpfi_init(zsqr);
    mpfi_init(half_sqrt_pi);
    mpfi_init(ktmp);
    mpfi_const_pi(half_sqrt_pi);
    mpfi_sqrt(half_sqrt_pi,half_sqrt_pi);
    mpfi_div_ui(half_sqrt_pi,half_sqrt_pi,2);

    mpfi_init(erfc_denom[0]);
    mpz_init(kbang);
    mpz_set_ui(kbang,1);
    mpfi_set(erfc_denom[0],half_sqrt_pi);
    for(i=1;i<ERFC_MAX_K+2;i++)
    {
        mpfi_init(erfc_denom[i]);
        mpz_mul_ui(kbang,kbang,i);
        mpfi_mul_z(erfc_denom[i],half_sqrt_pi,kbang);
        mpfi_mul_ui(erfc_denom[i],erfc_denom[i],i+i+1);
        if((i&1)==1)
            mpfi_neg(erfc_denom[i],erfc_denom[i]);
//	printf("\n%d ",i);mpfi_out_str(stdout,10,0,erfc_denom[i]);
    };

// erfc_demon[i] <- (-1)^i i! (2i+1) sqrt(Pi)/2

    
    mpfi_clear(half_sqrt_pi);
    mpz_clear(kbang);

    return(SUCCESS);    

};

void erfc_clear()
{
    int i;
    for(i=0;i<ERFC_MAX_K+2;i++)
        mpfi_clear(erfc_denom[i]);
    free(erfc_denom);
    mpfi_clear(ksum);
    mpfi_clear(z_2k_plus_1);
    mpfi_clear(zsqr);
    mpfi_clear(ktmp);
};

void mpfi_erfc(mpfi_ptr ans, mpfi_ptr z,int K)
{
    unsigned int k;
    if(K>ERFC_MAX_K)
    {
       printf("Error, too many iterations within mpfi_erfc requested, exiting.\n");
       exit(QUIT_CODE);
    };
    mpfi_sqr(zsqr,z);
    mpfi_set(z_2k_plus_1,z);
    if(mpfi_is_neg(z))
        mpfi_neg(z_2k_plus_1,z_2k_plus_1);
    mpfi_set_ui(ksum,0);
    for(k=0;k<=K;k++)
    {
        mpfi_div(ktmp,z_2k_plus_1,erfc_denom[k]);
        mpfi_add(ksum,ksum,ktmp);
        mpfi_mul(z_2k_plus_1,z_2k_plus_1,zsqr);
    };

// sort out the truncation error here

    mpfi_set_ui(ans,1);
    if(mpfi_is_neg(z))
        mpfi_add(ans,ans,ksum);
        else
        mpfi_sub(ans,ans,ksum);
};
