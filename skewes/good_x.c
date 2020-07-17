/*

File: good_x.c

Created: 12th March 2012

Version: Original

Last Modified: 

Dialect: C

Requires: GMP
          MPFR
          MPFI

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 1.0 Initial implementation

Uses the first 5000 zeros of zeta to find a "good" value of
x likely to yield a skewes point.


Build instructions: gcc -osum_li sum_li.c -O2 -lmpfi

By: DJ Platt
    Bristol University

Copyright 2012

The author is funded by the UK ESPRC. */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"
#include "../windowed/win_zeta.h"

inline void next_rho(mpfi_t del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}

int main()
{
  FILE *infile;
  if(!(infile=fopen("../../zeros1/zeros_14.dat","rb")))
    {
      printf("Failed to open file zeros_14.dat for binary input. Exiting.\n");
      exit(0);
    }
  mpfr_set_default_prec(100);
  init_in_bytes();
  mpfr_t t;mpfr_init(t);
  mpfr_t pi_t0;mpfr_init(pi_t0);mpfr_const_pi(pi_t0,GMP_RNDN); // will be pi/t0
  mpfi_t del_t;mpfi_init(del_t);
  mpfr_t tmp;mpfr_init(tmp);
  mpfr_t best_sum;mpfr_init(best_sum);mpfr_set_ui(best_sum,0,GMP_RNDN);
  mpfr_t this_sum;mpfr_init(this_sum);
  unsigned long int num_its;
  fread(&num_its,sizeof(long int),1,infile);
  if(num_its!=1)
    {
      printf("There should only be one block in this file.\n");
      exit(0);
    }
  double st[2];
  unsigned long int zs[2],z;
  fread(st,sizeof(double),2,infile);
  mpfr_set_d(t,st[0],GMP_RNDN);
  fread(zs,sizeof(long int),2,infile);
  mpfr_t *zeros=(mpfr_t *)malloc(sizeof(mpfr_t)*zs[1]);
  mpfr_t *numer=(mpfr_t *)malloc(sizeof(mpfr_t)*zs[1]);
  for(z=zs[0];z<zs[1];z++)
    {
      mpfr_init(zeros[z]);
      mpfr_init(numer[z]);
      next_rho(del_t,infile);
      if(mpfi_is_zero(del_t))
	{
	  printf("Two zeros 0 apart. Exiting.\n");
	  exit(0);
	}
      mpfi_get_left(tmp,del_t);
      mpfr_add(zeros[z],t,tmp,GMP_RNDN);
      mpfr_set(t,zeros[z],GMP_RNDN);
      mpfr_sqr(numer[z],zeros[z],GMP_RNDN);
      mpfr_add_ui(numer[z],numer[z],1,GMP_RNDN);
      mpfr_div(numer[z],zeros[z],numer[z],GMP_RNDN);
    }
  fclose(infile);

  mpfr_t best_w;mpfr_init(best_w);
  mpfr_div(pi_t0,pi_t0,zeros[0],GMP_RNDN);
  long unsigned int n,best_n;
  mpfr_t w;mpfr_init(w);
  for(n=0;n<3400;n+=2)
    {
      mpfr_mul_d(w,pi_t0,1.5+n,GMP_RNDN);
      mpfr_set_ui(this_sum,0,GMP_RNDN);
      for(z=1;z<zs[1];z++)
	{
	  mpfr_mul(tmp,zeros[z],w,GMP_RNDN);
	  mpfr_sin(tmp,tmp,GMP_RNDN);
	  mpfr_mul(tmp,tmp,numer[z],GMP_RNDN);
	  mpfr_add(this_sum,this_sum,tmp,GMP_RNDN);
	}
      if(mpfr_cmp(best_sum,this_sum)>0)
	{
	  best_n=n>>1;
	  mpfr_set(best_w,w,GMP_RNDN);
	  printf("Found a new best sum with n=%lu\nSum=",n);
	  mpfr_out_str(stdout,10,0,this_sum,GMP_RNDN);
	  printf("\n");
	  mpfr_set(best_sum,this_sum,GMP_RNDN);
	}
    }
  mpfr_sub(best_sum,best_sum,numer[0],GMP_RNDN);
  mpfr_mul_2ui(best_sum,best_sum,1,GMP_RNDN);
  printf("Sum=");
  mpfr_out_str(stdout,10,0,best_sum,GMP_RNDN);
  printf("\nat n=%lu\n",best_n);
  mpfr_exp(best_w,best_w,GMP_RNDN);
  printf("This was for x=");
  mpfr_out_str(stdout,10,0,best_w,GMP_RNDN);
  printf("\n");
  return(0);
}
