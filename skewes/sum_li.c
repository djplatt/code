/*

File: sum_li.c

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

Uses the first 5000 zeros of zeta to sum 8gam/(1+4gam^2) as
an approximation to the contribution each zero will make to
the cancellation of Pi_2(x) assuming the 2tsin(tlogx) components
all line up.

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
  mpfi_c_setup(100);
  mpfi_t t;mpfi_init(t);
  mpfi_t del_t;mpfi_init(del_t);
  mpfi_t res;mpfi_init(res);mpfi_set_ui(res,0);
  mpfi_t t2;mpfi_init(t2);
  init_in_bytes();
  unsigned long int num_its;
  fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  int calced=(1==0);
  double st[2],t0;
  unsigned long int zs[2],it,z;
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(&zs[0],sizeof(long int),1,infile);
      if(st[0]==0.0)
	continue;
      if(!calced)
	{
	  calced=(1==1);
	  t0=st[0];
	}
      fread(&zs[1],sizeof(long int),1,infile);
      //printf("Processing zero %ld to %ld=%ld in total.\nt0=%f\n",zs[0]+1,zs[1],zs[1]-zs[0],t0);mpfi_c_print_str("Taylor error=",taylor_err);exit(0);
      mpfi_set_d(t,st[0]);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile);
          if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_add(t,t,del_t); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  mpfi_add(del_t,t,pm1); // not exact
	  //mpfi_print_str("t   =",t);
	  //mpfi_print_str("del_t=",del_t);
	  mpfi_sqr(t2,del_t);
	  mpfi_mul_2ui(t2,t2,2);
	  mpfi_add_ui(t2,t2,1);
	  mpfi_div(del_t,del_t,t2);
	  mpfi_add(res,res,del_t);
	  if(mpfi_cmp_d(res,0.125)>0)
	    break;
	}
    }
  mpfi_mul_2ui(res,res,3);
  mpfi_print_str("used zeros up to height ",t);
  printf("We used the first %lu zeros\n",z);
  mpfi_print_str("Sigma 8gam/(1+4gam^2) is ",res);
}
