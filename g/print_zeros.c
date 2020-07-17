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

double H_MAX;       /* absolute upper limit for h_max otherwise taylor won't converge */
#define SUCCESS (0)




inline void next_rho(mpfi_t del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}

int main(int argc, char **argv)
{
  FILE *infile;
  if(argc!=1)
    {
      printf("Usage print_zeros <infile>\n");
    }
  if(!(infile=fopen(argv[1],"rb")))
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[4]);
      exit(0);
    }
  mpfi_c_setup(200);
  mpfi_t t,del_t;
  mpfi_init(t);mpfi_init(del_t);
  unsigned long int num_its,it,zs[2],n,z;
  double st[2],t0;
  fread(&num_its,sizeof(long int),1,infile);
  printf("Doing %ld iterations.\n",num_its);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(&zs[0],sizeof(long int),1,infile);
      if(st[0]==0.0)
	continue;
      if(it==0)
	{
	  t0=st[0];
	  n=zs[0]+1;
	}
      fread(&zs[1],sizeof(long int),1,infile);
      printf("Processing zero %ld to %ld=%ld in total.\n",zs[0]+1,zs[1],zs[1]-zs[0]);
      mpfi_set_d(t,st[0]);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile);
	  //mpfi_print_str("del_t=",del_t);
          if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  if(mpfi_is_neg(del_t))
	    {
	      printf("Negative del_t. Exiting.\n");
	      exit(0);
	    }
	  mpfi_add(t,t,del_t); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  printf("zero number %lu is at ",z);mpfi_print_str("",t);
	}
      // now do the G(end)-G(s_last)
      printf("st[1]=%12.0f\n",st[1]);
      mpfi_sub_d(del_t,t,st[1]);
      mpfi_neg(del_t,del_t);
      mpfi_print_str("Final delta t = ",del_t);
      mpfi_sub(del_t,del_t,del_t);
      //mpfi_print_str("t=",t);
      if(!mpfi_is_zero(del_t))
	printf("t was not exact.\n");
      mpfi_set_d(t,st[1]);
    }
  return(0);
}
