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
    Heilbronn Institute for Mathematical Research
    Bristol University

Copyright 2012

The author is a Heilbronn Research Fellow */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"
#include "../windowed/win_zeta.h"

#define N_ZEROS (4826908)

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
  mpfi_t del_t;mpfi_init(del_t);
  unsigned long int num_its,it;
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
  double *zeros=(double *)malloc(sizeof(double)*N_ZEROS);
  double *numer=(double *)malloc(sizeof(double)*N_ZEROS);
  mpfr_t tmp;mpfr_init(tmp);
  for(z=zs[0];z<zs[1];z++)
    {
      next_rho(del_t,infile);
      if(mpfi_is_zero(del_t))
	{
	  printf("Two zeros 0 apart. Exiting.\n");
	  exit(0);
	}
      mpfi_get_left(tmp,del_t);
      mpfr_add(t,t,tmp,GMP_RNDN);
      zeros[z]=mpfr_get_d(t,GMP_RNDN);
      numer[z]=zeros[z]/(zeros[z]*zeros[z]+1);
    }
  fclose(infile);

  if(!(infile=fopen("../../zeros1/zeros_5000.dat","rb")))
    {
      printf("Failed to open file zeros_14.dat for binary input. Exiting.\n");
      exit(0);
    }
  fread(&num_its,sizeof(long int),1,infile);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(zs,sizeof(long int),2,infile);
      for(z=zs[0];z<zs[1];z++)
	{
	  next_rho(del_t,infile);
	  if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_get_left(tmp,del_t);
	  mpfr_add(t,t,tmp,GMP_RNDN);
	  zeros[z]=mpfr_get_d(t,GMP_RNDN);
	  numer[z]=zeros[z]/(zeros[z]*zeros[z]+1);
	}
    }
  fclose(infile);

  if(!(infile=fopen("../../zeros1/zeros_26000.dat","rb")))
    {
      printf("Failed to open file zeros_26000.dat for binary input. Exiting.\n");
      exit(0);
    }
  fread(&num_its,sizeof(long int),1,infile);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(zs,sizeof(long int),2,infile);
      for(z=zs[0];z<zs[1];z++)
	{
	  next_rho(del_t,infile);
	  if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_get_left(tmp,del_t);
	  mpfr_add(t,t,tmp,GMP_RNDN);
	  zeros[z]=mpfr_get_d(t,GMP_RNDN);
	  numer[z]=zeros[z]/(zeros[z]*zeros[z]+1);
	}
    }
  fclose(infile);

  if(!(infile=fopen("../../zeros1/zeros_236000.dat","rb")))
    {
      printf("Failed to open file zeros_236000.dat for binary input. Exiting.\n");
      exit(0);
    }
  fread(&num_its,sizeof(long int),1,infile);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(zs,sizeof(long int),2,infile);
      for(z=zs[0];z<zs[1];z++)
	{
	  next_rho(del_t,infile);
	  if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_get_left(tmp,del_t);
	  mpfr_add(t,t,tmp,GMP_RNDN);
	  zeros[z]=mpfr_get_d(t,GMP_RNDN);
	  numer[z]=zeros[z]/(zeros[z]*zeros[z]+1);
	}
    }
  fclose(infile);

  if(!(infile=fopen("../../zeros1/zeros_446000.dat","rb")))
    {
      printf("Failed to open file zeros_4466000.dat for binary input. Exiting.\n");
      exit(0);
    }
  fread(&num_its,sizeof(long int),1,infile);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(zs,sizeof(long int),2,infile);
      for(z=zs[0];z<zs[1];z++)
	{
	  next_rho(del_t,infile);
	  if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_get_left(tmp,del_t);
	  mpfr_add(t,t,tmp,GMP_RNDN);
	  zeros[z]=mpfr_get_d(t,GMP_RNDN);
	  numer[z]=zeros[z]/(zeros[z]*zeros[z]+1);
	}
    }
  fclose(infile);



  double best_w,this_sum,best_sum;
  double pi_t0=M_PI/zeros[0];
  long unsigned int n,best_n;
  double w=pi_t0*1.5;
  pi_t0*=2;
  for(n=0;n<3400;n+=2,w+=pi_t0)
    {
      this_sum=0.0;
      for(z=1;z<zs[1];z++)
	this_sum+=numer[z]*sin(zeros[z]*w);

      if(best_sum>this_sum)
	{
	  best_n=n>>1;
	  best_w=w;
	  best_sum=this_sum;
	  printf("Found a new best sum with n=%lu\nSum=%f\n",n,(best_sum-numer[0])*2);
	}
    }
  best_sum-=numer[0];
  best_sum*=2;
  printf("Sum=%e\n",best_sum);
  printf("at n=%lu\n",best_n);
  printf("This was for x=%e\n",exp(best_w));
  return(0);
}
