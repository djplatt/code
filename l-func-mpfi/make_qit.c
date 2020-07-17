/*

File: make_qit.c

Created: 29th June 2010

Version: 1.0

Last Modified: 28th June 2010

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 0.0 Initial implementation

Build instructions: gcc make_qit.c -O2 -lmpfi -lmpfr -lgmp -lm

By: DJ Platt
    Bristol University

Copyright 2008,2009,2010.

This work is funded by the UK ESPRC. */


/*
creates a file of values of q^(it) for q=[0..100000]
t=n*dt*2^m for m=1..31
dt=gap size

  */

#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "../includes/mpfi.h"
#include "../includes/mpfi_io.h"

typedef struct {double left; double right;} int_double;
typedef struct {int_double real; int_double imag;} int_complex;

#include "../includes/upsamdefs.h"

#include "../includes/qit_struct.h"

#define dt (one_over_two_B/2.0)

void print_usage()
{
  printf("Usage:- make_qit <prec> <outfile>\n");
}

void print_int_complex_str(const char *str, const int_complex z)
{
  printf("%s [%20.18e,%20.18e]+i[%20.18e,%20.18e]\n",str,z.real.left,-z.real.right,z.imag.left,-z.imag.right);
}

#define MAX_ULP_ERR (3)
long int max_err_found=0;

long int max(long int x,long int y)
{
  if(x>=y)
    return(x);
  return(y);
}

long int check_acc1(int_double z)
{
  long int *l,*r,diff;
  double rgt;
  l=(long int *) &z.left;
  rgt=-z.right;
  r=(long int *) &rgt;
  diff=l[0]-r[0];
  if(diff<0)
    return(-diff);
  return(diff);
}

long int check_acc(int_complex z)
{
  return(max(check_acc1(z.real),check_acc1(z.imag)));
}


mpfr_t mx;

int_complex get_int_complex(mpfi_t real, mpfi_t imag)
{
  int_complex res;
  double x;
  mpfi_get_left(mx,real);
  res.real.left=mpfr_get_d(mx,GMP_RNDD);
  mpfi_get_right(mx,real);
  res.real.right=-mpfr_get_d(mx,GMP_RNDU);
  mpfi_get_left(mx,imag);
  res.imag.left=mpfr_get_d(mx,GMP_RNDD);
  mpfi_get_right(mx,imag);
  res.imag.right=-mpfr_get_d(mx,GMP_RNDU);
  max_err_found=max(max_err_found,check_acc(res));
  return(res);
}

int main(int argc, char **argv)
{

  long int prec,q,bit,i;
  qit_t qit;
  FILE *outfile;
  mpfi_t lnq,sx,cx,x;

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=3)
    {
	print_usage();
	return(1);
    }

    prec=atoi(argv[1]);
    if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    {
	print_usage();
	return(1);
    }
    mpfr_set_default_prec(prec);
    mpfr_init(mx);
    mpfi_init(lnq);
    mpfi_init(x);
    mpfi_init(sx);
    mpfi_init(cx);
    
    outfile=fopen(argv[2],"wb");
    if(!outfile)
      {
	printf("Failed to open file %s for output. Exiting.\n",argv[2]);
	return(1);
      }

    
    //q=0
    for(i=0;i<NUM_QIT_BITS;i++)
      {
	qit.bits[i].real.left=0.0;
	qit.bits[i].real.right=0.0;
	qit.bits[i].imag.left=0.0;
	qit.bits[i].imag.right=0.0;
      }
    fwrite(&qit,sizeof(qit_t),1,outfile);
    //q=1
    for(i=0;i<NUM_QIT_BITS;i++)
      {
	qit.bits[i].real.left=1.0;
	qit.bits[i].real.right=-1.0;
	qit.bits[i].imag.left=0.0;
	qit.bits[i].imag.right=-0.0;
      }
    fwrite(&qit,sizeof(qit_t),1,outfile);
    
    for(q=2;q<=MAX_Q;q++)
      {
	if((q%1000)==0)
	  printf("q=%ld.\n",q);
	mpfi_set_ui(lnq,q);
	mpfi_log(lnq,lnq);
	for(bit=0;bit<NUM_QIT_BITS;bit++)
	  {
	    mpfi_mul_d(x,lnq,dt*(1<<bit));
	    mpfi_sin(sx,x);
	    mpfi_cos(cx,x);
	    qit.bits[bit]=get_int_complex(cx,sx);
	  }
	fwrite(&qit,sizeof(qit_t),1,outfile);
      }
    fclose(outfile);
    printf("max ulp error found =%ld\n",max_err_found);
    printf("Completed.\n");
    return(0);
}
