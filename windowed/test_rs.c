//
// test_rs.c
//
// test new Riemann Siegel routine
//
// Vesrion 1.0 Initial implementation
//
// Created: 2 November 2010
// Last Modified: 2 November 2010
//
// DJ Platt
// University of Bristol
//

#include "stdlib.h"
#include "stdio.h"
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

void print_usage()
{
  printf("Usage:- test_rs <prec> <t0> <# points>. Exiting.\n");
  exit(0);
}

#define TO_DEL ((double) 41.0/1024.0)

int main(int argc, char **argv)
{
  int prec;
  double t0,t,nu;
  int num_t,n;
  mpfi_t *logs,*sqrts;
  mpfi_t res;

  if(argc!=4)
    print_usage();

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    print_usage();
  printf("Running at %d bits of precision.\n",prec);

  t0=atof(argv[2]);
  if(t0<1.0e6)
    {
      printf("Need t0>10^6 for lngamma and r_s. Exiting.\n");
      exit(0);
    }

  num_t=atoi(argv[3]);
  if(num_t<=0)
    print_usage();
    
  mpfi_c_setup(prec);

  nu=sqrt((t0+TO_DEL*(num_t-1))/2/M_PI);
  logs=(mpfi_t *) malloc(sizeof(mpfi_t)*nu);
  if(!logs)
    {
      printf("Failed to allocate memory for logs. Exiting.\n");
      exit(0);
    }
  sqrts=(mpfi_t *) malloc(sizeof(mpfi_t)*nu);
  if(!sqrts)
    {
      printf("Failed to allocate memory for logs. Exiting.\n");
      exit(0);
    }
  mpfi_init(res);

  for(n=0;n<nu;n++)
    {
      mpfi_init(logs[n]);
      mpfi_init(sqrts[n]);
      mpfi_set_ui(sqrts[n],1);
      mpfi_set_ui(logs[n],n+1);
      mpfi_sqrt(res,logs[n]);
      mpfi_div(sqrts[n],sqrts[n],res);
      mpfi_log(logs[n],logs[n]);
    }


  for(n=0,t=t0;n<num_t;n++,t+=TO_DEL)
    {
      mpfi_zeta_rs(res,t,logs,sqrts);
      printf("Z(%f)=",t);mpfi_printn(res,40);
    }
  
}
