/*

File: convert-zeta-file.c

Created: 18 August 2008

Version: 1.0

Last Modified: 

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: 

V 1.0 Initial implementation

Build instructions:
   gcc -oconvert-zeta-file convert-zeta-file.c -O2 -lmpfi -lmpfr -lgmp

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */


/*
Reads z_file.dat created by pari program z_file.gp and fixes it, writing
zeta_file.dat
*/



#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"

#define SUCCESS 0
#define QUIT_CODE 0
#define FAILURE 1
#define IZ_FILE "z_file.dat"
#define OZ_FILE "zeta_file.dat"
/* 
Format
<N> max number of taylor terms supported
<M> number of steps up the imaginary axis
<num> numerator of size of each step
<den>
<re0>   real part of zeta(0.5+0*i)
<im0>   imag
<re1>   real part of zeta(1.5+0*i)
<im1>
...
<im(N*(M+1))> imag part of zeta(N-0.5+step*(M+1)
*/




#define TRUE (1==1)
#define FALSE (0==1)

int main()
{
  FILE *ifile,*ofile;
  int N,M,n,m;
  double step;
  mpfi_t xi;
  mpfr_t xr;

  ifile=fopen(IZ_FILE,"r");
  if(ifile==NULL)
    {
      printf("Can't open %s. Exiting.\n",IZ_FILE);
      return(SUCCESS);
    };

  ofile=fopen(OZ_FILE,"w");
  if(ofile==NULL)
    {
      printf("Can't open %s. Exiting.\n",IZ_FILE);
      return(SUCCESS);
    };

  fscanf(ifile,"%d",&N);
  if(N<=0)
    {
      printf("Invalid value for N %d. Exiting.\n",N);
      return(SUCCESS);
    };

  printf("N set to %d.\n",N);

  fscanf(ifile,"%d",&M);
  if(M<=0)
    {
      printf("Invalid value for M %d. Exiting. \n",M);
      return(SUCCESS);
    };

  printf("M set to %d.\n",M);

  fscanf(ifile,"%40lf",&step);
  if(step<=0.0)
    {
      printf("Invalid value for step %lf. Exiting.\n",step);
      return(SUCCESS);
    };

  if(fprintf(ofile,"%d\n%d\n%20.18lf\n",N,M,step)<0)
    {
      printf("Error outputting N, M and step. Exiting.\n");
      return(0);
    };

  mpfr_set_default_prec(53);
  mpfi_init(xi);
  mpfr_set_default_prec(100);  /* Pari has 28 decimal digits. This is
				  about 30 so should suffice */
  mpfr_init(xr);

  for(m=0;m<=M;m++)
    for(n=0;n<(N<<1);n++)
      {
	if(!mpfr_inp_str(xr,ifile,10,GMP_RNDN))
	  {
	    printf("Error inputting Zeta Value N:%d M:%d. Exiting.\n",n,m);
	    return(0);
	  };
	mpfi_set_fr(xi,xr);
	if(!mpfi_out_str(ofile,10,0,xi))
	  {
	    printf("Error outputing interval N:%d M:%d. Exiting.\n",n,m);
	    return(0);
	  };
	fprintf(ofile,"\n");
      };

  mpfi_clear(xi);
  mpfr_clear(xr);
  close(ifile);
  close(ofile);
  return(SUCCESS);
};
	
