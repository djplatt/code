/*

File: make-data.c

Created: 28th June 2011

Version: 1.0

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1

V 0.0 Initial implementation

Build instructions: gcc make-data.c -O2 -lmpfr -lgmp -lm

By: DJ Platt
    Bristol University

Copyright 2008,2009,2010.

This work is funded by the UK ESPRC. */


/*

Reads input from Maple (or similar)
and converts to double prec intervals

*/

#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"

int main(int argc, char **argv)
{

  int prec,num_s,file_N,rn,no_gaps,i,j;
  double xd;
  mpfr_t x,two_gamma;
  FILE *infile,*outfile;

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=3)
    {
      printf("Usage make-data infile outfile\n");
      exit(0);
    }

    prec=100;
    if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    {
      printf("Bad precision %d\n",prec);
      exit(0);
    }
    
    mpfr_set_default_prec(prec);
    mpfr_init(x);
    mpfr_init(two_gamma);
    mpfr_const_euler(two_gamma,GMP_RNDN);
    mpfr_add(two_gamma,two_gamma,two_gamma,GMP_RNDN);

    infile=fopen(argv[1],"r");
    if(!infile)
      {
	printf("Failed to open file %s for input. Exiting.\n",argv[1]);
	exit(0);
      }

    outfile=fopen(argv[2],"wb");
    if(!infile)
      {
	printf("Failed to open file %s for binary output. Exiting.\n",argv[2]);
	exit(0);
      }

    fscanf(infile,"%d",&num_s);
    if(num_s!=1)
      {
	printf("Bad num_s %d. Exiting\n",num_s);
	exit(0);
      }
    fwrite(&num_s,sizeof(int),1,outfile);

    fscanf(infile,"%d",&file_N);
    if(file_N!=5)
      {
	printf("Bad file_N %d. Exiting\n",file_N);
	exit(0);
      }
    fwrite(&file_N,sizeof(int),1,outfile);

    fscanf(infile,"%d",&rn);
    if(rn!=2)
      {
	printf("Bad rn %d. Exiting\n",rn);
	exit(0);
      }
    fwrite(&rn,sizeof(int),1,outfile);

    fscanf(infile,"%d",&no_gaps);
    if(no_gaps!=4096)
      {
	printf("Bad no_gaps %d. Exiting\n",no_gaps);
	exit(0);
      }
    fwrite(&no_gaps,sizeof(int),1,outfile);

    for(i=0;i<=no_gaps;i++)
      {
	mpfr_inp_str(x,infile,10,GMP_RNDN);
	mpfr_sub(x,x,two_gamma,GMP_RNDN);
	xd=mpfr_get_d(x,GMP_RNDD);
	fwrite(&xd,sizeof(double),1,outfile);
	xd=mpfr_get_d(x,GMP_RNDU);
	fwrite(&xd,sizeof(double),1,outfile);
	for(j=1;j<file_N;j++)
	  {
	    mpfr_inp_str(x,infile,10,GMP_RNDN);
	    xd=mpfr_get_d(x,GMP_RNDD);
	    fwrite(&xd,sizeof(double),1,outfile);
	    xd=mpfr_get_d(x,GMP_RNDU);
	  fwrite(&xd,sizeof(double),1,outfile);
	  }
      }
    fclose(infile);
    fclose(outfile);
    exit(0);
}
