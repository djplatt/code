/*

File: test-lattice.c

Created: 28 January 2009

Version: 1.0

Last Modified: 28 January 2009

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1

Build instructions: gcc -otest-lattice test-lattice.c -O2 -lmpfr -lgmp -lm

By: DJ Platt
    Bristol University

Copyright 2009.

This work is funded by the UK ESPRC. */

/*

This file checks the lattice file produced by l-func-mpfi-1.0.c against
a reference file produced by (say) Maple. It is hard wired for

zfile:- text file containing real and imaginary parts of 
        Zeta(0,0.5+j+495*I,1+(NO_GAPS-i)/NO_GAPS)
        where i runs from START_ROW to END_ROW and j from 0 to N-1
        printed with enough accuarcy to swamp double. I used %30.28e
        in Maple 9.5 running with 30 digits.

lfile:- lattice file with height t=495 as first entry, NO_GAPS gaps,
        N_OUT=N, RN_OUT=1
*/


#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"

#define zfile "/cygdrive/c/data/zeta_495_2000_2096.dat"
#define lfile "/cygdrive/c/data/foo.dat"
#define N (15)
#define NO_GAPS (4096)
#define START_ROW (2000)
#define END_ROW (2096)
#define RE (1)
#define IM (0)

double re_lattice_l[N*(NO_GAPS+1)];
double im_lattice_l[N*(NO_GAPS+1)];
double re_lattice_r[N*(NO_GAPS+1)];
double im_lattice_r[N*(NO_GAPS+1)];

void check(char typ, double x, double left, double right, int row, int col)
{
  if((x<left)||(x>right))
    {
      if(typ==RE)
	printf("Re ");
      else
	printf("Im ");
      printf("Mismatch at row %d and column %d.\n",row,col);
      printf("%20.18e vs [%20.18e,%20.18e]\n",x,left,right);
      //      exit(0);
    }
}

int main()
{

  FILE *zeta_file;
  FILE *lattice_file;
  int i,j;
  mpfr_t x;
  double x1;

  mpfr_set_default_prec(200);
  mpfr_init(x);

  if(!(zeta_file=fopen(zfile,"r")))
    {
      printf("Failed to open zeta file.Exiting.\n");
      exit(0);
    }
  if(!(lattice_file=fopen(lfile,"rb")))
    {
      printf("Failed to open lattice file.Exiting.\n");
      exit(0);
    }

  for(i=0;i<22;i++)
    fread(&j,sizeof(int),1,lattice_file);
  for(i=0;i<N*(NO_GAPS+1);i++)
    {
      fread(&re_lattice_l[i],sizeof(double),1,lattice_file);
      fread(&re_lattice_r[i],sizeof(double),1,lattice_file);
      fread(&im_lattice_l[i],sizeof(double),1,lattice_file);
      fread(&im_lattice_r[i],sizeof(double),1,lattice_file);
    };

  //printf("[%20.18e,%20.18e] [%20.18e,%20.18e]\n",re_lattice_l[0],
  // re_lattice_r[0],im_lattice_l[0],im_lattice_r[0]);
  //exit(0);
  fclose(lattice_file);

  for(i=START_ROW;i<=END_ROW;i++)
    for(j=0;j<N;j++)
      {
	mpfr_inp_str(x,zeta_file,10,GMP_RNDD);
	x1=mpfr_get_d(x,GMP_RNDD);
	check(RE,x1,re_lattice_l[(i*N)+j],re_lattice_r[(i*N)+j],i,j);

	x1=mpfr_get_d(x,GMP_RNDU);
	check(RE,x1,re_lattice_l[(i*N)+j],re_lattice_r[(i*N)+j],i,j);

	mpfr_inp_str(x,zeta_file,10,GMP_RNDD);
	if((i==4091)&&(j==4))
	  {
	    mpfr_out_str(NULL,10,0,x,GMP_RNDN);
	    printf("\n");
	  }
	x1=mpfr_get_d(x,GMP_RNDD);
	check(IM,x1,im_lattice_l[(i*N)+j],im_lattice_r[(i*N)+j],i,j);

	x1=mpfr_get_d(x,GMP_RNDU);
	check(IM,x1,im_lattice_l[(i*N)+j],im_lattice_r[(i*N)+j],i,j);
	
      }
  fclose(zeta_file);
  return(0);
}
