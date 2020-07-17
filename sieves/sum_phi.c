/*


By: DJ Platt
    Bristol University

Copyright 2008,2009,2010,2011.

The author is funded by the UK ESPRC.

sum_phi.c

sum output from phi2.0.c 

output from mpfi_print is quite right
needs to be [ number , number ]
need a zero record at the end of the input file
 */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"

//#define infilename ("all_phi.txt")

int main(int argc, char ** argv)
{
  FILE *infile;
  if(argc!=2)
    exit(0);
  mpfi_c_setup(300);
  mpfi_t tmp;
  mpfi_t sum1;
  mpfr_t diam;

  infile=fopen(argv[1],"r");

  mpfi_init(tmp);
  mpfi_init(sum1);
  mpfi_set_ui(sum1,0);
  mpfr_init(diam);
  while(1==1)
    {
      mpfi_inp_str(tmp,infile,10); // real pt
      if(mpfi_is_zero(tmp))
	break;
      //mpfi_print_str("read=",tmp);
      //mpfi_diam_abs(diam,tmp);
      //printf("absolute diameter in term=%6.5e\n",mpfr_get_d(diam,GMP_RNDN));
      mpfi_add(sum1,sum1,tmp);
      //mpfi_diam_abs(diam,sum1);
      //printf("absolute diameter in sum=%6.5e\n",mpfr_get_d(diam,GMP_RNDN));
      //printf("absolute error in sum =%d\n",mpfi_abs_error(sum1));
    }
  mpfi_print_str("Sigma phi=",sum1);
  return(0);
}
