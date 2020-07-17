/*


By: DJ Platt
    Bristol University

Copyright 2008,2009,2010,2011.

The author is funded by the UK ESPRC.

sum_mpfi.c

sum a file of intervals
each line is [a,b]
with a,b in R
stops reading at [0,0]
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
  int line_count=0;
  FILE *infile;
  if(argc!=3)
    {
      printf("Usage:- %s <file> <prec>\n",argv[0]);
      exit(0);
    }
  mpfi_c_setup(atol(argv[2]));
  mpfi_t tmp;
  mpfi_t sum1;
  mpfr_t diam;

  infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Failed to open %s for reading. Exiting.\n",argv[1]);
      exit(0);
    }

  mpfi_init(tmp);
  mpfi_init(sum1);
  mpfi_set_ui(sum1,0);
  mpfr_init(diam);
  while(1==1)
    {
      if(mpfi_inp_str(tmp,infile,10)==0)
	break;
      line_count++;
      if(mpfi_is_zero(tmp))
	break;
      mpfi_add(sum1,sum1,tmp);
    }
  printf("%d records read Sum = ",line_count);
  mpfi_print_str("",sum1);
  return(0);
}
