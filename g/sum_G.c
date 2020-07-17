/*


By: DJ Platt
    Bristol University

Copyright 2008,2009,2010,2011.

The author is funded by the UK ESPRC. */

#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"

int main(int argc, char ** argv )
{
  FILE *infile;
  mpfi_c_setup(300);
  mpfi_t tmp;
  mpfi_c_t sum1;
  mpfi_c_t sum2;
  int i;

  if(argc!=2)
    {
      printf("Usage:- sum_G <file>\n");
      exit(0);
    }

  if(!(infile=fopen(argv[1],"r")))
    {
      printf("Failed to open file %s for input. Exiting.\n",argv[1]);
      exit(0);

    }

  mpfi_init(tmp);
  mpfi_c_init(sum1);
  mpfi_c_init(sum2);
  mpfi_c_zero(sum1);
  mpfi_c_zero(sum2);
  for(i=1;;i++)
    {
      mpfi_inp_str(tmp,infile,10); // real pt
      if(mpfi_is_zero(tmp))
	break;
      //mpfi_print_str("read ",tmp); 
      mpfi_add(sum1->re,sum1->re,tmp);
      //printf("%30.28e\n",mpfi_get_d(sum1->re));
      //mpfi_print_str("sum now ",sum1->re);
      mpfi_inp_str(tmp,infile,10); // im pt
      mpfi_add(sum1->im,sum1->im,tmp);
      mpfi_inp_str(tmp,infile,10); // re pt
      mpfi_add(sum2->re,sum2->re,tmp);
      mpfi_inp_str(tmp,infile,10); // im pt
      mpfi_add(sum2->im,sum2->im,tmp);
    }
  mpfi_c_print_str("n*G(1/2+infinity*i)-sum G(rho)-G(1/2+14i)=",sum1);
  mpfi_c_print_str("G(1/2+infinity*i)-G(1/2+14i)=",sum2);
  return(0);
}
