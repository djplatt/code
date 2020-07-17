//
// Example from Moore quoted in my thesis
//
// try running at 21,52,112,200 bits
//
// 333.75b^6+a^2(11a^2b^2-b^6-121b^4-2)+5.5b^8+a/2b
//
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "stdlib.h"

int main(int argc, char ** argv)
{
  if(argc!=2)
    {
      printf("Usage:- moore <prec in bits>\n");
      exit(0);
    }
  long int PREC=atol(argv[1]);
  mpfr_set_default_prec(PREC);

  mpfr_t a,a_2,b,b_2,b_4,b_6,b_8,a_2b,temp1,temp2,res;
  mpfr_init(a);mpfr_init(a_2);mpfr_init(b);mpfr_init(b_2);
  mpfr_init(b_4);mpfr_init(b_6);mpfr_init(b_8);
  mpfr_init(a_2b);mpfr_init(temp1);mpfr_init(temp2);
  mpfr_init(res);

  mpfr_set_d(a,77617.0,GMP_RNDN);
  mpfr_set_d(b,33096.0,GMP_RNDN);

  mpfr_mul(a_2,a,a,GMP_RNDN);
  mpfr_mul(b_2,b,b,GMP_RNDN);
  mpfr_mul(b_4,b_2,b_2,GMP_RNDN);
  mpfr_mul(b_6,b_4,b_2,GMP_RNDN);
  mpfr_mul(b_8,b_4,b_4,GMP_RNDN);
  mpfr_div(a_2b,a,b,GMP_RNDN);
  mpfr_div_ui(a_2b,a_2b,2,GMP_RNDN);

  mpfr_mul_d(res,b_6,333.75,GMP_RNDN);

  mpfr_mul_ui(temp1,a_2,11,GMP_RNDN);
  mpfr_mul(temp1,temp1,b_2,GMP_RNDN);
  mpfr_sub(temp1,temp1,b_6,GMP_RNDN);
  mpfr_mul_ui(temp2,b_4,121,GMP_RNDN);
  mpfr_sub(temp1,temp1,temp2,GMP_RNDN);
  mpfr_sub_ui(temp1,temp1,2,GMP_RNDN);
  mpfr_mul(temp1,temp1,a_2,GMP_RNDN);
  mpfr_add(res,res,temp1,GMP_RNDN);

  mpfr_mul_d(temp1,b_8,5.5,GMP_RNDN);
  mpfr_add(res,res,temp1,GMP_RNDN);

  mpfr_add(res,res,a_2b,GMP_RNDN);

  printf("Result at %lu bits is ",(long unsigned int) PREC);
  mpfr_out_str(stdout,10,50,res,GMP_RNDN);
  printf("\n");
  return(0);
}
