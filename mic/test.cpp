#include <stdlib.h>
#include <stdio.h>
#include "gmp.h"

int main(int argc, char ** argv)
{
  if(argc!=2)
    exit(0);

  mpz_t x;
  mpz_init(x);
  mpz_set_ui(x,1);
  mpz_mul_2exp(x,x,100);
  mpz_out_str(stdout,10,x);
}
