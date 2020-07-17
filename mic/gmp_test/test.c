#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"

int main(int argc, char ** argv)
{

  mpz_t x;
  mpz_init(x);
  int i,lim;
  lim=atol(argv[1]);
  for(i=0;i<lim;i++)
    {
      mpz_set_ui(x,1);
      mpz_mul_2exp(x,x,100);
    }
  int res;
  res=mpz_out_str(stdout,10,x);
  printf("\nres=%d\n",res);
  return(0);
}
