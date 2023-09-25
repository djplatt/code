#include "inttypes.h"
#include "acb_dirichlet.h"

#define N ((uint64_t) 20000) // this many zeros
#define PREC ((uint64_t) 3350) // 3,350 bits or about 1,000 dp

int main(int argc, char **argv)
{
  acb_t zeros[N]; // to hold N zeros, sigma+it
  uint64_t z;

  for(z=0;z<N;z++)
    acb_init(zeros[z]); // initialise zeros

  fmpz_t start;
  fmpz_init(start);
  fmpz_set_ui(start,1); // start from zero number 1

  // ask for N zeros to precision PREC starting at start
  acb_dirichlet_zeta_zeros((acb_ptr) zeros,start,N,PREC);

  // print out the imaginary part of the N'th zero
  printf("Im(rho_%lu) in ",N);
  arb_printd(acb_imagref(zeros[N-1]),100);
  printf("\n");

  // clear all data structures
  fmpz_clear(start);
  for(z=0;z<N;z++)
    acb_clear(zeros[z]);

  return 0;

}
  
