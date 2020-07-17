// modified from /projects/Zeros../buethe../source
// reads djp and buethe zeros one at a time and makes sure the latter
// contains the former
//

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <mpfr.h>
#include "arb.h"

#include "zero_stream.h"
#include "if.h"
#include "asm/config.h"
#include "djp_zeros.c"

static void Usage()
{
  printf("zeta-zeros N M zero-files-file\n");
  printf("prints M consecutive zeros of the Riemann zeta function\nbeginning with the N-th zero\n"); 
  exit(1);
}

int main(int argc, char **argv)
{
  char *zffn;
  if(argc!=3)
    {
      printf("Usage: %s <buthe list> <djp list>.\n",argv[0]);
      return 0;
    }
  zffn=argv[1];
  djp_fnames=argv[2];
  mpfr_t zero;
  mpfr_set_default_prec(200);
  mpfr_init(zero);
  if(zs_init(zffn)<0){
	complain("error reading zero-files-file\n");
  }
  long unsigned int n=1;
  arb_t djp_zero,buthe_zero;
  arb_init(djp_zero);
  arb_init(buthe_zero);
  uint64_t error_count=0;
#define MAX_ERRORS (1000)
  for(;;n++){
    if(djp_next_zero(djp_zero,200)<0)
      return 0;
    if(zs_get_next_zero(zero)<0)
  	complain("out of zeros\n");
    arb_set_interval_mpfr(buthe_zero,zero,zero,200);
    arb_add_error_2exp_si(buthe_zero,-65);
    /*
    if(n==1)
      {
	printf("%lu: Buthe zero ",n);arb_print(buthe_zero);
	printf(" djp zero ");arb_print(djp_zero);
	printf("\n");
      }
    */
    if(!arb_contains(buthe_zero,djp_zero))
      {
	error_count++;
	printf("%lu: Buthe zero ",n);arb_print(buthe_zero);
	printf(" does not contain djp zero ");arb_print(djp_zero);
	printf("\n");
	if(error_count>=MAX_ERRORS)
	  return 0;
      }


  }
  return 0;
}
