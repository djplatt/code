//
// File: cheb.cpp
//
// Find 
//
// Input: start and end x
//        start and end pi(x)
//

#include "primesieve.hpp"
#include "mpfi.h"
#include "mpfi_io.h"
#include "stdlib.h"

#define MAX_N (1000000000000000000L) // 1e18
long unsigned int x,pi_x,y,pi_y;
mpfi_t cheb,smallest_diff,tmp1,tmp2;

void init_callback()
{
  mpfi_init(cheb);
  mpfi_set_ui(cheb,0);
  mpfi_init(smallest_diff);
  mpfi_set_ui(smallest_diff,y);
  mpfi_init(tmp1);
  mpfi_init(tmp2);
}

void callback(uint64_t prime)
{
  pi_x++;
  mpfi_set_ui(tmp1,prime);
  mpfi_log(tmp2,tmp1);
  mpfi_add(cheb,cheb,tmp2);

  mpfi_sub_ui(tmp1,cheb,prime-x);
  mpfi_neg(tmp2,tmp1);
  if(mpfi_cmp(tmp2,smallest_diff)<0)
    mpfi_set(smallest_diff,tmp2);
}

int main(int argc,char **argv)
{
  // switch rounding mode for SSE
  // this makes int_doubles work
  printf("Command line:-");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=6)
    {
      printf("Usage:- %s <start> <end> <pi(start)> <pi(end)> <prec>.\n",argv[0]);
      exit(0);
    }
  x=atol(argv[1]);
  y=atol(argv[2]);
  pi_x=atol(argv[3]);
  pi_y=atol(argv[4]);

  if((x>=y)||(pi_x>=pi_y))
    {
      printf("Expect x<y and pi(x)<pi(y).\n");
      exit(0);
    }
  if(y>MAX_N)
    {
      printf("maximum y supported=%lu. Exiting.\n",MAX_N);
      exit(0);
    }
  if((x&1)||(y&1))
    {
      printf("specify start, end and start/e to be even to ensure they are composite.\n");
      exit(0);
    }
  mpfr_set_default_prec(atol(argv[5]));
  init_callback();

  //ps1.setSieveSize(atol(argv[7])); // how many Kbyte cache

  primesieve::callback_primes(x,y,callback);
  if(pi_x!=pi_y)
    printf("Pi(%lu)should be %lu not %lu\n",y,pi_y,pi_x);
  printf("cheb(%lu)-cheb(%lu)=",y,x);
  mpfi_out_str(stdout,10,0,cheb);
  printf("\nLargest deficiency between %lu and %lu was ",x,y);
  mpfi_out_str(stdout,10,0,smallest_diff);
  return(0);
}
