#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../G/mpfi_c.c"

#define N 11
#define G 2

#define SUCCESS 1
#define FAILURE 0


int vec_sort(mpfi_c_t *source_vec, mpfi_c_t *target_vec, 
	     unsigned int gen, unsigned int len)
{
  
  /* does a sort of vec by primitive root */
  /* question, could I do this in place?  */

  unsigned int ptr,power;

  power=1;
  for(ptr=0;ptr<(len-2);ptr++)
  {
    power=(power*gen)%len;
    /*    printf("%d ",power); */
    if(power==1)
      {
	printf("That was no primitive root!\n");
	return(FAILURE);
      };
    mpfi_c_set(target_vec[ptr],source_vec[power-1]);
  };
  mpfi_c_set(target_vec[len-2],source_vec[0]);

  return(SUCCESS);
};

int test_vec_sort()
{
  
  unsigned int ptr;
  mpfi_c_t *foo,*bar;
  mpfi_c_t nth_root;
  mpfi_t two_pi_by_n;

  foo=malloc(sizeof(mpfi_c_t)*(N-1));
  bar=malloc(sizeof(mpfi_c_t)*(N-1));

  mpfi_c_setup(53);
  mpfi_init(two_pi_by_n);
  mpfi_const_pi(two_pi_by_n);
  mpfi_mul_ui(two_pi_by_n,two_pi_by_n,2);
  mpfi_div_ui(two_pi_by_n,two_pi_by_n,N);

  mpfi_c_init(nth_root);
  mpfi_cos(nth_root->re,two_pi_by_n);
  mpfi_sin(nth_root->im,two_pi_by_n);

  for(ptr=0;ptr<N-1;ptr++)
    {
      mpfi_c_init(foo[ptr]);
      mpfi_c_init(bar[ptr]);
    };

  mpfi_c_set(foo[0],nth_root);
  for(ptr=1;ptr<N-1;ptr++)
    mpfi_c_mul(foo[ptr],foo[ptr-1],nth_root);

  if(!vec_sort(foo,bar,G,N))
    printf("Catastrophic error.\n");

  for(ptr=0;ptr<N-1;ptr++)
    mpfi_c_print(bar[ptr]);

  return(SUCCESS);
};


int main()
{
  test_vec_sort();
  return(SUCCESS);
};

