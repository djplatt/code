/*
Make up a degree 2 L-function by multiplying the quadratic character
mod 5 and 7 together.
*/

#define DIGITS 20
#include <acb_poly.h>
#include "glfunc.h"

// compute the Euler poly for p
// with L the product of non-principal characters mod 5 and 7
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec, void *param __attribute__((unused)))
{
  acb_poly_t p5;
  acb_poly_init(p5);
  acb_poly_one(p5);
  /*  if((p%5==1)||(p%5==4))
    acb_poly_set_coeff_si(p5,1,-1);
  if((p%5==2)||(p%5==3))
    acb_poly_set_coeff_si(p5,1,1);
  */
  if((p%5)==1)
    acb_poly_set_coeff_si(p5,1,-1);
  if((p%5)==2)
    acb_poly_set_coeff_si(p5,1,1);
  if((p%5)==3)
    acb_poly_set_coeff_si(p5,1,1);
  if((p%5)==4)
    acb_poly_set_coeff_si(p5,1,-1);

  acb_poly_t p3;
  acb_poly_init(p3);
  acb_poly_one(p3);
  if(p%3==1)
    acb_poly_set_coeff_si(p3,1,-1);
  if(p%3==2)
    acb_poly_set_coeff_si(p3,1,1);
  acb_poly_mul(poly,p5,p3,prec);
  acb_poly_clear(p5);
  acb_poly_clear(p3);
  if(p<20)
    {
      printf("poly for prime %d was ",p);
      acb_poly_printd(poly,20);
      printf("\n");
    }
}


int main (int argc, char**argv)
{
  printf("Command Line:- %s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");

  Lfunc_t L;
  double mus[]={0,1};
  Lerror_t ecode;

  // we have a degree 2 L-function with alg=anal so normalisation = 0.0
  L=Lfunc_init(2,3*5,0.0,mus,&ecode);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr,ecode);
    return 0;
  }

  ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, NULL);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr, ecode);
    return 0;
  }

  // do the computation
  ecode|=Lfunc_compute(L);
  if(fatal_error(ecode))
  {
    fprint_errors(stderr,ecode);
    return 0;
  }


  
  Lfunc_clear(L);

  // print any warnings collected along the way
  fprint_errors(stderr,ecode);

  return 0;
}


