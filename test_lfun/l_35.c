/*
Make up a degree 2 L-function by multiplying the quadratic characters
mod 5 and mod 7 and something else together.
*/
#include "errno.h"
#include <flint/acb_poly.h>
#include "glfunc.h"
#include "glfunc_internals.h"

#define DIGITS 20
// compute the Euler poly for p
// with L the product of non-principal characters mod 5 and 7
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec, void *param __attribute__((unused)))
{
  acb_poly_t p5;
  acb_poly_init(p5);
  acb_poly_one(p5);
  if((p%5==1)||(p%5==4))
    acb_poly_set_coeff_si(p5,1,-1);
  if((p%5==2)||(p%5==3))
    acb_poly_set_coeff_si(p5,1,1);
  acb_poly_t p7;
  acb_poly_init(p7);
  acb_poly_one(p7);
  if((p%7==1)||(p%7==2)||(p%7==4))
    acb_poly_set_coeff_si(p7,1,-1);
  if((p%7==3)||(p%7==5)||(p%7==6))
    acb_poly_set_coeff_si(p7,1,1);
  acb_poly_mul(poly,p5,p7,prec);
  acb_poly_clear(p5);
  acb_poly_clear(p7);
}

int main (int argc, char**argv)
{
  printf("Command Line:- %s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=1)
    {
      printf("Usage:- %s\n",argv[0]);
      return 0;
    }


  Lfunc_t L;
  double mus[]={0.0,1.0};


  Lerror_t ecode;
  L=Lfunc_init(2,35,0.0,mus,&ecode);
  Lfunc *lf=(Lfunc *) L;
	
  if(fatal_error(ecode))
    {
      fprintf(stderr,"Fatal error\n");
      fprint_errors(stderr,ecode);
      printf("\n");
      Lfunc_clear(L);
      fflush(stdout);
      return 0;
    }
  
  //lf->self_dual=YES; // product of two self duals is self dual

  //long int nmax=Lfunc_nmax(L);
  //printf("nmax=%ld\n",nmax);
  //Lfunc_reduce_nmax(L,3*nmax/4);

  ecode |= Lfunc_use_all_lpolys(L, lpoly_callback, NULL);

  if(fatal_error(ecode))
    {
      fprintf(stderr,"Fatal error\n");
      fprint_errors(stderr, ecode);
      printf("\n");
      fflush(stdout);
      Lfunc_clear(L);
      return 0;
    }
      
  // do the computation
  ecode|=Lfunc_compute(L);
  
  if(fatal_error(ecode))
    {
      fprintf(stderr,"Fatal error\n");
      fprint_errors(stderr,ecode);
      printf("\n");
      fflush(stdout);
      Lfunc_clear(L);
      return 0;
    }
  
  arb_t *zeros=(arb_t *)Lfunc_zeros(L,0);
  printf("First zero at:- ");
  arb_printd(zeros[0],35);
  printf("\n");
  fflush(stdout);
  
    Lplot_t *Lp=Lfunc_plot_data(L,0,32,32*5);
    for(uint64_t i=0;i<32*5;i++)
    printf("%f: %f\n",(double)i/5.0,Lp->points[i]);
  
  
  Lfunc_clear(L);
  
  // print any warnings collected along the way
  if(ecode!=ERR_SUCCESS)
    {
      fprintf(stderr,"Non-fatal error\n");
      fprint_errors(stderr,ecode);
    }
  
  return 0;
}


