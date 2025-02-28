/*
Make up a degree 2 L-function by multiplying one of the non-real primitive
characters mod 5 and the quadratic character mod 7 together.
*/

#define DIGITS 20
#define RAW false

#define __STDC_FORMAT_MACROS
#include <chrono>
#include <cstdint>
#include <cwctype>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <flint/fmpz.h>
//#include <flint/fmpzxx.h>
#include <flint/acb_poly.h>
#include "glfunc.h"
//#include "examples_tools.h"

//using flint::fmpzxx;
using std::cout;
using std::endl;
using std::int64_t;
using std::map;
using std::ostream;
using std::size_t;
using std::vector;

// compute the Euler poly for p
// with L the product of non-principal characters mod 5 and 7
void lpoly_callback(acb_poly_t poly, uint64_t p, int d __attribute__((unused)), int64_t prec, void *param __attribute__((unused)))
{
  acb_poly_t p5;
  acb_poly_init(p5);
  acb_poly_one(p5);

  
  if((p%5==1)||(p%5==4))
    acb_poly_set_coeff_si(p5,1,-1); // odd chi
  if((p%5==2)||(p%5==3))
    acb_poly_set_coeff_si(p5,1,1);
  
  /*
  // the Euler polynomials are (1-chi(p)p^-s)^-1
  acb_t i;
  acb_init(i);
  arb_set_ui(acb_imagref(i),1); // i
  if((p%5)==1)
    acb_poly_set_coeff_si(p5,1,-1);
  if((p%5)==2)
    acb_poly_set_coeff_acb(p5,1,i);
  if((p%5)==3)
    {
      acb_neg(i,i);
      acb_poly_set_coeff_acb(p5,1,i);
      acb_neg(i,i);
    }
  if((p%5)==4)
    acb_poly_set_coeff_si(p5,1,1);
  acb_clear(i);
  */
  
  acb_poly_t p7;
  acb_poly_init(p7);
  acb_poly_one(p7);
  if((p%7==1)||(p%7==2)||(p%7==4))
    acb_poly_set_coeff_si(p7,1,-1);
  if((p%7==3)||(p%7==5)||(p%7==6)) // even chi
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

  Lfunc_t L;
  double mus[]={0,1}; // one even, one odd
  Lerror_t ecode;

  // we have a degree 2 L-function with cond=5*7, alg=anal so
  // normalisation = 0.0
  L=Lfunc_init(2,5*7,0.0,mus,&ecode);
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


  // now extract some information
  printf("Order of vanishing = %" PRIu64 "\n",Lfunc_rank(L));
  printf("Sign = ");
  acb_printd(Lfunc_sign(L),DIGITS);
  printf("\n");
  if (RAW) cout<<"RAW: "<<Lfunc_sign(L) << endl;
  printf("First non-zero Taylor coeff = ");
  arb_printd(Lfunc_Taylor(L),DIGITS);
  printf("\n");
  if (RAW) cout<<"RAW: "<<Lfunc_Taylor(L) << endl;


  acb_t ctmp,ctmp1;
  acb_init(ctmp);acb_init(ctmp1);
  ecode|=Lfunc_special_value_choice(ctmp, ctmp1, L, 1, 0.0,true,true);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("Lam(1) = ");acb_printd(ctmp, DIGITS);printf("\n");
  printf("Lam'(1) = ");acb_printd(ctmp1, DIGITS);printf("\n");
  acb_div(ctmp1,ctmp1,ctmp,100);
  printf("Lam'/Lam(1) = ");acb_printd(ctmp1,DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp << endl;
  /*
  double eps=1.0/1024.0/1024.0/1024.0/1024.0/1024.0;
  ecode|=Lfunc_special_value_choice(ctmp1, NULL, L, 1.0+eps, 0.0,true);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("Lam(1+eps) = ");acb_printd(ctmp1, DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp1 << endl;
  
  int64_t prec=200;
  acb_sub(ctmp,ctmp1,ctmp,prec);
  acb_set_d(ctmp1,eps);
  acb_div(ctmp,ctmp,ctmp1,prec);
  printf("Lam'(1) circa ");acb_printd(ctmp, DIGITS);printf("\n");
  acb_clear(ctmp);acb_clear(ctmp1);
  */

  /*
  
  ecode|=Lfunc_special_value(ctmp, L, 1.0, 1.0);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  printf("L(1.0+i) = ");acb_printd(ctmp, DIGITS);printf("\n");
  if (RAW) cout<<"RAW: "<<ctmp << endl;



  printf("First 10 zeros\n");
  // we could use Lfunc_zeros(L, 1) for the dual L-function
  arb_srcptr zeros=Lfunc_zeros(L, 0);
  for(int i  = 0; i < 10; ++i) {
    printf("Zero %d = ", i);
    arb_printd(zeros+i, DIGITS);
    printf("\n");
    if (RAW) cout<<"RAW: "<<zeros + i<< endl;
  }

  printf("Z-plot in [0, 10]:\n");
  Lplot_t *Lpp=Lfunc_plot_data(L, 0, 10.0, 20);
  int z = 0;
  double zero_double = arf_get_d(arb_midref(zeros + z), ARF_RND_NEAR);
  for(size_t k=0; k < Lpp->n_points; ++k) {
    printf("%.2f\t%.2f\t", k*Lpp->spacing , Lpp->points[k]);
    int y = 30 + int(7.5*Lpp->points[k]);
    int zero = 30;
    // assuming 60 columns
    for(int i = 0; i < 61; ++i) {
      if(i == y) {
        printf("o");
      } else if (i == zero) {
        printf("|");
      } else if ( (i > zero and i < y) or (i < zero and i > y) ) {
        printf("-");
      } else {
        printf(" ");
      }
    }
    printf("\n");
    if(k*Lpp->spacing < zero_double and (k+1)*Lpp->spacing >= zero_double){
      printf("%.2f\tzero\t", zero_double);
      for(int i = 0; i < 30; ++i)
        printf(" ");
      printf("Z\n");
      zero_double = arf_get_d(arb_midref(zeros + ++z), ARF_RND_NEAR);
    }
  }
  
  //free memory
  Lfunc_clear_plot(Lpp);
  */
  Lfunc_clear(L);

  // print any warnings collected along the way
  fprint_errors(stderr,ecode);

  return 0;
}


