/*
The degree 2 L-function comprised of L3^2
Used to compute Lambda_3(1) and Lambda_3'(1)
Written to L3.dat
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
  acb_poly_one(poly);

  
  if(p%3==1) // (1-T)^2
    {
      acb_poly_set_coeff_si(poly,1,-2);
      acb_poly_set_coeff_si(poly,2,1);
      return;
    }
  if(p%3==2) // (1+T)^2
    {
      acb_poly_set_coeff_si(poly,1,2);
      acb_poly_set_coeff_si(poly,2,1);
      return;
    }
}


int main (int argc, char**argv)
{
  printf("Command Line:- %s",argv[0]);
  for(int i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");

  Lfunc_t L;
  double mus[]={1,1}; 
  Lerror_t ecode;

  // we have a degree 2 L-function with cond=5*7, alg=anal so
  // normalisation = 0.0
  L=Lfunc_init(2,9,0.0,mus,&ecode);
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

  acb_t ctmp,ctmp1;
  acb_init(ctmp);acb_init(ctmp1);
  ecode|=Lfunc_special_value(ctmp, L, 1, 0.0);
  printf("L(1) = ");acb_printd(ctmp,DIGITS);printf("\n");
  ecode|=Lfunc_special_value_choice(ctmp, ctmp1, L, 1, 0.0,true,true);
  if(fatal_error(ecode)) {
    fprint_errors(stderr,ecode);
    std::abort();
  }
  
  arb_t L1,L1_dash;
  arb_init(L1);arb_init(L1_dash);
  arb_sqrt(L1,acb_realref(ctmp),300);
  arb_div(L1_dash,acb_realref(ctmp1),L1,300);
  arb_mul_2exp_si(L1_dash,L1_dash,-1);

  FILE *ofile=fopen("L3.dat","w");
  arb_dump_file(ofile,L1);
  fprintf(ofile,"\n");
  arb_dump_file(ofile,L1_dash);
  fprintf(ofile,"\n");
  fclose(ofile);
  printf("%s\n",arb_dump_str(L1));
  printf("%s\n",arb_dump_str(L1_dash));
  
  printf("Lam(1) = ");acb_printd(ctmp, DIGITS);printf("\n");
  printf("Lam'(1) = ");acb_printd(ctmp1, DIGITS);printf("\n");
  printf("L3(1) = ");arb_printd(L1, DIGITS);printf("\n");
  printf("L3'(1) = ");arb_printd(L1_dash, DIGITS);printf("\n");
  
  //free memory
  arb_clear(L1);
  arb_clear(L1_dash);
  acb_clear(ctmp);acb_clear(ctmp1);
  Lfunc_clear(L);

  // print any warnings collected along the way
  fprint_errors(stderr,ecode);

  return 0;
}


