#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "math.h"
#include "inttypes.h"
#include "../includes/pi_x.h"

//#define DEBUG
//#define TIME

#ifdef DEBUG
#define ifdebug(x) x
#else
#define ifdebug(x)
#endif

#ifdef TIME
#include "time.h"
#endif

#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define bool int

#define INFILENAME "prime_powers.txt"

void fatal_error(const char *str)
{
  printf("%s. Exiting.\n",str);
  exit(0);
}

// convert a 128 bit unsigned int into mpfi_t
void mpfi_set_128(mpfi_ptr res, __uint128_t x)
{
  __uint64_t tmp=(__uint64_t) x; // just gets bottom 64 bits
  ifdebug(printf("bottom 64 bits = %lu\n",tmp));
  mpfi_set_ui(res,x>>64); // gets top 64 bits
  mpfi_mul_2ui(res,res,64);
  mpfi_add_ui(res,res,tmp);
}

// returns -0.5*phi(t)
void min_phi(mpfi_ptr res, mpfi_ptr t, mpfi_ptr x, mpfi_ptr root_2_lam)
{
  mpfi_div(res,t,x);
  mpfi_log(res,res);
  mpfi_mul(res,res,root_2_lam);
  mpfi_erf(res,res);
  mpfi_mul_d(res,res,0.25);
  mpfi_add_d(res,res,-0.25);
}
  
int main()
{
  FILE *infile;
  bigint xb=calc_x(),p2;
  ptype i,pow,p;
  mpfi_t x,root_2_lam,res,ph,t;
  if(!(infile=fopen(INFILENAME,"r")))
    fatal_error("Failed to open infile for non-binary read.");
  mpfi_c_setup(PREC);
  mpfi_init(root_2_lam);
  mpfi_set_d(root_2_lam,0.5);
  mpfi_sqrt(root_2_lam,root_2_lam);
  mpfi_div_d(root_2_lam,root_2_lam,LAMBDA);
  mpfi_init(x);
  mpfi_init(ph);
  mpfi_init(t);
  mpfi_init(res);
  mpfi_set_ui(res,0);
  mpfi_set_128(x,xb);
  printf("lambda=%30.28e\nx=",LAMBDA);
  print_bigint(xb);printf("\n");
  while(TRUE)
    {
      if(fscanf(infile,"%lu %lu\n",&pow,&p)==EOF)
	break;
      if(pow!=2)
	fatal_error("Only implemented for square powers.");
      p2=p;
      p2*=p2;
      mpfi_set_128(t,p2);
      min_phi(ph,t,x,root_2_lam);
      //mpfi_print_str("phi returned ",ph);
      mpfi_add(res,res,ph);
      if(p2<xb)
	mpfi_add_d(res,res,0.5);
    }
  mpfi_print_str("Sum over prime powers = ",res);
  return(0);
}

  


  
