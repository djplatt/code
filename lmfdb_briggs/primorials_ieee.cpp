/*

Compute successive primorials
and check Robin's inequality 
sigma(n)/n-exp(gamma)log log n

Uses int_double package

D.J.Platt 2019
*/
  
#include "stdio.h"
#include "stdlib.h"
#include "primesieve.h" // use Kim Walisch's primesieve
#undef ulong
#include "pari/pari.h"
#include "pari.c"
#undef ulong
#include "../includes/int_double14.2.h"


#define STACK_LEN (128)
int_double lhs[STACK_LEN];
uint64_t lhs_ptr=0;
int_double rhs[STACK_LEN];
uint64_t rhs_ptr=0;

int_double add_item(int_double *stack, int_double item, uint64_t k, uint64_t *ptr)
{
  stack[ptr[0]++]=item;
  while(k&1)
    {
      stack[ptr[0]-2]+=stack[ptr[0]-1];
      ptr[0]--;
      k>>=1;
    }
  int_double res=stack[0];
  for(uint64_t p=1;p<ptr[0];p++)
    res+=stack[p];
  return res;
}


// rhs contains log n
// lhs contains log sigma(n)/n
int check_Robin(int_double lhs, int_double rhs, int verbosep)
{
  static bool init=false;
  static int_double exp_gam;
  if(!init)
    {
      init=true;
      exp_gam=exp(d_gamma);
    }
  int_double lln=log(rhs); // log log n
  int_double eglln=lln*exp_gam;
  int_double sn_n=exp(lhs);
  int_double del=sn_n-eglln;
  if(verbosep)
    {
      //printf("log log n = ");print_int_double(lln);
      //printf("\nexp(gamma) log log n = ");print_int_double(eglln);
      //printf("\nsigma(n)/n = ");print_int_double(sn_n);
      printf("lhs/rhs=");print_int_double(sn_n/eglln);printf("\n");fflush(stdout);}
  if(del.right<=0.0)
    {
      printf("Robin check failed.\n");
      printf("log log n = ");print_int_double(lln);
      printf("\nexp(gamma) log log n = ");print_int_double(eglln);
      printf("\nsigma(n)/n = ");print_int_double(sn_n);
      printf("\nlhs-rhs=");print_int_double(del);printf("\n");fflush(stdout);
      return false;
    }
  return true;
}

int main(int argc, char **argv)
{
  // set sse registers to round down
  // needed for int_double
  _fpu_rndd();
  // echo the command line
  printf("Command line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  // check command line parameters
  if(argc!=2)
    {
      printf("Usage:- %s <n its>\n",argv[0]);
      return 0;
    }
  uint64_t n_its=atol(argv[1]); // how many iterations


  uint64_t p1,p2,e;

  // setup primesieve
  primesieve_iterator pit;
  primesieve_init(&pit);

  int_double llhs,rrhs;
  uint64_t tp;
  for(uint64_t it=0;it<n_its;it++)
    {
      tp=primesieve_next_prime(&pit);
      int_double ip=int_double(tp);
      llhs=add_item(lhs,log1p(1.0/ip),it,&lhs_ptr);
      rrhs=add_item(rhs,log(ip),it,&rhs_ptr);
      if((it&0xffffff)==0) 
	{
	  printf("Doing prime %lu\n",tp);
	  check_Robin(llhs,rrhs,true);
	}
	  check_Robin(llhs,rrhs,false);
    }  
  printf("Last prime used = %lu\n",tp);
  check_Robin(llhs,rrhs,true);

  return 0;
}

