/*

Compute R_t(N_n) = prod p<=p_n (1-p^-t)/(1-p^-1) / log log N_n
p_n is nth prime
N_n is p_n# product of first n primes.

Identify lowest n such that R_t(N_n)< delta * exp(gamma)

Uses int_double package

D.J.Platt 2018
*/
  
#include "stdio.h"
#include "stdlib.h"
#include "primesieve.h" // use Kim Walisch's primesieve
#undef ulong
#include "../includes/int_double14.0.h"


#define STACK_LEN (128)

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
int_double mul_item(int_double *stack, int_double item, uint64_t k, uint64_t *ptr)
{
  stack[ptr[0]++]=item;
  while(k&1)
    {
      stack[ptr[0]-2]*=stack[ptr[0]-1];
      ptr[0]--;
      k>>=1;
    }
  int_double res=stack[0];
  for(uint64_t p=1;p<ptr[0];p++)
    res*=stack[p];
  return res;
}

// log(1-p^-t)-log(1-p^-1)
int_double log_p_bit(uint64_t p, int_double log_p, int64_t t)
{
  int_double ln=log1p(-exp(-t*log_p));
  int_double ld=log1p(-1.0/int_double(p));
  return ln-ld;
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
  int_double delta;
  // check command line parameters
  if(argc==2)
    delta=0.0;
  else
    {
      if(argc==4)
	{
	  delta=atol(argv[2]);
	  delta/=atol(argv[3]);
	  delta=log(delta);
	  print_int_double_str("log(delta) = ",delta);
	}
      else
	{
	  printf("Usage:- %s <t> [<delta num> <delta den>]\n",argv[0]);
	  return 0;
	}
    }
  int64_t t=atol(argv[1]);

  // setup primesieve
  primesieve_iterator it;
  primesieve_init(&it);

  int_double num_stack[STACK_LEN];
  uint64_t num_ptr=0;
  int_double den_stack[STACK_LEN];
  uint64_t den_ptr=0;

  primesieve_next_prime(&it); // p=2
  int_double log2=log(int_double(2));
  add_item(num_stack,log_p_bit(2,log2,t),0,&num_ptr); // log(1-2^-t)-log(1-2^-1)
  add_item(den_stack,log2,0,&den_ptr); // log 2

  for(uint64_t k=1;;k++)
    {
      uint64_t p=primesieve_next_prime(&it); // next prime starting from 3
      int_double lp=log(int_double(p));
      int_double lnum=add_item(num_stack,log_p_bit(p,lp,t),k,&num_ptr); // log prod 
      int_double den=log(add_item(den_stack,lp,k,&den_ptr)); // log log p#
      int_double res=exp(lnum-d_gamma+delta)/den; // prod/log log p#/exp(gamma)
      if((k%10000000)==0)
	{
	  printf("p=%lu ",p);
	  print_int_double_str("res=",res);
	  print_int_double_str("log(log(p#))=",den);
	  fflush(stdout);
	}
      //printf("%lu ",p);print_int_double(lnum);printf(" ");print_int_double(eden);printf(" ");
      //print_int_double_str("res=",res);
      if(res.right>-1.0) // res < 1.0
	{
	  printf("Success with p=%lu\n",p);
	  print_int_double_str("res=",res);
	  print_int_double_str("log(log(p#))=",den);
	  break;
	}
    }
  return 0;
}

