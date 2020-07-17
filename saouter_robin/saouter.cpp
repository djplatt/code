#include "arb.h"
#include <primesieve.hpp>
#include "inttypes.h"

#define STACK_SIZE (1024)

arb_t sum,tmp1,tmp2,stack[STACK_SIZE];
int64_t prec,stack_ptr,stack_count;


void init_stack()
{
  for(uint64_t i=0;i<STACK_SIZE;i++)
    arb_init(stack[i]);
  stack_ptr=0;
  stack_count=0;
}

void add_to_stack(arb_t x,int64_t prec)
{
  if(stack_ptr==STACK_SIZE)
    {
      printf("ARB stack overflow. Exiting.\n");
      exit(0);
    }
  arb_set(stack[stack_ptr++],x);
  stack_count++;
  uint64_t i=stack_count;
  while((i&1)==0)
    {
      arb_add(stack[stack_ptr-2],stack[stack_ptr-2],stack[stack_ptr-1],prec);
      stack_ptr--;
      i>>=1;
    }
}

void return_stack(arb_t res, int64_t prec)
{
  arb_set(res,stack[0]);
  for(uint64_t i=1;i<stack_ptr;i++)
    arb_add(res,res,stack[i],prec);
}

void clear_stack()
{
  for(uint64_t i=0;i<STACK_SIZE;i++)
    arb_clear(stack[i]);
}

void callback(uint64_t p)
{
  arb_set_ui(tmp1,p);
  arb_inv(tmp2,tmp1,prec);
  arb_neg(tmp2,tmp2);
  arb_log1p(tmp1,tmp2,prec);
  add_to_stack(tmp1,prec);
}


int main(int argc, char **argv)
{
  printf("Command line:- %s ",argv[0]);
  for(uint64_t a=1;a<argc;a++)
    printf("%s ",argv[a]);
  printf("\n");
  if(argc!=4)
    {
      printf("Usage:- %s <start> <end> <prec>.\n",argv[0]);
      return 0;
    }

  arb_init(sum);
  arb_init(tmp1);arb_init(tmp2);
  prec=atol(argv[3]);
  init_stack();
  primesieve::callback_primes(atol(argv[1]),atol(argv[2]),callback);
  return_stack(sum,prec);
  clear_stack();
  printf("Sum returned ");arb_printd(sum,40);printf("\n");
  arb_clear(sum);
  arb_clear(tmp1);
  arb_clear(tmp2);
  return 0;
}
