// Use a simple double precision interval implementation
// to compute theta(n)=sum p<=n log(p)

#include "int_double.h"
#include "primesieve.hpp"

// to count the primes
uint64_t prime_count=0;

// a stack to do pairwise summation without storing
// all inputs. STACK_SIZE needs to be >= log_2 number of inputs
uint64_t stack_ptr=0; // points to top of stack
#define STACK_SIZE (100)
int_double stack[STACK_SIZE];

// this is called by primesieve once for every p
// NB NOT NECESSARILY IN ORDER
void callback(uint64_t p)
{
  prime_count++; // count it
  if(stack_ptr==STACK_SIZE) // stack is about to overflow
    {
      printf("Stack overflow. Exiting.\n");
      exit(0);
    }
  stack[stack_ptr++]=log(int_double(p)); // push this log onto the stack
  uint64_t pc=prime_count;
  while((pc&1)==0)
    {
      stack[stack_ptr-2]+=stack[stack_ptr-1]; // add top 2 stack into 2nd
      stack_ptr--; // pop first
      pc>>=1;
    }
}

int main(int argc, char **argv)
{
  if(argc!=2)
    return 0;
  uint64_t maxn=atol(argv[1]);
  primesieve::callback_primes(2,maxn,callback);
  for(uint64_t i=1;i<stack_ptr;i++) // sum everything left on the stack
    stack[0]+=stack[i];

  // print some stuff
  printf("pi(%lu)=%lu\n",maxn,prime_count);
  printf("theta(%lu) in ",maxn);
  print_int_double(stack[0]);
  printf("\n");
  return 0;
}
