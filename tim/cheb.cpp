//
// File: cheb.cpp
//
// Compute x-theta(x) over a range
// assumes x(x0)-theta(x0) is zero
//
// Input: start x0 end y0
//        pi(x0) pi(y0)
//

#include "../includes/int_double12.0.h"
#include "primesieve.hpp"

#define MAX_N (1000000000000000000L) // 1e18
long unsigned int pi_x;
int_double dcheb;
int_double smallest_dcheb;

uint64_t last_prime;

void callback(uint64_t prime)
{
  pi_x++;
  int64_t diff=prime-last_prime;
  int_double del=log(int_double(prime))-diff;
  dcheb-=del;
  if(dcheb.left<smallest_dcheb)
    smallest_dcheb=dcheb;
  last_prime=prime;
}

int main(int argc,char **argv)
{
  // switch rounding mode for SSE
  // this makes int_doubles work
  _fpu_rndd();
  dcheb=0.0;
  smallest_dcheb=0.0;
  printf("Command line:-");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=6)
    {
      printf("Usage:- %s <start> <end> <pi(start)> <pi(end)> <cache size>.\n",argv[0]);
      exit(0);
    }
  long unsigned int x=atol(argv[1]);
  last_prime=x;
  long unsigned int y=atol(argv[2]);
  pi_x=atol(argv[3]);
  long unsigned int pi_y=atol(argv[4]);

  if((x>=y)||(pi_x>=pi_y))
    {
      printf("Expect start<end and pi(start)<pi(end).\n");
      exit(0);
    }
  if(y>MAX_N)
    {
      printf("maximum y supported=%lu. Exiting.\n",MAX_N);
      exit(0);
    }
  if((x&1)||(y&1))
    {
      printf("specify start and to be even to ensure they are composite.\n");
      exit(0);
    }

  primesieve::callback_primes(x,y,callback); // callback called for every prime
                                             // in [x,y]
  dcheb+=(y-last_prime); // add the bit on at the end
  if(pi_x!=pi_y)
    printf("Pi(%lu)should be %lu not %lu\n",y,pi_y,pi_x);
  print_int_double_str("x-theta(x) increased by ",dcheb);
  printf("between %lu and %lu\n",x,y);
  print_int_double_str("Largest decrease in x-theta(x)",smallest_dcheb);

  return(0);
}
