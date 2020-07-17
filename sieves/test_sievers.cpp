#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"
#include "../includes/pi_x.h"
#include "../primegen/primegen.h"
#include "../primegen/primegen.c"
#include "../primegen/primegen_next.c"
#include "../primegen/primegen_init.c"
#include "../primesieve/src/soe/PrimeSieve.h"


primegen pg[1];

void sieve_it(uint64_t p)
{
  uint64_t q=primegen_next(pg);
  if(p!=q)
    {
      printf("Mismatch at prime %lu vs %lu\n",p,q);
      exit(0);
    }
}

int main()
{
  PrimeSieve ps;
  ps.setSieveSize(64);
  primegen_init(pg);
  ps.generatePrimes(2,1000000000000L,sieve_it);  
  return(0);
}
