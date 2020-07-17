#include "pari.h"
#include "math.h"

#define bitp ((unsigned long int)(1<<22)-1)

void print_usage()
{
  exit(0);
}


int main()
{
  long int n,count=0;
  double res=0.0,high_point=0.0;
  GEN pn;

  pari_init(1<<20,1<<20);
  pn=stoi(1);
  for(n=0;n<24;n++)
    {
      printf("multiplying %ld\n",n);
      pn=gmulgs(pn,10);
    }
  pn=gaddgs(pn,21); // pn=10^24+21=1 mod 31
  for(n=0;n<10000;n++)
    {
      if(isprime(pn))
	count++;
      pn=gaddgs(pn,30);
    }
  printf("%ld primes found\n",count);
  return(0);
}
