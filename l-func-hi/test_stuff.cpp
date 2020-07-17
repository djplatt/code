#include <iostream>
#include <cstdlib>
#include "characters.h"

using namespace std;

#define PREC (300)

void foo(uint64_t q)
{
  DirichletGroup G(q);

  printf("phi(%lu)=%lu\n",G.q,G.phi_q);
  printf("%lu\n",G.primes->at(0));

}

int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:= %s <q>.\n",argv[0]);
      exit(0);
    }

  dft_init(PREC);

  uint64_t q=atol(argv[1]);

  foo(q);

  DirichletGroup G(q);
  uint64_t val_at=0;
  for(uint64_t i=2;i<q;i++)
    if(G.is_coprime_to_q(i))
      {
	val_at=i;
	break;
      }

  for(uint64_t i=1;i<q;i++)
    {
      if(!G.is_coprime_to_q(i)) continue;
      DirichletCharacter c=G.character(i);
      //if(!c.is_primitive()) continue;
      printf("%lu ",i);
      if(c.is_even())
	printf("even ");
      else
	printf("odd ");
      
      if(c.is_primitive())
	printf("primitive     ");
      else
	printf("not primitive ");
      
      cout << c.value(val_at);
      printf("\n");
    }
  
  return(0);
}
