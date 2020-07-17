/*

compute int_X^2X |psi(x)-x|2 dx/X^2

*/
  
#include "stdio.h"
#include "stdlib.h"
#include "primesieve.h" // use Kim Walisch's primesieve
#include "math.h"

#define MAX_X (10)

double term(uint64_t a, double psa)
{
  double res=psa*(psa-2*a-1);
  //printf("term(%lu,%e) returning %e\n",a,psa,res);
  return res;
}

int main(int argc, char **argv)
{
  primesieve_iterator it;
  primesieve_init(&it);
  uint64_t tp=primesieve_next_prime(&it); // 2

  double *psi;
  psi=(double*)malloc(sizeof(double)*(MAX_X+1));
  for(uint64_t i=0;i<=MAX_X;i++)
    psi[i]=0.0;
  while(tp<=MAX_X)
    {
      double lp=log(tp);
      uint64_t pn=tp;
      while(pn<=MAX_X)
	{
	  psi[pn]=lp;
	  pn*=tp;
	}
      tp=primesieve_next_prime(&it);
    }
  // psi[i] contains lambda(i)

  for(uint64_t i=2;i<=MAX_X;i++)
    psi[i]=psi[i]+psi[i-1];
  // psi[i] now contains psi(i)

  //for(uint64_t x=0;x<=MAX_X;x++)
  //printf("psi(%lu) = %e\n",x,psi[x]);

  uint64_t X=1;
  double I1=term(1,psi[1]); // I1 now contains int_1^2 psi(t)^2 -2 t psi(t) dt 

  printf("%lu %e\n",X,7.0*X/3.0+I1/(X*X));
  for(X=2;X<=MAX_X/2;X++)
    {
      I1-=term(X-1,psi[X-1]);
      I1+=term(X+X-2,psi[X+X-2]);
      I1+=term(X+X-1,psi[X+X-1]);
      printf("%lu %e\n",X,7.0*X/3.0+I1/(X*X));  
    }

  return 0;
}

