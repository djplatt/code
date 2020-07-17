/*

compute int_X^2X |psi(x)-x|2 dx/X^2

*/
  
#include "stdio.h"
#include "stdlib.h"
#include "primesieve.h" // use Kim Walisch's primesieve
#include "math.h"

#define MAX_X (1000000000)

int main(int argc, char **argv)
{
  // echo the command line
  printf("Command line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  // check command line parameters
  primesieve_iterator it;
  primesieve_init(&it);
  uint64_t tp=primesieve_next_prime(&it); // 2

  double *psi,*psi2;
  psi=(double*)malloc(sizeof(double)*(MAX_X+1));
  psi2=(double*)malloc(sizeof(double)*(MAX_X+1));
  for(uint64_t i=0;i<=MAX_X;i++)
    psi[i]=0.0;
  psi2[0]=0.0;
  psi2[1]=0.0;

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
    {
      psi[i]=psi[i]+psi[i-1];
      psi2[i]=psi[i]*psi[i];
    }
  // psi[i] now contains psi(i)
  // psi2[i] contains psi^2(i)

  //for(uint64_t x=0;x<=MAX_X;x++)
  //printf("psi(%lu) = %e\n",x,psi[x]);


  for(int64_t i=0;i<=MAX_X;i++)
    psi[i]=psi[i]*(2.0*i+1.0);
  // psi[i] contains 2 int_{i}^{i+1} t*psi(t) dt


  //for(uint64_t x=0;x<=MAX_X;x++)
  //printf("Int %lu %lu 2*t*psi(t) = %e\n",x,x+1,psi[x]);
  //return 0;

  uint64_t X=1;
  double I1=psi[1]; // I1 now contains int_1^2 2*t*psi(t) dt 
  double I2=psi2[1]; // I2 contains int_1^2 psi^2(t) dt

  //for(uint64_t x=1;x<=MAX_X;x++)
  //printf("Int %lu %lu psi^2(t) = %e\n",x-1,x,psi2[x]);
  //return 0;

  //printf("%e %e ",I1,I2);
  printf("%lu %e\n",X,(I2-I1+7.0*X*X*X/3.0)/(X*X));
  //printf("int %lu %lu psi(t)^2 dt = %e\n",X,X+X,I2);
  //printf("int %lu %lu t*psi(t) dt = %e\n",X,X+X,I1);
  for(X=2;X<=MAX_X/2;X++)
    {
      I1-=psi[X-1];I1+=psi[X+X-2]+psi[X+X-1];
      I2-=psi2[X-1];I2+=psi2[X+X-2]+psi2[X+X-1];
      //printf("int %lu %lu psi(t)^2 dt = %e\n",X,X+X,I2);
      //printf("int %lu %lu 2*t*psi(t) dt = %e\n",X,X+X,I1);
      printf("%lu %e\n",X,(I2-I1+7.0*X*X*X/3.0)/(X*X));  
    }

  return 0;
}

