#include "stdio.h"
#include "inttypes.h"

#define PRIMORIAL ((long unsigned int) 2*3*5*7*11*13)
#define NEXT_PRIME ((long unsigned int) 17)


inline long unsigned int gcd (long unsigned int a, long unsigned int b)
/* Euclid algorithm gcd */
{
	long unsigned int c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

inline int co_prime(long unsigned int a, long unsigned int b)
{
	return(gcd(a,b)==1);
};



long unsigned int cross_count=0;


void cross_out(long unsigned int p)
{
  cross_count++;
}

int main()
{
  long unsigned int i,del=1,p,max;
  max=((long unsigned int) 100000000000);

  for(p=PRIMORIAL+1;p<max;p+=PRIMORIAL<<1)
    cross_out(p);

  for(p=NEXT_PRIME;p<max;p+=PRIMORIAL<<1)
    cross_out(p);

  for(del=NEXT_PRIME+2;del<PRIMORIAL;del+=2)
    if(co_prime(del,PRIMORIAL))
      for(p=del;p<max;p+=PRIMORIAL<<1)
	cross_out(p);

  printf("There were %lu crossings out\n",cross_count);
  return(0);
}
