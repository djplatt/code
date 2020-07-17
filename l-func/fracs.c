#include "stdio.h"
#include "stdlib.h"

#define Q (10000)
#define SQRT_Q (100)
#define SQRT_SQRT_Q (10)
#define N_PRIMES (25)
#define TRUE (1)
#define FALSE (0)
#define SUCCESS (1)
#define FAILURE (0)

unsigned char primes[SQRT_Q];
unsigned int prime_n[N_PRIMES];

void build_primes()
{
  unsigned int ptr,prime;

  for(ptr=1;ptr<SQRT_Q;ptr++)
    primes[ptr]=TRUE;
  primes[0]=FALSE;

  for(prime=2;prime<SQRT_SQRT_Q;prime++)
    if(primes[prime-1])
      for(ptr=prime*prime-1;ptr<SQRT_Q;ptr=ptr+prime)
	primes[ptr]=FALSE;
  ptr=0;
  for(prime=1;prime<SQRT_Q;prime++)
    if(primes[prime-1])
	prime_n[ptr++]=prime;
};

typedef struct {double val; unsigned int num; unsigned int den;} frac;


int gcd (unsigned int a, unsigned int b)
{
  unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

int co_prime(unsigned int a, unsigned int b)
{
  return(gcd(a,b)==1);
};

frac fracs[Q*Q];

inline int _frac_comp (frac a, frac b)
{
  if(a.val-b.val>0)
    return(1);
  return(-1);
};

int frac_comp(const void *a,const void *b)
{
  return(_frac_comp(*(frac*)a,*(frac*)b));
};  

int main(int argc, char **argv)
{
  unsigned int num,den,num1,den1,ptr,ptr2,pc;
  double max_diff;

  build_primes();
  /*
  for(ptr=0;ptr<N_PRIMES;ptr++)
      printf("%d is prime\n",prime_n[ptr]);


  for(ptr=1;ptr<11;ptr++)
    for(ptr2=1;ptr2<11;ptr2++)
      if(co_prime(ptr,ptr2))
	printf("%d and %d are co-prime.\n",ptr,ptr2);
  */

  //  fracs[0]=1.0;
  ptr=0;
  for(den=1;den<=Q;den++)
    for(num=1;num<=(den>>1);num++)
      if(co_prime(den,num))
	{
	  fracs[ptr].val=((double) num)/((double) den);
	  fracs[ptr].num=num;
	  fracs[ptr++].den=den;
	};
  printf("got %d distinct fractions.\n",ptr);
  /*
  for(num=0;num<ptr;num++)
    printf("%f\n",fracs[num]);
  */

  qsort(fracs,ptr,sizeof(frac),frac_comp);

  
  for(ptr2=1;ptr2<ptr;ptr2++)
    fracs[ptr2-1].val=fracs[ptr2].val-fracs[ptr2-1].val;
  
  qsort(fracs,ptr-1,sizeof(frac),frac_comp);

  printf("%12e %12\n",fracs[0].val,fracs[ptr-2].val);
  /*
  ptr2=0;
  for(num=1;num<50;num++)
    {
      den=0;
      for(;fracs[ptr2++].val<((double) num)/10000000.0;den++);
      printf(" %9.7f %9.7f %10d\n",(num-1)/10000000.0,num/10000000.0,den);
    };
  printf(">%9.7f %20d\n",50.0/10000000.0,ptr-ptr2);
  */

  for(pc=9900;pc<10000;pc++)
    printf("%6.2f %12e\n",pc/100.0,fracs[(int) (ptr/10000.0*pc)].val);
  printf("100.00 %12e\n",fracs[ptr-2].val);
  return(SUCCESS);
};

