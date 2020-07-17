// See Tim Browning email of 27/10/17
#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "../includes/co_prime.h" // defines gcd and co_prime for uint64_t

inline int64_t safe_abs(const int64_t x)
{
  if(x<0)
    return -x;
  else
    return x;
}

int main(int argc, char **argv)
{
  printf("Command line:- %s",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");

  if(argc!=3)
    {
      printf("Usage:- %s <t> <B>\n",argv[0]);
      return 1;
    }

  int64_t t=atol(argv[1]);
  int64_t t3=3*t;
  int64_t tt=t*t;
  int64_t B=atol(argv[2]);

  for(int64_t x=1;x<=B;x++)
    for(int64_t y=1;y<=B;y++)
      {
	int64_t gxyt=gcd(x,gcd(y,t));
	int64_t xy=x*y;
	int64_t t3=3*t;
	int64_t x2=2*x;
	int64_t y5=5*y;
	int64_t n1=t3-x2-y5;
	int64_t n2=t3+x2-y5;
	int64_t n3=t3-x2+y5;
	int64_t n4=t3+x2+y5;
	int64_t n,z,zz;
	if(n1!=0)
	  {
	    n=n1*tt;
	    if(n%xy==0)
	      {
		z=n/xy;
		zz=safe_abs(z);
		if((zz<=B)&&co_prime(gxyt,zz))
		  printf("1: %ld %ld %ld %ld\n",x,y,z,t);
	      }
	  }
	if(n2!=0)
	  {
	    n=n2*tt;
	    if(n%xy==0)
	      {
		z=n/xy;
		zz=safe_abs(z);
		if((zz<=B)&&co_prime(gxyt,zz))
		  printf("2: %ld %ld %ld %ld\n",-x,y,-z,t);
	      }
	  }
	if(n3!=0)
	  {
	    n=n3*tt;
	    if(n%xy==0)
	      {
		z=n/xy;
		zz=safe_abs(z);
		if((zz<=B)&&co_prime(gxyt,zz))
		  printf("3: %ld %ld %ld %ld\n",x,-y,-z,t);
	      }
	  }
	if(n4!=0)
	  {
	    n=n4*tt;
	    if(n%xy==0)
	      {
		z=n/xy;
		zz=safe_abs(z);
		if((zz<=B)&&co_prime(gxyt,zz))
		  printf("4: %ld %ld %ld %ld\n",-x,-y,z,t);
	      }
	  }
      }
}
