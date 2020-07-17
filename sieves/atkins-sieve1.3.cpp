// as yet this is bugged due to conversion
// of __int128_t to doubles which is a
// rounded conversion e.g.
// x=n;
// x=sqrt(x);
// may yield x^2 in [n-1,n+1] (I think)
// need to fix atkins_sieve(...) which
// assumes x^2 in [n-1,n]
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"
#include "../includes/pi_x.h"

#define primep(p,v) v[p/240]&mask1[p%240]
#define clear_prime(p,v) v[p/240]&=mask[p%240]
#define toggle_prime(p,v) v[p/240]^=mask1[p%240]

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}


inline ptype gcd (ptype a, ptype b)
// Euclid algorithm gcd
// best if a<=b ?
{
  unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
}


ptype mask1[240],mask[240];

void unset_primes(ptype *target, ptype len)
{
  ptype i;
  for(i=0;i<=len/240;i++)
    target[i]=0;
}
void set_primes(ptype *target, ptype len)
{
  ptype i;
  for(i=0;i<=len/240;i++)
    target[i]=0xFFFFFFFFFFFFFFFFL;
}

// return >=sqrt(x)
bigint sqrt_hi(bigint x)
{
  //printf("sqrt of ");print_bigint(x);
  bigint sqrt_x=sqrt((double) x);
  while(sqrt_x*sqrt_x>x)
    sqrt_x--;
  while(sqrt_x*sqrt_x<x)
    sqrt_x++;
  //printf("=");print_bigint(sqrt_x);printf("\n");
  return(sqrt_x);
}

void atkins1(long unsigned int delta, long unsigned int f, long unsigned int g, bigint L, long unsigned int B)
{
  long int x=f,x1=x;
  long unsigned int y0=g,y;
  bigint k0=(4*(bigint)f*(bigint)f+(bigint)g*(bigint)g-(bigint)delta)/(bigint) 60,k,k1=k0;
  bigint LB=L+B;
  /*
  while(k1<LB)
    {
      k1+=x1+x1+15;
      x1+=15;
      printf("x=%ld k0=",x1);print_bigint(k1);printf("\n");
    }
  printf("Starting with x=%ld k0=",x1);print_bigint(k1);printf("\n");
  */
  bigint m=ceil((double) (sqrt_hi(60*(LB-k0)+4*x*x)-x-x)/30.0);
  //printf("m=");print_bigint(m);printf("\n");
  k0+=m*(x+x)+15*m*m;
  x+=15*m;
  //printf("Starting with x=%ld k0=",x);print_bigint(k0);printf("\n");
  //exit(0);
  while(true)
    {
      if(x<15) return;
      x-=15;
      //if(x<0) return;
      k0-=x+x+15;
      while(k0<L)
	{
	  k0+=y0+15;
	  y0+=30;
	}
      k=k0;
      y=y0;
      while(k<LB)
	{
	  //printf("x=%ld,y=%lu,n=",x,y);print_bigint(k*60+delta);printf("\n");
	  k+=y+15;
	  y+=30;
	}
    }
  return;
}

int main(int argc, char **argv)
{
  long unsigned int f,g,delta;
  if(argc!=2)
    exit(0);
  int pow=atoi(argv[1])-5;
  bigint L=16667;
  for(int i=1;i<pow;i++) L*=10;
  long int B=(long unsigned int) 1<<33;
  for(f=1;f<=15;f++)
    for(g=1;g<=30;g++)
      {
	delta=(4*f*f+g*g)%60;
	if(gcd(delta,60)==1)
	  {
	    printf("Solutions to 4x^2+y^2 in[");
	    print_bigint(L*60);
	    printf(",");
	    print_bigint((L+B-1)*60);
	    printf("] with x=%ld (15) and y=%ld (30) and 4x^2+y^2=%ld (60)\n",f,g,delta);
	    atkins1(delta,f,g,L,B);
	    exit(0);
	  }
      }
  return(0);
}
