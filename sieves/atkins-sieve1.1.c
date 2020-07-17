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


// create base primes 2,3,5,...
inline void erat_sieve (ptype *erat, const ptype len)
{
  ptype i,j,sqrt_len;
  for(i=0,j=1;i<240;i++)
    if(gcd(i,30)==1)
      {
	mask1[i]=j;
	mask[i]=0xFFFFFFFFFFFFFFFFL;
	mask[i]^=j;
	j<<=1;
      }
  
  set_primes(erat,len);
  sqrt_len=sqrt(len);
  for(i=7;i<=sqrt_len;i+=2)
    if(primep(i,erat))
      for(j=i*i;j<=len;j+=i+i)
	clear_prime(j,erat);
}


// 
ptype floor_sqrt(bigint x)
{
  ptype y=sqrt(x);
  bigint y2=(bigint) y*y;
  while(y2<x)
    {
      y++;
      y2+=(y<<1)-1;
    }
  // y is now definately too big
  while(y2>x)
    {
      y--;
      y2-=(y<<1)+1;
    }
  return(y);
}

inline ptype ceil_sqrt(bigint x)
{
  ptype y=sqrt(x);
  bigint y2=(bigint) y*y;
  while(y2>x)
    {
      y--;
      y2-=(y<<1)+1;
    }
  // y is now definately too small
  while(y2<x)
    {
      y++;
      y2+=(y<<1)-1;
    }
  return(y);
}

void atkins_sieve(ptype *target, bigint A, ptype len)
{
  ptype x,y,xlim,count=0,ptr,rem;
  bigint fx2,y2,n,B=A+len-1;
  
  xlim=floor_sqrt((B-1)>>2);
  for(x=1,fx2=4;x<=xlim;x++,fx2+=(x<<3)-4)
    {
      if((x%1000000000)==0)
	printf("x=%ld\n",x);
      y=ceil_sqrt(A-fx2)|1; // y must be odd
      y2=(bigint)y*y;
      n=fx2+y2;
      while(n<=B)
	{
	  rem=n%12;
	  if((rem==1)||(rem==5))
	    {
	      ptr=n-A+1;
	      toggle_prime(ptr,target);
	    }
	  y+=2;
	  n+=(y<<2)-4;
	}
    }
}

int main()
{
  ptype *target;
  ptype i,j,target_len=(long unsigned int) 1<<35;
  bigint target_start;
  ptype count=0;

  for(i=0,target_start=1;i<24;i++,target_start*=10);
  while((target_start%30)!=1)
    target_start++;
  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len/240)))))
    fatal_error("Failed to allocate memory for target sieve.");

  unset_primes(target,target_len);

  atkins_sieve(target,target_start,target_len); // this leaves mp^2
  
  return(0);
}
