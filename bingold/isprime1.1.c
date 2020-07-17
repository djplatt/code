#include "stdlib.h"
#include "stdio.h"
#include "math.h"
typedef unsigned int u32;
typedef unsigned long u64;
typedef long i64;
typedef __uint128_t u128;

static u64 smallprime[]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};

// returns true if n is not divisible by any of the primes in the above
// list and passes a strong Fermat test to base 2
// assumes n is odd and 1 < n < 2^63
// From AR Booker
//
static inline int isprobableprime(u64 n) {
  //printf("ispp called with %lu\n",n);
	u64 q,one,neg1,x,y,nbar;
	int k;
	union { u128 d; u64 l[2]; } t;
#define mulredc(x,y) \
	(t.d=(u128)x*y,t.d+=(t.l[0]*nbar)*(u128)n,t.l[1]<n?t.l[1]:t.l[1]-n)
	/*
	for (k=0;k<sizeof(smallprime)/sizeof(smallprime[0]);k++)
		if (!(n % smallprime[k]))
			return (n == smallprime[k]);
	*/
	// compute nbar = -n^{-1} mod 2^64 by Newton's method
	k = -(int)n;
	k *= 2+k*(int)n; k *= 2+k*(int)n;
	k *= 2+k*(int)n; k *= 2+k*(int)n;
	nbar = (u64)k*(2+(u64)k*n);

	// compute Montgomery representations of 1 and -1
	one = 1+(~0ul)%n, neg1 = n-one;

	for (k=0,q=n>>1;!(q&1);k++) q>>=1;
	x = one+one; if (x >= n) x -= n;
	y = x, q>>=1;
	while (q) {
		x = mulredc(x,x);
		if (q & 1) y = mulredc(y,x);
		q>>=1;
	}
	if (y == one || y == neg1) return 1;
	while (--k >= 0) {
		y = mulredc(y,y);
		if (y == neg1) return 1;
	}
	return 0;
}

u64 masks[64],nmasks[64];

void init_masks()
{
  u64 i;
  for(i=0;i<64;i++)
    {
      masks[i]=1LL<<i;
      nmasks[i]=~masks[i];
    }
}
#define SIEVE_SIZE ((u64) 1L<<19)
#define set_bit(n,v) v[(n)>>6]|=masks[(n)&63] 
#define clear_bit(n,v) v[(n)>>6]&=(nmasks[(n)&63]);
#define is_bit(n,v) (v[(n)>>6]&masks[(n)&63])
#define isnt_bit(n,v) !is_bit(n,v)

u64 sieve[SIEVE_SIZE];

// n0 odd
void init_sieve(u64 n0)
{
  u64 i;
  for(i=0;i<SIEVE_SIZE;i++)
    sieve[i]=~(0L);
  for (i=0;i<sizeof(smallprime)/sizeof(smallprime[0]);i++)
    {
      u64 j=n0%smallprime[i];
      if(j)
	{
	  j=(smallprime[i]-j);
	  if(j&1)
	    j=(j+smallprime[i])>>1;
	  else
	    j>>=1;
	}
      while(j<(SIEVE_SIZE<<6))
	{
	  clear_bit(j,sieve);
	  j+=smallprime[i];
	}
    }
  for(i=0;i<SIEVE_SIZE<<6;i++)
    if(is_bit(i,sieve))
      if(!isprobableprime(n0+i*2))
	clear_bit(i,sieve);
}

#define DEL (1L<<20)
u64 prime_list[DEL/4];
u64 crossed[DEL/2];

// show that all even numbers in [N,N+DEL] can be written as p+q
// sieve starts at N0
void find_pairs(u64 N, u64 N0)
{
  printf("Looking for pairs from %lu to %lu with sieve from %lu to %lu\n",
	 N,N+DEL-2,N0,N0+(SIEVE_SIZE<<7)-1);
  //u64 prime_list[DEL/4];
  //u64 crossed[DEL/2];
  u64 N2=N>>1,sieve_start=(N2-(N0-1))>>1,sieve_end=sieve_start+DEL/4;
  u64 p_ptr=0,i,j;
  for(i=0;i<DEL/2;i++)
    crossed[i]=1;
  //printf("using sieve from %lu to %lu\n",sieve_start,sieve_end-1);
  u64 count=0;
  for(i=sieve_start;i<sieve_end;i++)
    if(is_bit(i,sieve))
      {
	prime_list[p_ptr]=N0+i*2;
	//printf("Found prime %lu\n",prime_list[p_ptr]);
	for(j=0;j<=p_ptr;j++)
	  {
	    u64 cp=(prime_list[p_ptr]+prime_list[j]-N)>>1;
	    if(crossed[cp])
	      {
		crossed[cp]=0;
		count++;
	      }
	  }
	p_ptr++;
      }
  printf("We found %lu primes and crossed out %lu pairs.\n",p_ptr,count);
  count=0;
  for(i=0;i<DEL/2;i++)
    if(crossed[i]) count++;
  printf("Couldn't do %lu on first pass.\n",count);
}

#define TARGET_START (1L<<61)
#define TARGET_END (TARGET_START+(1L<<30))
    
int main(int argc, char **argv)
{
  u64 N0;
  init_masks();
  for(N0=TARGET_START/2;N0<TARGET_END/2;N0+=SIEVE_SIZE<<7)
    {
      init_sieve(N0+1);
      u64 i,count=0;
      for(i=0;i<SIEVE_SIZE;i++)
	count+=__builtin_popcountl(sieve[i]);
      printf("There are %lu primes between %lu and %lu\n",count,N0,N0+(SIEVE_SIZE<<7)-1);
      u64 N;
      for(N=N0<<1;N<(N0<<1)+(SIEVE_SIZE<<7);N+=DEL)
	find_pairs(N,N0+1);
    }
  return(0);
}
