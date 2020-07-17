#include "stdlib.h"
#include "stdio.h"
#include "math.h"
typedef unsigned int u32;
typedef unsigned long u64;
typedef long i64;
typedef __uint128_t u128;

static u64 smallprime[]={3,5,7,11,13,17,19,23,29,31,37,41,43};


// returns true if n is not divisible by any of the primes in the above
// list and passes a strong Fermat test to base 2
// assumes n is odd and 1 < n < 2^63
// From AR Booker
// I do single GCD rather than many truial divisions
static inline int isprobableprime(u64 n) {
	u64 q,one,neg1,x,y,nbar;
	int k;
	union { u128 d; u64 l[2]; } t;
#define mulredc(x,y) \
	(t.d=(u128)x*y,t.d+=(t.l[0]*nbar)*(u128)n,t.l[1]<n?t.l[1]:t.l[1]-n)
		
	for (k=0;k<sizeof(smallprime)/sizeof(smallprime[0]);k++)
		if (!(n % smallprime[k]))
			return (n == smallprime[k]);
	/*
	if (!(n % 3)) return (n == 3);
	if (!(n % 5)) return (n == 5);
	if (!(n % 7)) return (n == 7);
	if (!(n % 11)) return (n == 11);
	if (!(n % 13)) return (n == 13);
	if (!(n % 17)) return (n == 17);
	if (!(n % 19)) return (n == 19);
	if (!(n % 23)) return (n == 23);
	if (!(n % 29)) return (n == 29);
	if (!(n % 31)) return (n == 31);
	if (!(n % 37)) return (n == 37);
	if (!(n % 41)) return (n == 41);
	if (!(n % 43)) return (n == 43);
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




/*
typedef struct{
  u64 prime;
  u64 word_offset;
  u64 bit_offset;
} prime_t;



#define set_bit(n,v) v[(n)>>6]|=masks[(n)&63] 
#define clear_bit(n,v) v[(n)>>6]&=(nmasks[(n)&63]);
#define is_bit(n,v) (v[(n)>>6]&masks[(n)&63])
#define isnt_bit(n,v) !is_bit(n,v)

#define SIEVE_LENGTH (1024)
void binary_sieve()
{
  u64 masks[64];
  u64 nmasks[64];

  u64 sieve[SIEVE_LENGTH];
  u64 prime_masks[64][SIEVE_LENGTH+1];
  u64 i,j,k,lim=sqrt(SIEVE_LENGTH*64*2+1);
  for(i=0;i<64;i++)
    {
      masks[i]=1LL<<i;
      nmasks[i]=~masks[i];
      //printf("mask[%lu]=%lx %lx\n",i,masks[i],nmasks[i]);
    }
  for(i=0;i<SIEVE_LENGTH;i++)
    {
      sieve[i]=0;
      for(j=0;j<64;j++)
	prime_masks[j][i]=~0;
    }
  for(k=1;k<64;k++)
    for(i=0;i<k;i++)
      clear_bit(i,prime_masks[k]);
  //printf("limit=%lu\n",lim);
  for(i=3;i<=lim;i+=2)
    {
      //printf("i=%lu\n",i);
      if(is_bit((i-1)>>1,prime_masks[0]))
	for(j=i*i;j<SIEVE_LENGTH*64*2+1;j+=(i<<1))
	  {
	    //printf("clearing %lu bit ptr %lu mask %lx\n",j,(j-3)>>1,masks[((j-3)>>1)&63]);
	    for(k=0;k<64;k++)
	      {
		//printf("before clear %lu we have %lx\n",j,prime_masks[k][((j-3)>>1)>>6]);
		clear_bit(((j-3)>>1)+k,prime_masks[k]);
	      }
	  }
    }

  u64 pcount=0;
  for(i=0;i<SIEVE_LENGTH;i++)
    pcount+=__builtin_popcountl(prime_masks[0][i]);
  printf("Found %lu primes.\n",pcount);

  prime_t *primes=(prime_t *)malloc(sizeof(prime_t)*pcount);
  if(!primes)
    {
      printf("Error allocating memory for primes. Exiting.\n");
      exit(0);
    }

  u64 p_ptr=0;
  for(i=0;i<SIEVE_LENGTH;i++)
    for(j=0;j<64;j++)
      if(prime_masks[0][i]&masks[j])
	{
	  primes[p_ptr].word_offset=i;
	  primes[p_ptr].bit_offset=j;
	  primes[p_ptr].prime=(((i<<6)+j)<<1)+3;
	  p_ptr++;
	}

  u64 sieve_ptr=0,b=0;

  
  u64 n=1000000; // start of sieve
  while(n<2000000)
    {
      if(sieve[sieve_ptr]==(~0)) // all these have been done
	{
	  sieve[sieve_ptr]=0;
	  sieve_ptr++;
	  n+=128;
	  b=0;
	  if(sieve_ptr==SIEVE_LENGTH)
	    sieve_ptr=0;
	  continue;
	}
      while(sieve[sieve_ptr]&masks[b]) b++;
      // n+2b hasn't been crossed off
      // find a p such that n+2b-p is prime
      for(i=0;i<pcount;i++)
	if(isprobableprime(n+(b<<1)-primes[i].prime))
	  break;
      printf("Using prime %lu for n=%lu\n",primes[i].prime,n+(b<<1));
      printf("sieve[0]=%x\n",sieve[0]);
      u64 mask_no=primes[i].bit_offset;
      u64 pp=primes[i].word_offset;

      for(i=sieve_ptr;(i<SIEVE_LENGTH)&&(pp<SIEVE_LENGTH);i++,pp++)
	  sieve[i]|=prime_masks[mask_no][pp];
      if(pp<SIEVE_LENGTH)
	for(i=0;pp<SIEVE_LENGTH;i++,pp++)
	  sieve[i]|=prime_masks[mask_no][pp];
      printf("sieve[0]=%x\n",sieve[0]);
      b++;
      
    }
  
}
*/

#define MAX_PRIME (9781) // largest prime needed according to Oliviera e Silva
#define SQRT_MP (97) // no need to go further than this
#define N_PRIMES (1205) // excludes 2
#define SIEVE_LENGTH ((u64) 1<<20)
typedef unsigned char u8;
void binary_sieve(u64 n0, u64 no_its)
{
  u8 sieve[SIEVE_LENGTH];
  /*
  u8 *sieve=(u8 *)malloc(sizeof(u8)*SIEVE_LENGTH);
  if(!sieve)
    {
      printf("Fatal error allocating memory for sieve. Exiting.\n");
      exit(0);
    }
  */
  u8 primes[MAX_PRIME+1];
  u64 lprimes[N_PRIMES];
  u8 gaps[N_PRIMES-1];
  u64 i,j,k;
  for(i=3;i<=MAX_PRIME;i+=2)
    primes[i]=1;
  for(i=3;i<=SQRT_MP;i+=2)
    if(primes[i])
      for(j=i*i;j<=MAX_PRIME;j+=(i<<1))
	primes[j]=0;
  for(i=0,j=3;i<N_PRIMES;j+=2)
    if(primes[j])
      lprimes[i++]=j;
  for(i=0;i<N_PRIMES-1;i++)
    gaps[i]=(lprimes[i+1]-lprimes[i])>>1;

  u64 it;
  for(it=0;it<no_its;it++,n0+=(SIEVE_LENGTH<<1))
    {
      for(i=0;i<SIEVE_LENGTH;i++)
	sieve[i]=1;
      
      for(i=0;i<SIEVE_LENGTH;i++)
	if(sieve[i])
	  for(j=0;;j++)
	    if(isprobableprime((i<<1)+n0-lprimes[j]))
	      {
		//printf("Base prime %lu\n",(i<<1)+n0-lprimes[j]);
		for(k=i;(j<N_PRIMES-1)&&(k<SIEVE_LENGTH);k+=gaps[j++])
		  {
		    //printf("small prime %lu - killing %lu\n",lprimes[j],n0+(k<<1));
		    sieve[k]=0;
		  }
		break;
	      }
    }
}
    
#define ITS ((u64) (1LL<<32)/SIEVE_LENGTH)
int main(int argc, char **argv)
{
  u64 n0,i;
  n0=1LL<<62;
  printf("Starting at %lu\n",n0);
  binary_sieve(n0,ITS);
  printf("Finished at %lu\n",n0+SIEVE_LENGTH*ITS*2-2);
  return(0);
}
