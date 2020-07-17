// proth1.2.cpp
//
// look for a ladder of proth primes
// no more than STEP_SIZE apart.
//

#include "stdio.h"
#include "inttypes.h"
#include "stdlib.h"
#include "gmp.h"
/*
#include "../primegen/primegen.h"
#include "../primegen/primegen.c"
#include "../primegen/primegen_next.c"
#include "../primegen/primegen_init.c"
*/

#define STEP_SIZE (4000000000000000000L)
uint64_t SMALL_PRIME_LIMIT;

uint64_t *ismall_primes;

mpz_t log2_tmp,one,*small_primes,pp,pe,pq;

bool isprime (uint64_t p)
{
  if((p&1)==0)
    return(false);
  uint64_t j=3;
  while(j*j<p)
    if((p%j)==0)
      return(false);
    else
      j+=2;
  return(true);
}

void init_proth()
{
  mpz_init(log2_tmp);
  mpz_init(one);
  mpz_set_ui(one,1);
  uint64_t i;
  ismall_primes=(uint64_t *) malloc(sizeof(uint64_t)*SMALL_PRIME_LIMIT);
  small_primes=(mpz_t *)malloc(sizeof(mpz_t)*SMALL_PRIME_LIMIT);
  //primegen pg;
  //primegen_init(&pg);
  //uint64_t p=primegen_next(&pg);
  ismall_primes[0]=2;
  mpz_init(small_primes[0]);
  mpz_set_ui(small_primes[0],2);
  uint64_t p=3;
  for(i=1;i<SMALL_PRIME_LIMIT;i++)
    {
      while(!isprime(p))
	p+=2;
      ismall_primes[i]=p;
      mpz_init(small_primes[i]);
      mpz_set_ui(small_primes[i],ismall_primes[i]);
    }
  mpz_init(pp);mpz_init(pe);mpz_init(pq);
}

bool proth_p (uint64_t h, uint64_t n)
{
  uint64_t ptr;
  /* 
  for(ptr=1;ptr<SMALL_PRIME_LIMIT;ptr++) // skip 2
    if(mpz_fdiv_ui(pp,ismall_primes[ptr])==0)
      return(false);
  */

  __uint128_t pp128=h;
  pp128<<=n;
  pp128++;
  for(ptr=1;ptr<SMALL_PRIME_LIMIT;ptr++) // can't be even
    if((pp128%ismall_primes[ptr])==0)
      return(false);
  /*
// This is slower than trial divison by first 30 or so primes
  if(mpz_gcd_ui(NULL,pp,(uint64_t) 16294579238595022365L)!=0)
    return(false);
  */
  mpz_set_ui(pe,h);
  mpz_mul_2exp(pe,pe,n-1);
  mpz_add(pp,pe,pe);
  mpz_add_ui(pp,pp,1);
  for(ptr=0;ptr<SMALL_PRIME_LIMIT;ptr++)
    if(mpz_jacobi(small_primes[ptr],pp)==-1)
      break;
  if(ptr==SMALL_PRIME_LIMIT)
    return(false);
  mpz_powm(pq,small_primes[ptr],pe,pp);
  mpz_add_ui(pq,pq,1);
  return(mpz_cmp(pq,pp)==0);
}

int main(int argc, char **argv)
{

  if(argc!=3)
    exit(0);


  uint64_t n=atol(argv[1]);
  uint64_t h0=1LL<<(n-2);
  uint64_t h1=(h0<<2)-1;
  h0--;

  if(h1>=1000000000LL)
    h0=h1-1000000000LL;
  else
    h0=0;
  printf("Hardwired to do 1e9 steps h=[%lu,%lu].\n",h0,h1);
  SMALL_PRIME_LIMIT=atol(argv[2]);
  if(SMALL_PRIME_LIMIT<2)
    SMALL_PRIME_LIMIT=2;
  mpz_t h0z,h1z;
  mpz_init(h0z);mpz_init(h1z);
  mpz_set_ui(h0z,h0);
  mpz_set_ui(h1z,h1);

  init_proth();

  mpz_t two_n;
  mpz_init(two_n);
  mpz_mul_2exp(two_n,one,n);
  printf("n=%lu, 2^n=",n);mpz_out_str(NULL,10,two_n);printf("\n");

  printf("h0=");mpz_out_str(NULL,10,h0z);printf("\n");
  printf("h1=");mpz_out_str(NULL,10,h1z);printf("\n");

  mpz_t start,finish;
  mpz_init(start);mpz_init(finish);
  mpz_mul_2exp(start,h0z,n);
  mpz_add_ui(start,start,1);
  mpz_mul_2exp(finish,h1z,n);
  mpz_add_ui(finish,finish,1);

  printf("start=");mpz_out_str(NULL,10,start);printf("\n");
  printf("finish=");mpz_out_str(NULL,10,finish);printf("\n");


  printf("Searching for Proth primes...\n");
  uint64_t last_h=0,this_h;
  uint64_t first_h=STEP_SIZE>>(n+1);
  uint64_t del_h=STEP_SIZE>>n;
  for(this_h=first_h;this_h<=h1-h0;)
    {
      if(!proth_p(this_h+h0,n))
	this_h--;
      else
	{
	  //printf("Proth prime found at %lu*2^%lu+1.\n",this_h+h0,n);
	  last_h=this_h;
	  this_h=this_h+del_h;
	}
      if(this_h==last_h)
	{
	  printf("Failed to find Proth prime in region. Exiting.\n");
	  last_h+=del_h;
	  this_h=last_h+del_h;
	  //exit(0);
	}
    }

  this_h=h1-h0;
  while(last_h+first_h+h0<h1)
    {
      if(!proth_p(this_h+h0,n))
	this_h--;
      else
	{
	  printf("Final Proth prime found at %lu*2^%lu+1.\n",this_h+h0,n);
	  break;
	}
    }

  if(this_h+h0<h1-first_h)
    printf("Last h of %lu not close enough to end h %lu.\n",last_h,h1);

  return (0);
}
