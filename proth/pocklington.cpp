// pocklington.cpp
//
#include "stdio.h"
#include "inttypes.h"
#include "stdlib.h"
#include "gmp.h"
#include "../primegen/primegen.h"
#include "../primegen/primegen.c"
#include "../primegen/primegen_next.c"
#include "../primegen/primegen_init.c"


#define STEP_SIZE (4000000000000000000L)


int main(int argc, char **argv)
{

  if(argc!=5)
    exit(0);


  uint64_t n=atol(argv[1]);
  uint64_t h0=atol(argv[2]);
  uint64_t h1=atol(argv[3]);
  PRIME_LIMIT=atol(argv[4]);

  mpz_t h0z,h1z;
  mpz_init(h0z);mpz_init(h1z);
  mpz_set_ui(h0z,h0);
  mpz_set_ui(h1z,h1);

  init_proth();

  mpz_t two_n;
  mpz_init(two_n);
  mpz_mul_2exp(two_n,one,n);
  if(mpz_cmp(h1z,two_n)>=0)
    {
      printf("h1 exceeds 2^n. Exiting.\n");
      exit(0);
    }
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


  uint64_t bit_array_len=((h1-h0)>>6)+1;

  uint64_t *bits=(uint64_t *)malloc(bit_array_len<<3);
  if(!bits)
    {
      printf("Failed to allocate memory for bit array. Exiting.\n");
      exit(0);
    }

  uint64_t i;
  for(i=0;i<bit_array_len;i++)
    bits[i]=0xFFFFFFFFFFFFFFFFL;

  primegen pg;
  uint64_t p;
  primegen_init(&pg);
  p=primegen_next(&pg); // skip 2

  printf("enumerating primes to %ld\n",PRIME_LIMIT);
  for(p=primegen_next(&pg);p<PRIME_LIMIT;p=primegen_next(&pg))
    {
      // this could be done in integers
      uint64_t r1=mpz_fdiv_ui(start,p);
      uint64_t r2=mpz_fdiv_ui(two_n,p);
      uint64_t k=solve_mod(r1,r2,p);
      uint64_t ptr=k;
      while(ptr<=(h1-h0))
	{
	  //printf("Clearing ptr %lu\n",ptr);
	  clear_bit(ptr,bits);
	  ptr+=p;
	}
    }

  printf("Searching for Proth primes...\n");
  uint64_t ptr=0;
  uint64_t last_h=0,this_h;
  uint64_t first_h=STEP_SIZE>>(n+1);
  uint64_t del_h=STEP_SIZE>>n;
  for(this_h=first_h;this_h<=h1-h0;)
    {
      while((!test_bit(this_h,bits))&&(this_h>last_h))
	this_h--;
      if(this_h==last_h)
	{
	  printf("Failed to find proth prime after h=%lu. Exiting.\n",last_h);
	  exit(0);
	}
      //printf("Proth testing for h=%lu\n",this_h+h0);
      if(!proth_p(this_h+h0,n))
	this_h--;
      else
	{
	  printf("Proth prime found at %lu*2^%lu+1.\n",this_h+h0,n);
	  last_h=this_h;
	  this_h=this_h+del_h;
	}
    }

  this_h=h1-h0;
  while(last_h+first_h+h0<h1)
    {
      while((!test_bit(this_h,bits))&&(this_h>last_h))
	this_h--;
      if(this_h==last_h)
	{
	  printf("Failed to find Final proth prime after h=%lu. Exiting.\n",last_h+h0);
	  exit(0);
	}
      //printf("Proth testing for h=%lu\n",this_h+h0);
      if(!proth_p(this_h+h0,n))
	this_h--;
      else
	{
	  printf("Final Proth prime found at %lu*2^%lu+1.\n",this_h+h0,n);
	  break;
	}
    }

  if(this_h+h0<h1-first_h)
    printf("Last h of %lu not close enough to end h %lu.\n",last_h+h0,h1);

  return (0);
}
