//
// File: bestcase.cpp
//
// Code to do Tim Trudgian's counting
//
// DJ Platt 
// Setup
//
// we are given a positive odd integer d>1 and 
// we have K variables, v_1 to v_k.
// We let the v_i range from 1 to ord_d(2) and count
// how many ways we can get each of the d possible values for N 
// N=2^{v_1}+2^{v_2}+...+2^{v_K} mod d
// We return the minimum found.
// 
// Assume the number of ways distribute perfectly (the best case)
// so each squarefree odd d makes a contribution of
//
// 1/(d*k(d))
//
#include "stdlib.h"
#include "stdio.h"
#include "malloc.h"
#include "inttypes.h"
#include "math.h"
//#include "fmpz.h"
#include "gmp.h"
#include "mpfi.h"
#include "mpfi_io.h"

#define MAX_D ((uint64_t) 100000000) 
#define MAX_FACS ((uint64_t) 7) // 2*3*5*7*11*13*17>MAX_D

mpfi_t one;

// structure to hold factor information
typedef struct
{
	uint64_t pr;   // = 0 if no primitive root
	uint64_t phi;  // Euler phi(q)
	uint64_t num_facs;  // number of prime power factors
	uint64_t facs[MAX_FACS];  // the factors p^n
	uint64_t primes[MAX_FACS]; // the prime factors p
} factor;

inline uint64_t gcd (uint64_t a, uint64_t b)
/* Euclid algorithm gcd */
{
	uint64_t c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

inline int co_prime(uint64_t a, uint64_t b)
{
	return(gcd(a,b)==1);
};

inline int phi(int p, int p_n)
{
  return(p_n-p_n/p);
}

uint64_t pow_mod(uint64_t a, uint64_t b, uint64_t m)
{
  uint64_t a_pw=a,pw=b,res=1;
  while(true)
    {
      //printf("pw=%ld a_pw=%ld res=%ld\n",pw,a_pw,res);
      if(pw&1)
	res=(res*a_pw)%m;
      pw>>=1;
      if(pw==0)
	return(res);
      a_pw=(a_pw*a_pw)%m;
    }
}
      
unsigned long int pr(unsigned long int i, factor *factors)
{
  int phi=factors[i].phi;
  for(int p=2;p<i;p++)
    {
      if(gcd(p,i)!=1)
	continue;
      bool good=true;
      for(int j=0;j<factors[phi].num_facs;j++)
	{
	  if(pow_mod(p,phi/factors[phi].primes[j],i)!=1)
	    continue;
	  good=false;
	  break;
	}
      if(good)
	return(p);
    }
}

bool make_factors(factor *factors, int q_end)
{
  printf("Making factor database.\n");
  int *primes,max_p=floor(sqrt(q_end));
  primes=(int *) malloc(sizeof(int)*(q_end+1));
  for(int i=2;i<=q_end;i++)
    primes[i]=0;
  for(int i=4;i<=q_end;i+=2)
    primes[i]=2;
  for(int i=3;i<=max_p;i+=2)
    if(primes[i]==0)
      for(int j=i*i;j<=q_end;j+=i)
	if(primes[j]==0)
	  primes[j]=i;
  printf("Prime sieve completed.\n");
  // now each entry primes[i] is 0 if i prime, else = smallest prime factor
  factors[3].num_facs=1;
  factors[3].primes[0]=3;
  factors[3].facs[0]=3;
  factors[4].num_facs=1;
  factors[4].primes[0]=2;
  factors[4].facs[0]=4;

  for(int f=5;f<=q_end;f++)
    if(primes[f]==0) // a prime
      {
	factors[f].num_facs=1;
	factors[f].primes[0]=f;
	factors[f].facs[0]=f;
      }
    else
      {
	factors[f].primes[0]=primes[f];
	int n=f/primes[f];
	if(factors[n].primes[0]==primes[f])
	  {
	    factors[f].num_facs=factors[n].num_facs;
	    factors[f].facs[0]=primes[f]*factors[n].facs[0];
	    for(int i=1;i<factors[n].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i];
		factors[f].facs[i]=factors[n].facs[i];
	      }
	  }
	else
	  {
	    factors[f].num_facs=factors[n].num_facs+1;
	    factors[f].facs[0]=primes[f];
	    factors[f].primes[0]=primes[f];
	    for(int i=1;i<factors[f].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i-1];
		factors[f].facs[i]=factors[n].facs[i-1];
	      }
	  }
      }
  free(primes);
  printf("Factors computed.\n");
  // now calculate phi(f)
  for(int i=3;i<=q_end;i++)
    {
      factors[i].phi=1;
      for(int j=0;j<factors[i].num_facs;j++)
	factors[i].phi*=phi(factors[i].primes[j],factors[i].facs[j]);
    }
  printf("phi computed.\n");

  //now do the prim roots
  factors[3].pr=2;
  factors[4].pr=3;
  //uint64_t max_pr=3;
  for(int i=5;i<=q_end;i++)
    {
      if(((factors[i].num_facs==1)&&(factors[i].primes[0]!=2))|| // p^m, p an odd prime
	 ((factors[i].num_facs==2)&&(factors[i].facs[0]==2)))    // 2p^m
	factors[i].pr=pr(i,factors);
      else
	factors[i].pr=0;
    }
  printf("pr's computed.\n");
  return(true);
}

bool squarefree(uint64_t d, factor *factors)
{
  for(uint64_t i=0;i<factors[d].num_facs;i++)
    if(factors[d].facs[i]!=factors[d].primes[i])
      return(false);
  return(true);
}


//
// kd(p)=(p-2)^{-1}
// multiplicative
// (kd(squarefull)=0 but we only pass in squarefree d)
// (kd(even)=0 but we don't pass in even d)
void kd(uint64_t d, factor *factors, mpfi_ptr res)
{
  mpfi_set_ui(res,1);
  for(uint64_t i=0;i<factors[d].num_facs;i++)
    mpfi_mul_ui(res,res,factors[d].primes[i]-2);
  mpfi_div(res,one,res);
}


//
// sums 2 kd(d) ed^{-k) H(d,k,N)
// for d in [d0..d1]
//
int main(int argc, char **argv)
{

  if(argc!=3)
    {
      printf("Usage:- %s <min_d> < max_d>\n",argv[0]);
      exit(0);
    }

  uint64_t d0=atol(argv[1]);
  uint64_t d1=atol(argv[2]);

  if((d0<3)||(d1<d0)||(d1<5)||(d1>MAX_D))
    {
      printf("min_d must be >2, max_d>min_d, max_d>4, max_d<=%lu. Exiting.\n",MAX_D);
      exit(0);
    }

  factor *factors=(factor *)malloc(sizeof(factor)*(d1+1));

  if(!make_factors(factors,d1))
    {
      printf("Error building factor database. Exiting.\n");
      exit(0);
    }

  mpfr_set_default_prec(200);
  mpfi_init(one);
  mpfi_set_ui(one,1);
  mpfi_t this_kd,temp_i,run_tot;
  mpfi_init(this_kd);
  mpfi_init(temp_i);
  mpfi_init(run_tot);
  mpfi_set_ui(run_tot,0);

  for(uint64_t d=d0;d<=d1;d++)
    {
      if(!(d&1))
	continue;
      if(!squarefree(d,factors))
	continue;
      kd(d,factors,this_kd);
      mpfi_div_ui(temp_i,this_kd,d);
      mpfi_add(run_tot,run_tot,temp_i);
    }
  mpfi_mul_ui(run_tot,run_tot,2);
  printf("\nFinal = ");
  mpfi_out_str(stdout,10,0,run_tot);
  printf("\n");

  return(0);
}
