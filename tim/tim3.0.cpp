//
// File: tim3.0.cpp
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
// We return the min found.
//
// version 3.0 exploits residue classes
//
#include "stdlib.h"
#include "stdio.h"
#include "malloc.h"
#include "inttypes.h"
#include "math.h"
#include "gmp.h"
#include "mpfi.h"
#include "mpfi_io.h"

#define MAX_D ((uint64_t) 100000) 
#define MAX_FACS ((uint64_t) 7) // 2*3*5*7*11*13*17>MAX_D

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

// rotate W_spare right by r places and add it into W component wise
// r is in [0,d)
void rotate_and_add(mpz_t *W, mpz_t *W_spare, uint64_t r, uint64_t d)
{
  uint64_t ptr=d-r;
  for(uint64_t i=0;i<d;i++)
    {
      mpz_add(W[i],W[i],W_spare[ptr++]);
      if(ptr==d)
	ptr=0;
    }
}

// see how many times sum_{i=1}^k 2^v_i hits each residue class of d
// where each v_i spans [1..ed]
// return minimum found
void doit(uint64_t k, uint64_t d, uint64_t ed, mpz_t *W, mpz_t *Wk)
{
  for(uint64_t i=0;i<d;i++)
    mpz_set_ui(W[i],0);
  // set up W_1

  for(uint64_t i=1,ptr=2*k;i<=ed;i++)
    {
      ptr=(ptr+(1<<i))%d;
      mpz_add_ui(W[ptr],W[ptr],1);
    }

  // now do W_2..W_K
  for(uint64_t i=2;i<=k;i++)
    {
      // copy W to Wk
      for(uint64_t j=0;j<d;j++)
	mpz_set(Wk[j],W[j]);
      //
      // now rotate Wk by 2,6,14,30... and add into W
      //
      for(uint64_t rl=2,drl=4,e=1;e<ed;e++)
	{
	  rotate_and_add(W,Wk,rl,d);
	  rl+=drl;
	  rl%=d;
	  drl<<=1;
	  drl%=d;
	}
    }
  /*  
  for(uint64_t i=0;i<d;i++)
    {
      printf("%lu:%lu ",d,i);
      mpz_out_str(NULL,10,W[i]);
      printf("\n");
    }
  */
}

//
// kd(p)=(p-2)^{-1}
// multiplicative
// (kd(squarefull)=0 but we only pass in squarefree d)
// (kd(even)=0 but we don't pass in even d)
mpfi_t one;
void kd(uint64_t d, factor *factors, mpfi_ptr res)
{
  mpfi_set_ui(res,1);
  for(uint64_t i=0;i<factors[d].num_facs;i++)
    mpfi_div_ui(res,res,factors[d].primes[i]-2);
}

// compute smallest n>0 s.t. 2^n = 1 mod d
// don't call this with even d!
uint64_t ord2(uint64_t d)
{
  uint64_t a=2,res=1;
  while(a!=1)
    {
      a<<=1;res++;
      if(a>=d)
	a-=d;
    }
  return(res);
}

//
// sums 2 kd(d) ed^{-k) H(d,k,N)
// for d in [d0..d1]
//
int main(int argc, char **argv)
{

  printf("Command line:-");
  for(uint64_t i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=5)
    {
      printf("Usage:- %s <K> <min_d> < max_d> <rclass>\n",argv[0]);
      exit(0);
    }

  if(atol(argv[1])<1)
    {
      printf("K must be > 0. Exiting.\n");
      exit(0);
    }

  uint64_t k=atol(argv[1]);
  uint64_t d0=atol(argv[2]);
  uint64_t d1=atol(argv[3]);
  uint64_t rclass=atol(argv[4]);
  mpz_t *W,*Wk,this_W;
  W=(mpz_t *)malloc(sizeof(mpz_t)*d1);
  Wk=(mpz_t *)malloc(sizeof(mpz_t)*d1);
  for(uint64_t i=0;i<d1;i++)
    {
      mpz_init(W[i]);
      mpz_init(Wk[i]);
    }
  mpz_init(this_W);
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
  mpfi_t *results=(mpfi_t *)malloc(sizeof(mpfi_t)*rclass);
  for(uint64_t i=0;i<rclass;i++)
    {
      mpfi_init(results[i]);
      mpfi_set_ui(results[i],0);
    }
  mpfi_t this_kd,temp_i,run_tot;
  mpfi_init(this_kd);
  mpfi_init(temp_i);
  mpz_t edk;
  mpz_init(edk);

  for(uint64_t d=d0,count=0;d<=d1;d++)
    {
      if(!(d&1))
	continue;
      if(!squarefree(d,factors))
	continue;
      count++;
      kd(d,factors,this_kd);
      uint64_t ed=ord2(d);
      //printf("Running with K=%lu d=%lu and Ord_%lu(2)=%lu.\n",k,d,d,ed);
      mpz_set_ui(edk,ed);
      for(uint64_t i=1;i<k;i++)
	mpz_mul_ui(edk,edk,ed);
      //printf("%lu^%lu=",ed,k);mpz_out_str(NULL,10,edk);printf("\n");
      doit(k,d,ed,W,Wk);
 
      // W now contains the counts
      uint64_t thisgcd=gcd(d,rclass);
      if(thisgcd==1) // nothing to be gained
	{
	  mpz_set(this_W,W[0]); // find smallest W
	  for(uint64_t i=1;i<d;i++)
	    if(mpz_cmp(W[i],this_W)<0)
	      mpz_set(this_W,W[i]);
	  //printf("d=%lu smallest W=",d);
	  //mpz_out_str(NULL,10,this_W);
	  //printf("\n");
	  mpfi_mul_z(temp_i,this_kd,this_W); // compute contribution
	  mpfi_div_z(temp_i,temp_i,edk);
	  mpfi_mul_2exp(temp_i,temp_i,1);
	  for(uint64_t i=0;i<rclass;i++)
	    mpfi_add(results[i],results[i],temp_i); // add it to every res class
	}
      else
	{
	  for(uint64_t i=0;i<thisgcd;i++)
	    {
	      uint64_t j=i;
	      mpz_set(this_W,W[j]);
	      j+=thisgcd;
	      while(j<d)
		{
		  if(mpz_cmp(W[j],this_W)<0)
		    mpz_set(this_W,W[j]);
		  j+=thisgcd;
		}
	      //printf("d=%lu,%lu smallest W=",d,i);
	      //mpz_out_str(NULL,10,this_W);
	      //printf("\n");

	      mpfi_mul_z(temp_i,this_kd,this_W);
	      mpfi_div_z(temp_i,temp_i,edk);
	      mpfi_mul_2exp(temp_i,temp_i,1);
	      j=i;
	      while(j<rclass)
		{
		  mpfi_add(results[j],results[j],temp_i);
		  j+=thisgcd;
		}
	    }
	}
      printf("after d=%lu results contains: ",d);
      mpfi_set(temp_i,results[0]);
      for(uint64_t i=1;i<rclass;i++)
	if(mpfi_cmp(results[i],temp_i)<0)
	  mpfi_set(temp_i,results[i]);
      mpfi_out_str(stdout,10,0,temp_i);
      printf("\n");
    }

  for(uint64_t i=1;i<rclass;i++)
    if(mpfi_cmp(results[i],results[0])<0)
      mpfi_set(results[0],results[i]);

  mpfi_add_ui(results[0],results[0],2);
  printf("Smallest was ");
  mpfi_out_str(stdout,10,0,results[0]);
  printf("\n");
  return(0);
}
