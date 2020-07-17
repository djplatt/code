/*
 *
 * Count how many triples of the form {a,b,d} with a<b<d, b^6<=d<=b^7.7
 * can be extended to {a,b,c,d} with b<c<d
 * We expect the answer to be zero
 */
#include "inttypes.h"
#include "stdio.h"
#include "stdlib.h"
#include "gmp.h"

#define MAX_B (1300000000L)
#define NPRIMES (65228333L)
uint64_t primes[NPRIMES];
#define PRIME_FILE "primes.txt"

void read_primes()
{
  FILE *pfile=fopen(PRIME_FILE,"r");
  uint64_t ptr;
  for(ptr=0;ptr<NPRIMES;ptr++)
    if(fscanf(pfile,"%lu\n",primes+ptr)!=1)
      {
	printf("Error reading primes.\n");
	exit(0);
      }
}


#define MAX_FACS (17) // 2*3*..59>49e18

// structure to hold factor information
typedef struct
{
  //uint64_t pr;   // = 0 if no primitive root
  //uint64_t phi;  // Euler phi(q)
	uint64_t num_facs;  // number of prime power factors
	uint64_t primes[MAX_FACS]; // the prime factors p of n
	uint64_t powers[MAX_FACS];  // the powers of p
	uint64_t facs[MAX_FACS]; // the p^power
} factor_t;

void print_factors(factor_t *f, uint64_t n)
{
  printf("%lu factors as\n",n);
  uint64_t i=0;
  for(i=0;i<f->num_facs;i++)
    printf("  %lu^%lu\n",f->primes[i],f->powers[i]);
}

void factor_me(factor_t *factors, uint64_t n)
{
  //printf("In factor_me\n");fflush(stdout);
  uint64_t prime=0,ptr=0;
  factors->num_facs=0;
  while(prime<NPRIMES)
    {
      //printf("prime=%lu\n",primes[prime]);
      if((n%primes[prime])==0)
	{
	  factors->num_facs++;
	  factors->primes[ptr]=primes[prime];
	  factors->powers[ptr]=1;
	  factors->facs[ptr]=primes[prime];
	  n/=primes[prime];
	  while((n>1)&&((n%primes[prime])==0))
	    {
	      n/=primes[prime];
	      factors->powers[ptr]++;
	      factors->facs[ptr]*=primes[prime];
	    }
	  if(n==1)
	    return;
	  ptr++;
	}
      prime++;
    }
}


mpz_t tmp1,tmp2,tmp3,ab1,r,t,s,d,upb,lwb;

mpz_t lwbc,upbc,tc,sc,rc,c;

uint64_t ncalls=0;

uint64_t abd_quick(uint64_t a,uint64_t b)
{
  printf("trying %lu %lu ",a,b);
  mpz_out_str(NULL,10,d);
  printf("\n");
  ncalls++;
  mpz_mul_ui(tmp1,d,b);
  mpz_add_ui(tmp1,tmp1,1);
  mpz_sqrt(rc,tmp1); // rc=sqrt(bd+1)
  //printf("rc=");mpz_out_str(NULL,10,rc);printf("\n");
  mpz_sub_ui(tmp1,rc,b);
  mpz_mul(tmp2,tmp1,tmp1);
  mpz_sub_ui(tmp2,tmp2,1);
  mpz_div_ui(c,tmp2,b);
  /*
  mpz_set_ui(tc,1); // starting column vector (tc,sc) is (1,1) or (1,-1)
  mpz_set_si(sc,-1);

  mpz_mul(tmp1,tc,rc);
  mpz_mul(tmp2,sc,d);
  mpz_add(tmp3,tmp1,tmp2);
  mpz_mul_ui(tmp1,tc,b);
  mpz_mul(tmp2,rc,sc);
  mpz_add(sc,tmp1,tmp2); // sc=tc*b+sc*rc
  mpz_set(tc,tmp3); // tc=tc*rc+sc*d
  mpz_mul(tmp1,sc,sc);
  mpz_sub_ui(tmp1,tmp1,1);
  mpz_div_ui(c,tmp1,b); // c=(sc^2-1)/b=(tc^2-1)/d
  */
  //printf("   trying with c=");
  //mpz_out_str(NULL,10,c);
  //printf("\n");
  mpz_mul_ui(tmp1,c,a);
  mpz_add_ui(tmp1,tmp1,1);
  if(mpz_perfect_square_p(tmp1)) // check ac+1 is a square
    {
      printf("%lu %lu ",a,b);
      mpz_out_str(NULL,10,c);
      printf(" ");
      mpz_out_str(NULL,10,d);
      printf(" works\n");
      return(1);
    }
  else
    return(0);
}
  

// Given a,b and d, find all c in (b,d) such that {a,b,c,d} is a D. quadruple
uint64_t abd(uint64_t a,uint64_t b, int64_t ss)
{
  if(ss>0)
    printf("trying (+) %lu %lu ",a,b);
  else
    printf("trying (-) %lu %lu ",a,b);
  mpz_out_str(NULL,10,d);
  printf("\n");
  ncalls++;
  uint64_t nsols=0;
  uint64_t lwbc=b+1;
  mpz_mul_ui(tmp1,d,b);
  mpz_add_ui(tmp1,tmp1,1);
  mpz_sqrt(rc,tmp1); // rc=sqrt(bd+1)
  mpz_sub_ui(upbc,d,1);
  mpz_set_ui(tc,1); // starting column vector (tc,sc) is (1,1) or (1,-1)
  mpz_set_si(sc,ss);
  while(1)
    {
      mpz_mul(tmp1,tc,rc);
      mpz_mul(tmp2,sc,d);
      mpz_add(tmp3,tmp1,tmp2);
      mpz_mul_ui(tmp1,tc,b);
      mpz_mul(tmp2,rc,sc);
      mpz_add(sc,tmp1,tmp2); // sc=tc*b+sc*rc
      mpz_set(tc,tmp3); // tc=tc*rc+sc*d
      mpz_mul(tmp1,sc,sc);
      mpz_sub_ui(tmp1,tmp1,1);
      mpz_div_ui(c,tmp1,b); // c=(sc^2-1)/b=(tc^2-1)/d
      if(mpz_cmp_ui(c,lwbc)>=0) // c big enough
	{
	  if(mpz_cmp(c,upbc)<=0) // c not too big
	    {
	      printf("   trying with c=");
	      mpz_out_str(NULL,10,c);
	      printf("\n");
	      mpz_mul_ui(tmp1,c,a);
	      mpz_add_ui(tmp1,tmp1,1);
	      if(mpz_perfect_square_p(tmp1)) // check ac+1 is a square
		{
		  printf("%lu %lu ",a,b);
		  mpz_out_str(NULL,10,c);
		  printf(" ");
		  mpz_out_str(NULL,10,d);
		  printf(" works\n");
		  nsols++;
		}
	    }
	  else // c too big
	    break;
	}
    }
  return(nsols);
}

// given (a,b) such that ab+1=r^2, find all d such that {a,b,d} is a triple
// with d in [b^6,b^7.7]
// NB Will fail to find all solutions if 
uint64_t abp(uint64_t a,uint64_t b,int64_t ss)
{
  uint64_t nsols=0;
  mpz_set_ui(t,1);
  mpz_set_si(s,ss); // vector [t,s]=[1,1] or [1,-1]
  mpz_ui_pow_ui(lwb,b,5);
  mpz_ui_pow_ui(upb,b,8); // quicker than b^7.7 ?
  printf("b^8=");mpz_out_str(NULL,10,upb);printf("\n");
  while(1)
    {
      mpz_mul(tmp1,t,r);
      mpz_mul_ui(tmp2,s,b);
      mpz_add(tmp3,tmp1,tmp2); // t*r+s*b
      mpz_mul_ui(tmp1,t,a);
      mpz_mul(tmp2,r,s);
      mpz_add(s,tmp1,tmp2); // t*a+r*s
      mpz_set(t,tmp3);
      mpz_mul(tmp1,s,s);
      mpz_sub_ui(tmp1,tmp1,1);
      mpz_div_ui(d,tmp1,a); //d=(s^2-1)/a=(t^2-1)/b
      //printf("found %lu %lu ",a,b);
      //mpz_out_str(NULL,10,d);
      //printf("\n");

      if(mpz_cmp(d,lwb)>=0)
	{
	  if(mpz_cmp(d,upb)<=0) // we have {a,b,d} with d in range
	    nsols+=abd_quick(a,b); // try to extend to {a,b,c,d}
	  // abd(a,b,1) would produce a c larger than d
	  else
	    break; // d too big, so break
	}
    }
  return(nsols);
}

uint64_t abcalls=0;
uint64_t ab(uint64_t a, uint64_t b, uint64_t ur)
{
  abcalls++;
  printf("%lu %lu %lu\n",a,b,ur);

  //printf("in ab with a=%lu b=%lu ur=%lu\n",a,b,ur);fflush(stdout);
  mpz_set_ui(r,ur);
  return(abp(a,b,1)+abp(a,b,-1));
}

      
int main(int argc, char **argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <r0> <r1>.\n",argv[0]);
      exit(0);
    }

  int64_t r0=atol(argv[1]);
  int64_t r1=atol(argv[2]);
  read_primes();
  //mpz_inits(tmp1,tmp2,tmp3,ab1,r,t,s,d,upb,lwb,NULL);
  
  mpz_init(tmp1);
  mpz_init(tmp2);
  mpz_init(tmp3);
  mpz_init(ab1);
  mpz_init(r);
  mpz_init(s);
  mpz_init(t);
  mpz_init(d);
  mpz_init(upb);
  mpz_init(lwb);
  mpz_init(upbc);
  mpz_init(lwbc);
  mpz_init(rc);
  mpz_init(sc);
  mpz_init(tc);
  mpz_init(c);  

  factor_t f;
  

  uint64_t ur,uab1,nsols=0;
  for(ur=r0;ur<=r1;ur++)
    {
      uab1=ur*ur-1;
      factor_me(&f,uab1);
      //print_factors(&f,uab1);
      uint64_t a=1,b=uab1;
      uint64_t this_pow[MAX_FACS];
      uint64_t i;
      for(i=0;i<f.num_facs;i++)
	this_pow[i]=0;
      while(1==1)
	{
	  if((b<=MAX_B)&&(b>a)&&(b<a+a)) nsols+=ab(a,b,ur); 
	  this_pow[0]++;
	  uint64_t ptr=0;
	  while(this_pow[ptr]>f.powers[ptr])
	    {
	      this_pow[ptr]=0;
	      a/=f.facs[ptr];
	      ptr++;
	      if(ptr==f.num_facs)
		break;
	      this_pow[ptr]++;
	    }
	  if(ptr==f.num_facs) break;
	  a*=f.primes[ptr];
	  b=uab1/a;
	}
    }

  printf("We found %lu solutions with r in [%ld,%ld].\nWe tried %lu {a,b} pairs\nWe tried %lu {a,b,d} triples.\n",nsols,r0,r1,abcalls,ncalls);
  return(0);
}
