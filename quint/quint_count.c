/*
 *
 * Count how many double of the form {a,b} with 1<=a<a+2<b<a+a
 * have ab+1 a perfect square
 */
#include "inttypes.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

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


#define MAX_B (5688611L)
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

// read factorisations from factors.gp
int old_read_factorisation(FILE *infile, uint64_t *r, factor_t *f)
{
  if(fscanf(infile,"%lu %lu",r,&(f->num_facs))!=2)
    return(0);
  uint64_t ptr;
  for(ptr=0;ptr<f->num_facs;ptr++)
    {
      if(fscanf(infile,"%lu %lu",f->primes+ptr,f->powers+ptr)!=2)
	{
	  printf("Failed to parse factorisation of %lu. Exiting.\n",r[0]);
	  exit(0);
	}
      uint64_t i;
      f->facs[ptr]=f->primes[ptr];
      for(i=1;i<f->powers[ptr];i++)
	f->facs[ptr]*=f->primes[ptr];
    }
  fscanf(infile,"\n");
  return(1);
}

// read factorisations from factors.gp

int read_factorisation(FILE *infile, uint64_t *r, factor_t *f)
{
  if(fscanf(infile,"%lu %lu",r,&(f->num_facs))!=2)
    return(0);
  uint64_t fptr,pptr,last_p=0,p,pows,ffacs=f->num_facs;
  for(fptr=0,pptr=0;fptr<ffacs;fptr++)
    {
      if(fscanf(infile,"%lu %lu",&p,&pows)!=2)
	{
	  printf("Failed to parse factorisation of %lu. Exiting.\n",r[0]);
	  exit(0);
	}

      if(p==last_p)
	{
	  f->powers[pptr-1]+=pows;
	  uint64_t i;
	  for(i=0;i<pows;i++)
	    f->facs[pptr-1]*=p;
	  f->num_facs--;
	}
      else
	{
	  f->primes[pptr]=p;
	  f->powers[pptr]=pows;
	  f->facs[pptr]=p;
	  uint64_t i;
	  for(i=1;i<pows;i++)
	    f->facs[pptr]*=p;
	  last_p=p;pptr++;
	}
    }
  fscanf(infile,"\n");
  return(1);
}
      

void print_factors(factor_t *f, uint64_t n)
{
  printf("%lu factors as\n",n);
  uint64_t i=0;
  for(i=0;i<f->num_facs;i++)
    printf("  %lu^%lu\n",f->primes[i],f->powers[i]);
}


int main(int argc, char **argv)
{
  if(argc>2)
    {
      printf("Usage:- %s [r factor file].\n",argv[0]);
      exit(0);
    }

  FILE *rfile;
  if(argc==2)
    {
      rfile=fopen(argv[1],"r");
      if(!rfile)
	{
	  printf("Error opening %s for input. Exiting.\n");
	  exit(0);
	}
    }
  else
    rfile=stdin;


  uint64_t r0=0,r1;

  uint64_t ur,uab1,nsols=0;
  factor_t f;

  while(read_factorisation(rfile,&ur,&f))
    {
      //print_factors(&f,ur*ur-1);
      if(r0==0) r0=ur;
      r1=ur;
      uab1=ur*ur-1;
      uint64_t a=1,b=uab1;
      uint64_t this_pow[MAX_FACS];
      uint64_t i;
      for(i=0;i<f.num_facs;i++)
	this_pow[i]=0;
      while(1==1)
	{
	  if((b<=MAX_B)&&(b>a+2)&&(b<a+a)) nsols++;
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

  printf("We found %lu {a,b} pairs for r in [%lu,%lu]\n.",nsols,r0,r1);
  return(0);
}
