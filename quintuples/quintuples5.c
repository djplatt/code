/*
 *
 * Count how many triples of the form {a,b,d} with a<b<d, b^6<=d<=b^7.7
 * can be extended to {a,b,c,d} with b<c<d
 * We expect the answer to be zero
 */
#include "inttypes.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "gmp.h"

#define MAX_B (1300000000L)
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

mpz_t tmp1,tmp2,tmp3,ab1,r,t,s,d,upb,lwb;

mpz_t lwbc,upbc,tc,sc,rc,c;


uint64_t abcalls=0;

// c <-- (x^2-1)/a
void set_c(mpz_t x, uint64_t a, mpz_t c)
{
  mpz_mul(tmp1,x,x);
  mpz_sub_ui(tmp1,tmp1,1);
  mpz_div_ui(c,tmp1,a);
}

#define C_LIM (20)
mpz_t cs[C_LIM];

//
void next_xy(mpz_t s, mpz_t t, uint64_t a, uint64_t b, uint64_t ur)
{
  mpz_mul_ui(tmp1,t,ur);
  mpz_mul_ui(tmp2,s,b);
  mpz_add(tmp3,tmp1,tmp2); // t*r+s*b
  mpz_mul_ui(tmp1,t,a);
  mpz_mul_ui(tmp2,s,ur);
  mpz_add(s,tmp1,tmp2); // t*a+r*s
  mpz_set(t,tmp3);
  //printf("x,y set to ");mpz_out_str(NULL,10,s);printf(",");
  //mpz_out_str(NULL,10,t);printf("\n");
}

uint64_t try_it(uint64_t a, uint64_t b, uint64_t ur, int64_t x, int64_t y)
{
  printf("In try_it with %lu %lu %lu %ld %ld\n",a,b,ur,x,y);
  uint64_t count=0;
  mpz_ui_pow_ui(lwb,b,5);
  mpz_ui_pow_ui(upb,b,8); // quicker than b^7.7 ?
  mpz_set_si(s,x);
  mpz_set_si(t,y);
  uint64_t cptr=0;
  set_c(s,a,cs[0]);
  while(mpz_cmp(cs[cptr],lwb)<0)
    {
      printf("found candidate c=");
      mpz_out_str(NULL,10,cs[cptr]);
      printf("\n");
      next_xy(s,t,a,b,ur);
      cptr++;
      if(cptr==C_LIM)
	{
	  printf("Fatal error. Must increase C_LIM. Exiting.\n");
	  exit(0);
	}
      set_c(s,a,cs[cptr]);
    }
  mpz_set(d,cs[cptr]);
  while(mpz_cmp(d,upb)<=0)
    {
      uint64_t i;
      for(i=0;i<cptr;i++)
	{
	  if(mpz_cmp_ui(cs[i],b)<=0)
	    continue;
	  
	  printf("testing r=%lu a=%lu b=%lu c=",ur,a,b);
	  mpz_out_str(NULL,10,cs[i]);
	  printf(" d=");
	  mpz_out_str(NULL,10,d);
	  printf("\n");
	  
	  mpz_mul(tmp1,d,cs[i]);
	  mpz_add_ui(tmp1,tmp1,1);
	  if(mpz_perfect_square_p(tmp1))
	    {
	      	      
	      printf("solution found with r=%lu a=%lu b=%lu c=",ur,a,b);
	      mpz_out_str(NULL,10,cs[i]);
	      printf(" d=");
	      mpz_out_str(NULL,10,d);
	      printf("\n");
	      
	      count++;
	    }
	}
      next_xy(s,t,a,b,ur);
      set_c(s,a,d);
    }
  return(count);
}

uint64_t brute_force(uint64_t a, uint64_t b, uint64_t ur)
{
  abcalls++;
  uint64_t count=0;
  double S=b;
  S/=a;
  S=(S-1)*(ur-1);
  S=sqrt(S);
  uint64_t L1=floor(S),D=a*b,N=b*(b-a),y;
  count+=try_it(a,b,ur,1,1);
  count+=try_it(a,b,ur,1,-1);
  //return(count);
  for(y=0;y<=L1;y++)
    {
      uint64_t Sq=N+D*y*y;
      uint64_t x=sqrt((double) Sq);
      if((x*x==Sq)&&(x%b==0))
	{
	  count+=try_it(a,b,ur,x/b,y);
	  count+=try_it(a,b,ur,x/b,-y);
	  printf("%lu %lu %lu (%lu,%lu)\n",ur,a,b,x/b,y);
	}
    }
  return(count);
}
/*
uint64_t brute_force(uint64_t a, uint64_t b)
{
  abcalls++;
  uint64_t count=0;
  double S=b;
  S/=a;
  S=(S-1)*(ur-1);
  S=sqrt(S);
  uint64_t L1=floor(S),D=a*b,N=b*(b-a),y;
  count+=try_it(a,b,ur,1,1);
  count+=try_it(a,b,ur,1,-1);
  //return(count);
  for(y=0;y<=L1;y++)
    {
      uint64_t Sq=N+D*y*y;
      uint64_t x=sqrt((double) Sq);
      if((x*x==Sq)&&(x%b==0))
	{
	  count+=try_it(a,b,ur,x/b,y);
	  count+=try_it(a,b,ur,x/b,-y);
	  printf("%lu %lu %lu (%lu,%lu)\n",ur,a,b,x/b,y);
	}
    }
  return(count);
}
*/


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

  uint64_t i;
  for(i=0;i<C_LIM;i++)
    mpz_init(cs[C_LIM]);  

  //brute_force(21,37,223);exit(0);


  uint64_t r0=0,r1;

  uint64_t ur,uab1,nsols=0;
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
	  if((b<=MAX_B)&&(b>a+2)&&(b<a+a)) nsols+=brute_force(a,b,ur); 
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

  printf("We found %lu solutions for r in [%lu,%lu].\nWe tried %lu {a,b} pairs\n",nsols,r0,r1,abcalls);
  return(0);
}
