#include "acb_poly.h"
#include "inttypes.h"

typedef struct{
  uint64_t conductor;
  uint64_t num_g_factors;
  uint64_t *g_factors;
  uint64_t denom;
  uint64_t num_polys;
  acb_poly_t *polys; // the inverted defining polynomials
  uint64_t num_primes;
  uint64_t *primes; // the primes
  uint64_t *poly_ptr; // which defining poly to use for this prime
  uint64_t poly_len; // number of coefficients carried in polys
  uint64_t num_coeffs;
  acb_t *coeffs; // the Dirichlet coefficients
} dirichlet_t;

void output_arf(arf_t x, uint64_t prec)
{
  static bool init=false;
  static fmpz_t m,e;
  if(!init)
    {
      init=true;
      fmpz_init(m);
      fmpz_init(e);
    }
  arf_get_fmpz_2exp(m,e,x);
  fmpz_print(m);
  printf(" ");
  fmpz_print(e);
}
  
// output arb as a double
void output_arb(arb_t x, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  static arf_t tmp1;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arf_init(tmp1);
    }
  arb_get_mid_arb(tmp,x);
  arb_get_ubound_arf(tmp1,tmp,prec);
  //printf("outputting centre = ");arf_printd(tmp1,10);printf("\n");
  printf("%20.18e",arf_get_d(tmp1,ARF_RND_NEAR));
}

void output_acb(acb_t z,uint64_t prec)
{
  output_arb(acb_realref(z),prec);
  printf(" ");
  output_arb(acb_imagref(z),prec);
  printf("\n");
}

void output_dirichlet_coefficients(dirichlet_t *d, uint64_t prec)
{
  for(uint64_t n=1;n<=d->num_coeffs;n++)
    output_acb(d->coeffs[n],prec);
}

bool sieve(dirichlet_t *d, uint64_t prec)
{
  acb_t tmp;
  acb_init(tmp);

  d->coeffs=(acb_t *)malloc(sizeof(acb_t)*(d->num_coeffs+1));
  for(uint64_t n=1;n<=d->num_coeffs;n++)
    {
      acb_init(d->coeffs[n]);
      acb_one(d->coeffs[n]);
    }
  for(uint64_t n=0;n<d->num_primes;n++)
    {
      uint64_t p=d->primes[n],pn=p,pnn=pn*p;
      uint64_t pow=1;
      //printf("Doing %lu\n",p);
      while(pn<=d->num_coeffs)
	{
	  uint64_t polyn=d->poly_ptr[n]-1; // we are 0 based, file is 1 based
	  acb_poly_get_coeff_acb(tmp,d->polys[polyn],pow);
	  //printf("doing %lu^%lu=%lu with ",p,pow,pn);acb_printd(tmp,10);printf("\n");
	  uint64_t ptr=pn;
	  while(ptr<=d->num_coeffs)
	    {
	      acb_mul(d->coeffs[ptr],d->coeffs[ptr],tmp,prec);
	      ptr+=pn;
	      if((ptr%pnn)==0)
		ptr+=pn;
	    }
	  pow++;
	  pn=pnn;
	  pnn*=p;
	}
    }

  acb_clear(tmp);
  return(true);
}
	      
void make_poly(acb_poly_t res, uint64_t num, uint64_t den, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2;
  static acb_t tmp;
  if(!init)
    {
      init=true;
      acb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
    }
  acb_poly_set_coeff_si(res,0,1);
  arb_set_ui(tmp1,num);
  arb_div_ui(tmp2,tmp1,den,prec);
  arb_mul_2exp_si(tmp2,tmp2,1);
  arb_sin_cos_pi(acb_imagref(tmp),acb_realref(tmp),tmp2,prec);
  acb_neg(tmp,tmp);
  acb_poly_set_coeff_acb(res,1,tmp);
}

bool parse_input(dirichlet_t *d, uint64_t prec)
{
  if(scanf("%lu\n",&d->conductor)!=1)
    {
      printf("Bad conductor. Exiting.\n");
      return(false);
    }
  uint64_t n_gamma;
  if(scanf("%lu\n",&d->num_g_factors)!=1)
    {
      printf("Bad number of gamma factors. Exiting.\n");
      return(false);
    }
  d->g_factors=(uint64_t *)malloc(sizeof(uint64_t)*d->num_g_factors);
  for(uint64_t n=0;n<d->num_g_factors;n++)
    if(scanf("%lu\n",d->g_factors+n)!=1)
      {
	printf("Bad gamma factor. Exiting.\n");
	return(false);
      }
  if(scanf("%lu\n",&d->denom)!=1)
    {
      printf("Bad denominator. Exiting.\n");
      return(false);
    }
  if(scanf("%lu\n",&d->num_polys)!=1)
    {
      printf("Bad number of polynomials. Exiting.\n");
      return(false);
    }
  d->polys=(acb_poly_t *)malloc(sizeof(acb_poly_t)*d->num_polys);
  acb_poly_t tmp;
  acb_poly_init(tmp);
  for(uint64_t n=0;n<d->num_polys;n++)
    {
      acb_poly_init(d->polys[n]);
      acb_poly_one(d->polys[n]);
      while(true)
	{
	  int64_t numer;
	  if(scanf("%ld",&numer)!=1)
	    {
	      printf("Bad numerator in polynomial. Exiting.\n");
	      return(false);
	    }
	  if(numer==-1) // it was []
	    break;
	  make_poly(tmp,numer,d->denom,prec);
	  acb_poly_mul(d->polys[n],d->polys[n],tmp,prec);
	}
    }
  if(scanf("%lu\n",&d->num_primes)!=1)
    {
      printf("Bad number of primes. Exiting.\n");
      return(false);
    }
  //printf("Number of primes=%lu\n",d->num_primes);
  d->primes=(uint64_t *)malloc(sizeof(uint64_t)*d->num_primes);
  d->poly_ptr=(uint64_t *)malloc(sizeof(uint64_t)*d->num_primes);
  for(uint64_t n=0;n<d->num_primes;n++)
    if(scanf("%lu %lu\n",d->primes+n,d->poly_ptr+n)!=2)
      {
	printf("Error reading prime/ptr pair. Exiting.\n");
	return(false);
      }
  // the input is 1 based, we are 0 based
  //for(uint64_t n=0;n<d->num_primes;n++)
  //d->poly_ptr[n]--;
  // how many coefficients do we need
  d->num_coeffs=d->primes[d->num_primes-1];
  d->poly_len=ceil(log((double)d->num_coeffs)/log(2.0));
  for(uint64_t n=0;n<d->num_polys;n++)
    {
      acb_poly_inv_series(d->polys[n],d->polys[n],d->poly_len,prec);
      //printf("Poly %lu\n",n);
      //acb_poly_printd(d->polys[n],10);
      //printf("\n");
    }
  return(true);
}




int main(int argc, char **argv)
{
  dirichlet_t d;
  if(!parse_input(&d,200))
    exit(0);
  if(!sieve(&d,200))
    exit(0);
  output_dirichlet_coefficients(&d,200);
  return(0);
}
