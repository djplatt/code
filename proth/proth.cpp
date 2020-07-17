#include "stdlib.h"
#include "gmp.h"

#define NUM_TRIES (10)

int main(int argc, char **argv)
{

  if(argc!=2)
    exit(0);

  unsigned long int as[NUM_TRIES]={2,3,5,7,11,13,17,19,23,29};
  mpz_t a,p,e,res,twon,limit;
  mpz_init(a);
  mpz_init(p);
  mpz_init(e);
  mpz_init(res);
  mpz_init(twon);
  mpz_init(limit);

  mpz_set_ui(p,2);
  mpz_pow_ui(twon,p,atol(argv[1]));
  mpz_pow_ui(limit,p,128);

  for(unsigned long int i,j=1;;j+=2)
    {
      if(mpz_cmp_ui(twon,j)<=0)
	break;
      mpz_mul_ui(p,twon,j);
      mpz_div_ui(e,p,2);
      mpz_add_ui(p,p,1);
      if(mpz_cmp(limit,p)<=0)
	break;
      for(i=0;i<NUM_TRIES;i++)
	{
	  mpz_set_ui(a,as[i]);
	  mpz_powm(res,a,e,p);
	  mpz_add_ui(res,res,1);
	  if(mpz_cmp(res,p)==0)
	    {
	      mpz_out_str(NULL,10,p);printf(" is prime, used Proth base %lu\n",as[i]);
	      /*
	      if(mpz_probab_prime_p(p,10)==0)
		printf("But gmp thinks its composite.\n");
	      else
		printf("and gmp agrees.\n");
	      */
	      break;
	    }
	}
      /*
      if(i==NUM_TRIES)
	{
	  mpz_out_str(NULL,10,p);printf(" is not a Proth prime.\n");
	  if(mpz_probab_prime_p(p,10)==0)
	    printf("gmp thinks its composite.\n");
	  else
	    printf("gmp thinks it is (probably) prime.\n");
	}
      */
    }
  
  return (0);
}
