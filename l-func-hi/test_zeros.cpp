#include <iostream>
#include <cstdlib>
#include "characters.h"

#define PREC (300)
#define OP_ACC (101)

using namespace std;

uint64_t get_u64(FILE *infile)
{
  uint64_t res;
  if(fread(&res,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading u64. Exiting.\n");
      exit(0);
    }
  return(res);
}
uint32_t get_u32(FILE *infile)
{
  uint32_t res;
  if(fread(&res,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Error reading u32. Exiting.\n");
      exit(0);
    }
  return(res);
}
uint8_t get_u8(FILE *infile)
{
  uint8_t res;
  if(fread(&res,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Error reading u8. Exiting.\n");
      exit(0);
    }
  return(res);
}

int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>.\n",argv[0]);
      exit(0);
    }

  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Error opening %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  uint64_t q=get_u64(infile);
  printf("q=%lu\n",q);

  dft_init(PREC);
  DirichletGroup G(q);

  uint64_t index;
  mpz_t t,dt;
  mpz_init(t);
  mpz_init(dt);
  if(q==3)
    {
      mpz_set_ui(t,8);
      mpz_mul_2exp(t,t,OP_ACC);
    }
  else
    mpz_set_ui(t,0);
  mpfr_t(rt);
  mpfr_init(rt);
  while(fread(&index,sizeof(uint64_t),1,infile)==1)
    {
      printf("Index=%lu\n",index);
      uint64_t num_zeros=get_u64(infile);
      if(num_zeros==0)
	{
	  printf("No zeros for this index. Continuing.\n");
	  continue;
	}
      printf("There are %lu zeros\n",num_zeros);
      uint64_t conj_index=InvMod(index,q);
      uint64_t a;uint32_t b;uint8_t c;
      if(conj_index==index) // real character
	{
	  printf("Real Character.\n");
	  for(uint64_t n=0;n<num_zeros;n++)
	    {
	      a=get_u64(infile);
	      b=get_u32(infile);
	      c=get_u8(infile);

	      mpz_set_ui(dt,c);
	      mpz_mul_2exp(dt,dt,32);
	      mpz_add_ui(dt,dt,b);
	      mpz_mul_2exp(dt,dt,64);
	      mpz_add_ui(dt,dt,a);
	      mpz_add(t,t,dt);
	      mpfr_set_z(rt,t,GMP_RNDN);
	      mpfr_div_2ui(rt,rt,OP_ACC,GMP_RNDN);
	      mpfr_out_str(stdout,10,0,rt,GMP_RNDN);
	      printf("\n");
	    }
	  mpz_set_ui(t,0);
	}
      else
	{
	  printf("Complex character pair.\n");
	  for(uint64_t n=0;n<=num_zeros;n++)
	    {
	      a=get_u64(infile);
	      b=get_u32(infile);
	      c=get_u8(infile);
	      if((a==0)&&(b==0)&&(c==0))
		{
		  printf("End of first character.\n");
		  mpz_set_ui(t,0);
		  continue;
		}
	      mpz_set_ui(dt,c);
	      mpz_mul_2exp(dt,dt,32);
	      mpz_add_ui(dt,dt,b);
	      mpz_mul_2exp(dt,dt,64);
	      mpz_add_ui(dt,dt,a);
	      mpz_add(t,t,dt);
	      mpfr_set_z(rt,t,GMP_RNDN);
	      mpfr_div_2ui(rt,rt,OP_ACC,GMP_RNDN);
	      mpfr_out_str(stdout,10,0,rt,GMP_RNDN);
	      printf("\n");
	    }
	  mpz_set_ui(t,0);

	}
    }
	      
  return(0);
}
