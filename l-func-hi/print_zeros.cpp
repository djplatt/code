//
// print_zeros.cpp
//
// print zeros from file
// Just to demonstrate file format

#include <iostream>
#include <cstdlib>
#include "inttypes.h"
#include "mpfi.h"
#include "mpfi_io.h"

using namespace std;

#define ZPREC (102) // -log_2 working precision of zeros

#define fname (argv[1])
#define SUCCESS (0)
#define FAILURE (-1)

int in_bytes(mpfi_ptr t, FILE *infile)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
  int res;

  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    return(FAILURE);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    return(FAILURE);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    return(FAILURE);
  mpfi_set_ui(t,c);
  mpfi_mul_2ui(t,t,32);
  mpfi_add_ui(t,t,b);
  mpfi_mul_2ui(t,t,64);
  mpfi_add_ui(t,t,a);
  mpfi_div_2ui(t,t,ZPREC-1);
  return(SUCCESS);
}


// read the gap to the next zero into del_t
inline int next_delta(mpfi_ptr del_t, FILE *infile)
{
  return(in_bytes(del_t,infile));
}


int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>.\n",argv[0]);
      exit(0);
    }

  FILE *zfile=fopen(fname,"r");

  if(!zfile)
    {
      printf("Failed to open %s for binary input. Exiting.\n",fname);
      exit(0);
    }

  mpfr_set_default_prec(150);
  mpfi_t rho1,rho,del,pm1;
  mpfi_init(rho1);mpfi_init(rho);mpfi_init(del);
  mpfi_t tmp;
  mpfi_init(tmp);
  mpfi_init(pm1);
  mpfi_set_ui(pm1,1);
  mpfi_neg(tmp,pm1);
  mpfi_put(pm1,tmp);
  mpfi_div_2ui(pm1,pm1,ZPREC);
  mpfi_clear(tmp);

  uint64_t q,index,num_zeros;

  if(fread(&q,sizeof(uint64_t),1,zfile)!=1)
    {
      fprintf(stderr,"Error reading q from %s. Exiting.\n",fname);
      exit(0);
    }

  printf("Modulus %lu\n",q);

  while(fread(&index,sizeof(uint64_t),1,zfile)==1)
    {
      printf("Index %lu\n",index);
      if(fread(&num_zeros,sizeof(uint64_t),1,zfile)!=1)
	{
	  fprintf(stderr,"Error reading num_zeros for index %lu file %s. Exiting.\n",index,fname);
	  exit(0);
	}
      // all zeros files start from 0, except q=3 which starts from 8
      if(q==3)
	mpfi_set_ui(rho,8);
      else
	mpfi_set_ui(rho,0);
      for(uint64_t z=0;z<num_zeros;z++)
	if(next_delta(del,zfile)!=SUCCESS)
	  {
	    fprintf(stderr,"Error reading zero number %lu for index %lu from file %s. Exiting.\n",z,index,fname);
	    exit(0);
	  }
	else
	  {
	    mpfi_add(rho,rho,del); // this is exact
	    mpfi_add(rho1,rho,pm1); // this is +/- 2^{-ZPREC}
	    printf("%lu ",z);
	    mpfi_out_str(stdout,10,0,rho1);
	    printf("\n");
	  }
    }
  mpfi_clear(del);mpfi_clear(rho);mpfi_clear(rho1);
  mpfi_clear(pm1);
  return(0);
}
