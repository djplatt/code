#include "arb.h"
#include "inttypes.h"
#include "stdlib.h"

#define SUCCESS ((int) 0)
#define FAILURE ((int) -1)
#define PREC ((int64_t) 120)

int read_delta(arb_t res, FILE *infile,int64_t prec)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
 if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    return(FAILURE);
 if(fread(&b,sizeof(uint32_t),1,infile)!=1)
   return(FAILURE);
 if(fread(&c,sizeof(uint8_t),1,infile)!=1)
   return(FAILURE);
 arb_set_ui(res,c);
 arb_mul_2exp_si(res,res,32);
 arb_add_ui(res,res,b,prec);
 arb_mul_2exp_si(res,res,64);
 arb_add_ui(res,res,c,prec);
 arb_mul_2exp_si(res,res,-102);
 return(SUCCESS);
}

int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>.\n",argv[0]);
      exit(0);
    }

  FILE *zfile=fopen(argv[1],"r");

  if(!zfile)
    {
      printf("Failed to open %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  arb_t rho1,rho,del,pm1;
  arb_init(rho1);arb_init(rho);arb_init(del);arb_init(pm1);

  arb_zero(pm1);
  arb_set_ui(del,1);
  arb_add_error(pm1,del);
  arb_mul_2exp_si(pm1,pm1,-102);

  uint64_t q,index,num_zeros,z;

  if(fread(&q,sizeof(uint64_t),1,zfile)!=1)
    {
      fprintf(stderr,"Error reading q from %s. Exiting.\n",argv[1]);
      exit(0);
    }


  while(fread(&index,sizeof(uint64_t),1,zfile)==1)
    {
      if(fread(&num_zeros,sizeof(uint64_t),1,zfile)!=1)
	{
	  fprintf(stderr,"Error reading num_zeros for index %lu file %s. Exiting.\n",index,argv[1]);
	  exit(0);
	}
      // all zeros files start from 0, except q=3 which starts from 8
      if(q==3)
	arb_set_ui(rho,8);
      else
	arb_set_ui(rho,0);

      uint64_t count=0;
      for(z=0;z<num_zeros;z++)
	if(read_delta(del,zfile,PREC)!=SUCCESS)
	  {
	    fprintf(stderr,"Error reading zero number %lu for index %lu from file %s. Exiting.\n",z,index,argv[1]);
	    exit(0);
	  }
	else
	  {
	    arb_add(rho,rho,del,PREC); // this is exact
	    arb_add(rho1,rho,pm1,PREC); // this is +/- 2^{-ZPREC}
	    //arb_printd(rho1,20);printf("\n");
	    arb_sub_ui(rho1,rho1,1,PREC);
	    if(arb_is_negative(rho1))
	      count++;
	    else
	      {
		if(!arb_is_positive(rho1))
		  {
		    fprintf(stderr,"Zero straddling 1. Exiting.\n");
		    exit(0);
		  }
	      }
	  }
      printf("%lu %lu %lu\n",q,index,count);
    }
  return(0);
}
