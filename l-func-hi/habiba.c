//
// habiba.c
//

#include "stdio.h"
#include "stdlib.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"

#define OP_ACC (101)
#include "inttypes.h"

inline void mpfi_print(mpfi_ptr x)
{
    mpfi_out_str(stdout,10,0,x);
    printf("\n");
};

inline void mpfi_print_str(const char *str, mpfi_ptr x)
{
  printf("%s",str);
  mpfi_print(x);
}

mpfi_t pm1;
void init_in_bytes()
{
  mpfi_t tmp;
  mpfi_init(tmp);
  mpfi_init(pm1);
  mpfi_set_ui(pm1,1);
  mpfi_neg(tmp,pm1);
  mpfi_put(pm1,tmp);
  mpfi_div_2ui(pm1,pm1,OP_ACC+1);
  mpfi_clear(tmp);
}

void in_bytes(mpfi_ptr t, FILE *infile)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
  int res;

  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  mpfi_set_ui(t,c);
  mpfi_mul_2ui(t,t,32);
  mpfi_add_ui(t,t,b);
  mpfi_mul_2ui(t,t,64);
  mpfi_add_ui(t,t,a);
  mpfi_div_2ui(t,t,OP_ACC);
}


// read the next imaginary part of rho into del_t
inline void next_rho(mpfi_t del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}

void read_null(FILE *infile)
{
  uint64_t a;
  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading null zero a. Exiting.\n");
      exit(0);
    }
  if(a!=0)
    {
      printf("Error reading null zero a neq 0. Exiting.\n");
      exit(0);
    }
  uint32_t b;
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Error reading null zero b. Exiting.\n");
      exit(0);
    }
  if(b!=0)
    {
      printf("Error reading null zero b neq 0. Exiting.\n");
      exit(0);
    }
  uint8_t c;
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Error reading null zero c. Exiting.\n");
      exit(0);
    }
  if(c!=0)
    {
      printf("Error reading null zero c neq 0. Exiting.\n");
      exit(0);
    }
}

int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>.\n",argv[0]);
      exit(0);
    }

  mpfr_set_default_prec(200);
  init_in_bytes();
  mpfi_t sum,del_t,t,rho,one,recip;
  mpfi_init(sum);
  mpfi_init(del_t);
  mpfi_init(t);
  mpfi_init(rho);
  mpfi_init(one);
  mpfi_set_ui(one,1);
  mpfi_init(recip);

  FILE *infile=fopen(argv[1],"r");

  if(!infile)
    {
      printf("Failed to open %s for input. Exiting.\n",argv[1]);
      exit(0);
    }


  uint64_t q;

  if(fread(&q,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading q from %s. Exiting.\n",argv[1]);
      exit(0);
    }

  uint64_t index;
  while(fread(&index,sizeof(uint64_t),1,infile)==1)
    {
      uint64_t num_zeros;
      if(fread(&num_zeros,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading num_zeros for index %lu file %s. Exiting.\n",index,argv[1]);
	  exit(0);
	}
      mpfi_set_ui(sum,0);
      mpfi_set_ui(t,0);
      uint64_t i,count1=0,count2=0;
      int c1=(1==1),c2=(1==1);
      for(i=0;i<num_zeros;i++)
	{
	  next_rho(del_t,infile);
	  mpfi_add(t,t,del_t);
	  mpfi_add(rho,t,pm1);
	  //mpfi_print_str("",rho);	  

	  mpfi_div(recip,one,rho);
	  mpfi_add(sum,sum,recip);
	  if(c1)
	    {
	      int cmp=mpfi_cmp_ui(rho,1);
	      if(cmp<0)
		count1++;
	      if(cmp==0)
		{
		  printf("Error, we have a zero straddling 1. Exiting.\n");
		  exit(0);
		}
	      if(cmp>0)
		c1=(1==0);
	    }
	  if(c2)
	    {
	      int cmp=mpfi_cmp_ui(rho,2);
	      if(cmp<0)
		count2++;
	      if(cmp==0)
		{
		  printf("Error, we have a zero straddling 2. Exiting.\n");
		  exit(0);
		}
	      if(cmp>0)
		c2=(1==0);
	    }
	}
      int cmp=mpfi_cmp_ui(rho,200);
      if(cmp==0)
	{
	  mpfi_print_str("Error, last zero straddles 200. Exiting. ",rho);
	  exit(0);
	}
      if(cmp>0)
	{
	  printf("q: %lu index: %lu %lu %lu %lu ",q,index,count1,count2,num_zeros);
	  mpfi_print(sum);
	  mpfi_print_str("Error, last zero exceeds 200. Exiting. ",rho);
	  exit(0);
	}


      printf("q: %lu index: %lu %lu %lu %lu ",q,index,count1,count2,num_zeros);
      mpfi_print(sum);
    }
  return(0);
}
