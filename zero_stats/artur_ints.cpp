// find where sum 1/gamma crosses integers
// reads a b e from stdin and starts at [a,b]*2^e
#include "stdio.h"
#include "stdlib.h"
#include "acb.h"
#include "../trace/quad.h"
#include "../includes/arb_sum.h"

#define OP_ACC (101)
#include "inttypes.h"

arb_t pm1; // +/- 1 in the last place

// initialise some variables
void init_in_bytes()
{
  arb_t tmp;
  arb_init(tmp);
  arb_init(pm1);
  arb_set_ui(tmp,1);
  arb_mul_2exp_si(tmp,tmp,-OP_ACC-1);
  arb_zero(pm1);
  arb_add_error(pm1,tmp);
  arb_clear(tmp);
}

// read a 13 byte number from file
// structured 8,4,1
// read as if its exact
void in_bytes(arb_ptr t, FILE *infile, int64_t prec)
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
  arb_set_ui(t,c);
  arb_mul_2exp_si(t,t,32);
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64);
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-OP_ACC);
}


// read the next imaginary part of rho into del_t
inline void next_rho(arb_t del_t, FILE *infile, int64_t prec)
{
  in_bytes(del_t,infile,prec);
}


// Yannick Saouter
inline void do_rho(arb_t g1, arb_t gamma, int64_t prec)
{
  //printf("gamma = ");arb_printd(gamma,30);printf("\n");
  arb_inv(g1,gamma,prec);
}

int main(int argc, char **argv)
{
  int i;
  int64_t prec;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=4)
    {
      printf("Fatal error in main: usage %s <prec> <zeros file> <next int>. Exiting.\n",argv[0]);
      exit(0);
    }
  prec=atol(argv[1]);
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }
  uint64_t ni=atol(argv[3]);
  fmpz_t a,b,e;
  fmpz_init(a);
  fmpz_init(b);
  fmpz_init(e);
  arb_t tot1;arb_init(tot1);
  arb_sum(tot1,prec,a,b,e);
  fmpz_clear(a);
  fmpz_clear(b);
  fmpz_clear(e);

  printf("Start with sum 1/gamma = ");arb_printd(tot1,20);printf("\n");
  arb_t gamma,g1,del_t,t,tmp;
  arb_init(gamma);arb_init(del_t);arb_init(t);
  arb_init(g1);arb_init(tmp);



  
  init_in_bytes();


  long int num_its,it,z;
  double st[2];
  long int zs[2];
  fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile); // starting/ending t, exact
      fread(&zs[0],sizeof(long int),1,infile); // starting zero number
      if(st[0]==0.0)
	continue;
      fread(&zs[1],sizeof(long int),1,infile); // ending zero number

      arb_set_d(gamma,st[0]);
      arb_set(t,gamma);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile,prec); // distance to next gamma
          if(arb_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  arb_add(t,t,del_t,prec); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  arb_add(gamma,t,pm1,prec); // not exact

	  // compute Chris's sum at this gamma
	  do_rho(g1,gamma,prec);
	  // add the results to our running totals
	  arb_add(tot1,tot1,g1,prec);
	  arb_sub_ui(tmp,tot1,ni,prec);
	  if(arb_is_positive(tmp))
	    {
	      printf("Passed through %lu with zero %lu\n",ni,z);
	      ni++;
	      continue;
	    }
	  if(arb_is_negative(tmp))
	    continue;
	  printf("Error at zero %lu, can't tell which side of %lu sum is at. Exiting.\n",z,ni);
	  exit(0);	  
	}
    }
  return 0;
}
