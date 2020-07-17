#include "stdio.h"
#include "stdlib.h"
#include "acb.h"
#include "../trace/quad.h"

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
inline void do_rho(arb_t g1, arb_t g3, arb_t g4, arb_t g5, arb_t g6, arb_t gamma, int64_t prec)
{
  //printf("gamma = ");arb_printd(gamma,30);printf("\n");
  arb_inv(g1,gamma,prec);
  arb_div(g4,g1,gamma,prec);
  arb_div(g3,g4,gamma,prec);
  arb_div(g4,g3,gamma,prec);
  arb_div(g5,g4,gamma,prec);
  arb_div(g6,g5,gamma,prec);
}

int main(int argc, char **argv)
{
  int i;
  int64_t prec=200;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Fatal error in main: usage %s <prec> <zeros file>. Exiting.\n",argv[0]);
      exit(0);
    }
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }

  arb_t gamma,g1,g3,g4,g5,g6,tot1,tot3,tot4,tot5,tot6,del_t,t;
  arb_init(gamma);arb_init(del_t);arb_init(t);
  arb_init(g1);arb_init(tot1);
  arb_init(g3);arb_init(tot3);
  arb_init(g4);arb_init(tot4);
  arb_init(g5);arb_init(tot5);
  arb_init(g6);arb_init(tot6);



  
  init_in_bytes();

  arb_set_ui(tot1,0);
  arb_set_ui(tot3,0);
  arb_set_ui(tot4,0);
  arb_set_ui(tot5,0);
  arb_set_ui(tot6,0);

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
	  do_rho(g1,g3,g4,g5,g6,gamma,prec);
	  // add the results to our running totals
	  arb_add(tot1,tot1,g1,prec);
	  arb_add(tot3,tot3,g3,prec);
	  arb_add(tot4,tot4,g4,prec);
	  arb_add(tot5,tot5,g5,prec);
	  arb_add(tot6,tot6,g6,prec);
	  
	}
    }
  printf("term by term we got\n1     ");
  arb_printd(tot1,20);
  fmpz_t a,b,e;
  fmpz_init(a);
  fmpz_init(b);
  fmpz_init(e);
  arb_get_interval_fmpz_2exp(a,b,e,tot1);
  printf(" ");
  fmpz_print(a);printf(" ");
  fmpz_print(b);printf(" ");
  fmpz_print(e);printf(" ");
  printf("\n");
  printf("3     ");
  arb_printd(tot3,20);
  arb_get_interval_fmpz_2exp(a,b,e,tot3);
  printf(" ");
  fmpz_print(a);printf(" ");
  fmpz_print(b);printf(" ");
  fmpz_print(e);printf(" ");
  printf("\n");
  printf("4     ");
  arb_printd(tot4,20);
  arb_get_interval_fmpz_2exp(a,b,e,tot4);
  printf(" ");
  fmpz_print(a);printf(" ");
  fmpz_print(b);printf(" ");
  fmpz_print(e);printf(" ");
  printf("\n");
  printf("5     ");
  arb_printd(tot5,20);
  arb_get_interval_fmpz_2exp(a,b,e,tot5);
  printf(" ");
  fmpz_print(a);printf(" ");
  fmpz_print(b);printf(" ");
  fmpz_print(e);printf(" ");
  printf("\n");
  printf("6     ");
  arb_printd(tot6,20);
  arb_get_interval_fmpz_2exp(a,b,e,tot6);
  printf(" ");
  fmpz_print(a);printf(" ");
  fmpz_print(b);printf(" ");
  fmpz_print(e);printf(" ");
  printf("\n");
  return 0;
}
