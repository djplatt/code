#include "acb_dirichlet.h"

#include "stdio.h"
#include "stdlib.h"

#define OP_ACC (101)
#include "inttypes.h"

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


// read the delta between zeros into del_t
inline void next_rho(arb_t del_t, FILE *infile, int64_t prec)
{
  in_bytes(del_t,infile,prec);
}

inline void do_rho(arb_t gamma, int64_t prec)
{
  //printf("Here I would do something with the zero at ");
  //arb_printd(gamma,20);
  //printf("\n");
  static acb_struct *zeta;
  static acb_t s;
  static arb_t smallest,x;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(x);
      arb_init(smallest);
      arb_set_ui(smallest,1000);
      zeta=(acb_struct *)malloc(sizeof(acb_t)*2);
      acb_init(zeta);acb_init(zeta+1);acb_init(s);
      arb_set_d(acb_realref(s),0.5);
    }

  arb_set(acb_imagref(s),gamma);

  acb_dirichlet_zeta_jet(zeta,s,0,2,40);

  if((!arb_contains_zero(acb_realref(zeta)))||(!arb_contains_zero(acb_imagref(zeta))))
    {
      printf("Error, no zero of zeta at ");
      acb_printd(s,40);
      printf("\n");
      exit(0);
    }

  acb_abs(x,zeta+1,100);
  if(arb_overlaps(x,smallest))
    {
      printf("Smallest so far and this zeta' overlap.\n");
      acb_printd(s,40);
      printf("\n");
      exit(0);
    }
  if(arb_lt(x,smallest))
    {
      printf("New smallest |zeta'| = ");
      arb_printd(x,40);
      printf("\nAt gamma = ");
      arb_printd(gamma,40);
      printf("\n");
      arb_set(smallest,x);
    }

  
}

int main(int argc, char **argv)
{
  uint64_t i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Fatal error in %s: usage %s <prec> <zeros file>. Exiting.\n",argv[0],argv[0]);
      exit(0);
    }
  int64_t prec=atol(argv[1]);
  printf("ARB working precision set to %ld\n",prec);
  FILE *infile=fopen(argv[2],"r");
  if(!infile)
    {
      printf("Fatal error in %s: failed to open zeros file %s for binary input. Exiting.\n",argv[0],argv[2]);
      exit(0);
    }

  arb_t gamma,pm1,del_t,t,max_so_far;
  acb_t res;
  arb_init(gamma);arb_init(pm1);arb_init(del_t);arb_init(t);
  arb_init(max_so_far);
  acb_init(res);
  arb_set_ui(gamma,1);
  arb_mul_2exp_si(gamma,gamma,-OP_ACC-1); // 2^{-102}
  arb_zero(pm1);
  arb_add_error(pm1,gamma);

  long int num_its,it,z;
  double st[2];
  long int zs[2],n_zeros=0;
  fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile); // starting/ending t, exact
      fread(&zs[0],sizeof(long int),1,infile); // starting zero number
      if(st[0]==0.0)
	continue;
      fread(&zs[1],sizeof(long int),1,infile); // ending zero number
      //printf("processing zeros %ld to %ld inclusive\n",zs[0]+1,zs[1]);
      n_zeros+=zs[1]-zs[0];
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

	  //printf("Zeros %lu ",z);arb_printd(gamma,20);printf(" ");
	  do_rho(gamma,prec);
	}
    }
  return 0;
}


