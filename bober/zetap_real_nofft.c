#include "flint/acb_dirichlet.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define OP_ACC (101)
#include "inttypes.h"

#define ARB_PREC (140) // enough for 102 bits after dp of 10^11 size thing

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
void next_rho(arb_t del_t, FILE *infile, int64_t prec)
{
  in_bytes(del_t,infile,prec);
}



uint64_t count=0;
uint64_t neg_count=0;

// compute 1/zeta'(s) for s=1/2+i gamma
// checks that 1/2+i gamma is a zero of zeta
void do_rho(arb_t gamma, arb_t t0, int64_t prec)
{
  static acb_struct *zeta;
  static acb_t s,ctmp1,z;
  static bool init=false;
  if(!init)
    {
      init=true;
      zeta=(acb_struct *)malloc(sizeof(acb_t)*2);
      acb_init(zeta);acb_init(zeta+1);acb_init(s);acb_init(ctmp1);acb_init(z);
      arb_set_d(acb_realref(s),0.5);
    }

  
  arb_set(acb_imagref(s),gamma);
  acb_dirichlet_zeta_jet(zeta,s,0,2,prec/4);
  if((!arb_contains_zero(acb_realref(zeta)))||(!arb_contains_zero(acb_imagref(zeta))))
    {
      printf("Error, no zero of zeta at ");
      acb_printd(s,40);
      printf("\n");
      exit(0);
    }
  if(!arb_is_positive(acb_realref(zeta+1)))
    count++;
  else
    neg_count++;
}
  
int main(int argc, char **argv)
{
  uint64_t i;
  uint64_t prec=ARB_PREC;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  fflush(stdout);
  if(argc!=3)
    {
      printf("Fatal error in %s: usage %s <zeros file> <t0>. Exiting.\n",argv[0],argv[0]);
      exit(0);
    }
  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Fatal error in %s: failed to open zeros file %s for binary input. Exiting.\n",argv[0],argv[1]);
      exit(0);
    }

  int64_t t0=atol(argv[2]);
  if(t0<=0)
    {
      printf("Need t0 > 0. Exiting.\n");
      exit(0);
    }
  arb_t a_t0;
  arb_init(a_t0);
  arb_set_ui(a_t0,t0);
  
  arb_t gamma,pm1,del_t,t;
  //acb_t res;

  
  arb_init(gamma);arb_init(pm1);arb_init(del_t);arb_init(t);
  //acb_init(res);
  arb_set_ui(gamma,1);
  arb_mul_2exp_si(gamma,gamma,-OP_ACC-1); // 2^{-102}
  arb_zero(pm1);
  arb_add_error(pm1,gamma);

  long int num_its,it,z;
  double st[2];
  long int zs[2],n_zeros=0;
  int rval;
  rval=fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  bool done=false;
  for(it=0;it<num_its;it++)
    {
      if(done)
	break;
      rval=fread(st,sizeof(double),2,infile); // starting/ending t, exact
      rval=fread(zs,sizeof(long int),1,infile); // starting zero number
      if(st[0]==0.0)
	{
	  printf("Iteration %lu empty.\n",it);
	  continue;
	}
      rval=fread(zs+1,sizeof(long int),1,infile); // ending zero number
      //printf("processing zeros %ld to %ld = %ld inclusive\n",zs[0]+1,zs[1],zs[1]-zs[0]);
      n_zeros+=zs[1]-zs[0];
      arb_set_d(gamma,st[0]);
      arb_set(t,gamma);
      //printf("doing t from %f to %f zeros from %ld to %ld\n",st[0],st[1],zs[0],zs[1]);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  //printf("%ld\n",z);
	  next_rho(del_t,infile,prec); // distance to next gamma
          if(arb_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  arb_add(t,t,del_t,prec); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  arb_sub_ui(gamma,t,t0,prec);
	  if(arb_is_positive(gamma)) // gamma exceeds t0
	    {
	      done=true;
	      break;
	    }

	  arb_add(gamma,t,pm1,prec); // not exact
	  do_rho(gamma,a_t0,prec); // do something with this zero
	}
    }
  printf("There were %lu zeros where Re zeta'(rho)<0.\n",count);
  printf("There were %lu zeros where Re zeta'(rho)>=0.\n",neg_count);
  printf("Ratio %f\n",(double) count/ (double) (n_zeros));
  return 0;
  printf("\nWe processed %ld zeros.\n",n_zeros);
  
  return 0;
}


