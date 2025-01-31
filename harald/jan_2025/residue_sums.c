#include "flint/acb_dirichlet.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

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


void eta(arb_t res, arb_t t, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_neg(tmp,t);
  arb_add_ui(res,tmp,2,prec);
}

void et_gam(arb_t res, arb_t gamma, uint64_t t0_by_2, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_div_ui(tmp,gamma,t0_by_2,prec);
  eta(res,tmp,prec);
}



arb_t sum1,sum2,sum3,sum1a,sum2a,sum3a,sum1b,sum2b,sum3b;

// compute 1/zeta'(s) for s=1/2+i gamma
// checks that 1/2+i gamma is a zero of zeta
void do_rho(arb_t gamma, uint64_t t0, int64_t prec)
{
  //printf("Here I would do something with the zero at ");
  //arb_printd(gamma,20);
  //printf("\n");
  static acb_struct *zeta;
  static acb_t s;
  static arb_t arho,ares,tmp1,tmp2,tmp3;
  static bool init=false;
  static bool under=true;
  if(!init)
    {
      init=true;
      zeta=(acb_struct *)malloc(sizeof(acb_t)*2);
      acb_init(zeta);acb_init(zeta+1);acb_init(s);
      arb_set_d(acb_realref(s),0.5);
      arb_init(arho);arb_init(ares);arb_init(tmp1);arb_init(tmp2);arb_init(tmp3);
    }

  if(under)
    {
      arb_sub_ui(tmp1,gamma,t0/2,prec);
      if(arb_is_positive(tmp1))
	under=false;
    }

  arb_set(acb_imagref(s),gamma);

  acb_dirichlet_zeta_jet(zeta,s,0,2,100);

  if((!arb_contains_zero(acb_realref(zeta)))||(!arb_contains_zero(acb_imagref(zeta))))
    {
      printf("Error, no zero of zeta at ");
      acb_printd(s,40);
      printf("\n");
      exit(0);
    }

  acb_abs(tmp1,zeta+1,prec);
  arb_inv(ares,tmp1,prec); // |1/zeta'|=abs(residue)
  acb_abs(tmp1,s,prec);
  arb_inv(arho,tmp1,prec); // |1/rho|

  arb_mul(tmp1,ares,arho,prec); // |res|/|rho|
  arb_add(sum1,sum1,tmp1,prec);
  if(under)
    arb_add(sum1a,sum1a,tmp1,prec);
  else
    {
      et_gam(tmp2,gamma,t0/2,prec);
      arb_mul(tmp3,tmp2,tmp1,prec);
      arb_add(sum1b,sum1b,tmp3,prec);
    }

  arb_mul(tmp1,tmp1,arho,prec); // |res|/|rho|^2
  arb_add(sum2,sum2,tmp1,prec);
  if(under)
    arb_add(sum2a,sum2a,tmp1,prec);
  else
    {
      arb_mul(tmp3,tmp2,tmp1,prec);
      arb_add(sum2b,sum2b,tmp3,prec);
    }

  arb_mul(tmp1,tmp1,arho,prec); // |res|/|rho|^3
  arb_add(sum3,sum3,tmp1,prec);
  if(under)
    arb_add(sum3a,sum3a,tmp1,prec);
  else
    {
      arb_mul(tmp3,tmp2,tmp1,prec);
      arb_add(sum3b,sum3b,tmp3,prec);
    }

  /*
  acb_inv(res,zeta+1,prec); // 1/zeta'(s) = residue at 1/2+i gamma

  arb_dump_file(outfile,gamma);
  fprintf(outfile,"\n");  
  arb_dump_file(outfile,acb_realref(res));
  fprintf(outfile,"\n");
  arb_dump_file(outfile,acb_imagref(res));
  fprintf(outfile,"\n");
  */
}

int main(int argc, char **argv)
{
  uint64_t i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  fflush(stdout);
  if(argc!=4)
    {
      printf("Fatal error in %s: usage %s <prec> <zeros file> <t0>. Exiting.\n",argv[0],argv[0]);
      exit(0);
    }
  int64_t prec=atol(argv[1]);
  //printf("ARB working precision set to %ld\n",prec);
  FILE *infile=fopen(argv[2],"r");
  if(!infile)
    {
      printf("Fatal error in %s: failed to open zeros file %s for binary input. Exiting.\n",argv[0],argv[2]);
      exit(0);
    }

  int64_t t0=atol(argv[3]);
  
  arb_t gamma,pm1,del_t,t;
  acb_t res;

  arb_init(sum1);arb_init(sum1a);arb_init(sum1b);
  arb_init(sum2);arb_init(sum2a);arb_init(sum2b);
  arb_init(sum3);arb_init(sum3a);arb_init(sum3b);


  arb_init(gamma);arb_init(pm1);arb_init(del_t);arb_init(t);
  acb_init(res);
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
  for(it=0;it<num_its;it++)
    {
      rval=fread(st,sizeof(double),2,infile); // starting/ending t, exact
      rval=fread(zs,sizeof(long int),1,infile); // starting zero number
      if(st[0]==0.0)
	{
	  printf("Iteration %lu empty.\n",it);
	  continue;
	}
      rval=fread(zs+1,sizeof(long int),1,infile); // ending zero number
      //printf("processing zeros %ld to %ld inclusive\n",zs[0]+1,zs[1]);
      n_zeros+=zs[1]-zs[0];
      arb_set_d(gamma,st[0]);
      arb_set(t,gamma);
      //printf("doing t from %f to %f zeros from %ld to %ld\n",st[0],st[1],zs[0],zs[1]);
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
	  do_rho(gamma,t0,prec); // do something with this zero
	}
    }
  printf("sum 1 ");arb_printd(sum1,20);
  printf("\nsum 1a ");arb_printd(sum1a,20);
  printf("\nsum 1b ");arb_printd(sum1b,20);
  printf("\nsum 2 ");arb_printd(sum2,20);
  printf("\nsum 2a ");arb_printd(sum2a,20);
  printf("\nsum 2b ");arb_printd(sum2b,20);
  printf("\nsum 3 ");arb_printd(sum3,20);
  printf("\nsum 3a ");arb_printd(sum3a,20);
  printf("\nsum 3b ");arb_printd(sum3b,20);
  printf("\nWe processed %ld zeros.\n",n_zeros);
  arb_dump_file(stdout,sum1);printf("\n");
  arb_dump_file(stdout,sum1a);printf("\n");
  arb_dump_file(stdout,sum1b);printf("\n");
  arb_dump_file(stdout,sum2);printf("\n");
  arb_dump_file(stdout,sum2a);printf("\n");
  arb_dump_file(stdout,sum2b);printf("\n");
  arb_dump_file(stdout,sum3);printf("\n");
  arb_dump_file(stdout,sum3a);printf("\n");
  arb_dump_file(stdout,sum3b);printf("\n");

  
  return 0;
}


