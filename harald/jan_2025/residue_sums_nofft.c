#include "flint/acb_dirichlet.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define OP_ACC (101)
#include "inttypes.h"

#define ARB_PREC (140) // enough for 102 bits after dp of 10^10 size thing

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


void varphi(acb_t res, acb_t x, int64_t prec)
{
  static bool init=false;
  static arb_t pi2;
  static acb_t ctmp1,ctmp2;
  
  if(!init)
    {
      init=true;
      arb_init(pi2);
      arb_const_pi(pi2,prec);
      arb_mul_2exp_si(pi2,pi2,-1); // pi/2
      acb_init(ctmp1);acb_init(ctmp2);
    }

  acb_mul_arb(ctmp1,x,pi2,prec);
  acb_coth(ctmp2,ctmp1,prec); // coth(x pi/2)
  acb_mul(res,ctmp2,ctmp1,prec); // x pi/2 cot(x pi/2)
}

arb_t max_zetap,min_zetap,max_gamma,min_gamma;

void max_min_zetap(arb_t zetap, arb_t gamma, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(max_zetap); // init to 0
      arb_init(min_zetap);
      arb_set_ui(min_zetap,UWORD_MAX); // very large
      arb_init(max_gamma);
      arb_init(min_gamma);
    }
  //printf("zetap = ");arb_printd(zetap,20);printf("\n");
  arb_sub(tmp,max_zetap,zetap,prec);
  if(arb_is_negative(tmp))
    {
      arb_set(max_zetap,zetap);
      arb_set(max_gamma,gamma);
      return;
    }
  else
    {
      if(!arb_is_positive(tmp))
	{
	  printf("Warning. zetap = ");arb_printd(zetap,20);
	  printf(" and max_zetap = ");arb_printd(max_zetap,20);
	  printf(" overlap at gamma = ");arb_printd(gamma,20);printf("\n");
	}
    }
  arb_sub(tmp,zetap,min_zetap,prec);
  if(arb_is_negative(tmp))
    {
      arb_set(min_zetap,zetap);
      arb_set(min_gamma,gamma);
    }
  else
    {
      if(!arb_is_positive(tmp))
	{
	  printf("Warning. zetap = ");arb_printd(zetap,20);
	  printf(" and min_zetap = ");arb_printd(min_zetap,20);
	  printf(" overlap at gamma = ");arb_printd(gamma,20);printf("\n");
	}
    }
}  


arb_t sum1,sum2,sum3,sum4,sum5;//,sum1a,sum2a,sum3a,sum1b,sum2b,sum3b;

// compute 1/zeta'(s) for s=1/2+i gamma
// checks that 1/2+i gamma is a zero of zeta
void do_rho(arb_t gamma, arb_t t0, int64_t prec)
{
  static acb_struct *zeta;
  static acb_t s,ctmp1,z;
  static arb_t tmp1,t1,t2,varp,inv_zetap,den,inv_s_abs;
  static bool init=false;
  if(!init)
    {
      init=true;
      zeta=(acb_struct *)malloc(sizeof(acb_t)*2);
      acb_init(zeta);acb_init(zeta+1);acb_init(s);acb_init(ctmp1);acb_init(z);
      arb_set_d(acb_realref(s),0.5);
      arb_init(tmp1);arb_init(varp);
      arb_init(inv_s_abs);arb_init(inv_zetap);arb_init(den);
      arb_init(t1);
      arb_init(t2);
      arb_mul_2exp_si(t1,t0,-1);
      arb_div_ui(t2,t0,10,prec);
    }

  prec=prec/2;
  
  arb_set(acb_imagref(s),gamma);
  acb_abs(tmp1,s,prec);
  arb_inv(inv_s_abs,tmp1,prec); // |1/rho|
  acb_dirichlet_zeta_jet(zeta,s,0,2,prec);
  //printf("zeros at ");arb_printd(gamma,20);
  // acb_abs(tmp1,zeta,prec);
  // printf("\n|zeta(t)| = ");arb_printd(tmp1,20);
  //printf("\n|zeta'(t)| = ");arb_printd(tmp1,20);printf("\n");

  if((!arb_contains_zero(acb_realref(zeta)))||(!arb_contains_zero(acb_imagref(zeta))))
    {
      printf("Error, no zero of zeta at ");
      acb_printd(s,40);
      printf("\n");
      exit(0);
    }
  acb_abs(tmp1,zeta+1,prec);
  max_min_zetap(tmp1,gamma,prec);
  arb_inv(inv_zetap,tmp1,prec); // 1/|zeta'|
  arb_add(sum4,sum4,inv_zetap,prec);
  arb_mul(den,inv_zetap,inv_s_abs,prec); // 1/|zeta'(rho) rho|
  arb_add(sum5,sum5,den,prec);
  
  acb_sub_ui(ctmp1,s,1,prec); // rho-1
  acb_div_onei(ctmp1,ctmp1); // (rho-1)/i
  acb_div_arb(z,ctmp1,t0,prec); // (rho-1)/(t0 i)
  
  // first with t0 into sum1
  varphi(ctmp1,z,prec);
  acb_abs(varp,ctmp1,prec); // |varphi(z)|
  //printf("|varphi| = ");arb_printd(varp,20);printf("\n");
  arb_mul(tmp1,varp,den,prec); // 1/|zeta' rho|
  arb_add(sum1,sum1,tmp1,prec);

  arb_sub(tmp1,gamma,t1,prec);
  if(arb_is_positive(tmp1))
    return;

  // now with t1 into sum2
  acb_mul_2exp_si(z,z,1); // (rho-1)/(t0/2 i)
  varphi(ctmp1,z,prec);
  acb_abs(varp,ctmp1,prec); // |varphi(z)|
  arb_mul(tmp1,varp,den,prec); // 1/|zeta' rho|
  arb_add(sum2,sum2,tmp1,prec);

  arb_sub(tmp1,gamma,t2,prec);
  if(arb_is_positive(tmp1))
    return;
  acb_mul_ui(z,z,5,prec); // (rho-1)/(t0/10 i)
  varphi(ctmp1,z,prec);
  acb_abs(varp,ctmp1,prec); // |varphi(z)|
  arb_mul(tmp1,varp,den,prec); // 1/|zeta' rho|
  arb_add(sum3,sum3,tmp1,prec);
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

  arb_init(sum1);//arb_init(sum1a);arb_init(sum1b);
  arb_init(sum2);//arb_init(sum2a);arb_init(sum2b);
  arb_init(sum3);//arb_init(sum3a);arb_init(sum3b);
  arb_init(sum4);// sum 1/zeta'(rho)
  arb_init(sum5);// sum 1/|rho zeta'(rho)|
  
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
  printf("sum (t0=%f) ",(double) t0);arb_printd(sum1,20);
  //printf("\nsum 1a ");arb_printd(sum1a,20);
  //printf("\nsum 1b ");arb_printd(sum1b,20);
  printf("\nsum (t0=%f) ",(double) t0 / 2.0);arb_printd(sum2,20);
  //printf("\nsum 2a ");arb_printd(sum2a,20);
  //printf("\nsum 2b ");arb_printd(sum2b,20);
  printf("\nsum (t0=%f) ",(double) t0 / 10.0);arb_printd(sum3,20);
  printf("\nsum 1/|zeta'| ");arb_printd(sum4,20);
  printf("\nsum 1/|rho zeta'| ");arb_printd(sum5,20);
  printf("\nMax zeta' = ");arb_printd(max_zetap,20);
  printf(" seen at ");arb_printd(max_gamma,20);printf("\n");
  printf("Min zeta' = ");arb_printd(min_zetap,20);
  printf(" seen at ");arb_printd(min_gamma,20);printf("\n");
  printf("\n");

  //printf("\nsum 3a ");arb_printd(sum3a,20);
  //printf("\nsum 3b ");arb_printd(sum3b,20);
  printf("\nWe processed %ld zeros.\n",n_zeros);
  arb_dump_file(stdout,sum1);printf("\n");
  //arb_dump_file(stdout,sum1a);printf("\n");
  //arb_dump_file(stdout,sum1b);printf("\n");
  arb_dump_file(stdout,sum2);printf("\n");
  //arb_dump_file(stdout,sum2a);printf("\n");
  //arb_dump_file(stdout,sum2b);printf("\n");
  arb_dump_file(stdout,sum3);printf("\n");
  arb_dump_file(stdout,sum4);printf("\n");
  arb_dump_file(stdout,sum5);printf("\n");
  //arb_dump_file(stdout,sum3a);printf("\n");
  //arb_dump_file(stdout,sum3b);printf("\n");

  
  return 0;
}


