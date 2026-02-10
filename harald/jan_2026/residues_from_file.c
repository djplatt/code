#include "flint/acb_dirichlet.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define OP_ACC (101)
#include "inttypes.h"

#define ARB_PREC (140) // enough for 102 bits after dp of 10^10 size thing



// Pi/(2T) coth (Pi s/(2T))
void W_t(arb_t res, acb_t s, arb_t delta, int64_t prec)
{
  static bool init=false;
  static acb_t ctmp1,ctmp2;
  
  if(!init)
    {
      init=true;
      acb_init(ctmp1);acb_init(ctmp2);
    }

  acb_mul_arb(ctmp2,s,delta,prec); // delta s
  acb_coth(ctmp1,ctmp2,prec); // coth(delta s)
  acb_abs(res,ctmp1,prec); // |coth(delta s)|
  arb_mul(res,res,delta,prec); // |delta blah|
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


arb_t sum1,sum2,sum3,sum4,sum5,sum6;//,sum1a,sum2a,sum3a,sum1b,sum2b,sum3b;

void do_sum4(arb_t zeta_bit, acb_t rho, int64_t prec)
{
  
  static bool init=false;
  static arb_t abs_rho,tmp1;
  if(!init)
    {
      init=true;
      arb_init(abs_rho);
      arb_init(tmp1);
    }
  
  
  acb_abs(abs_rho,rho,prec);
  arb_div(tmp1,zeta_bit,abs_rho,prec);
  arb_add(sum4,sum4,tmp1,prec);
  return;

  /*
  arb_log(tmp1,acb_imagref(rho),prec);
  arb_neg(tmp1,tmp1);
  arb_mul_2exp_si(tmp1,tmp1,-2);
  arb_exp(T_1_4,tmp1,prec); // gamma^(-1/4)

  arb_mul(tmp1,sum4,T_1_4,prec); // value for C
  arb_sub(tmp2,tmp1,best_C,prec);
  if(!arb_is_negative(tmp2))
    {
      arb_set(best_C,tmp1);
      arb_set(best_C_gamma,acb_imagref(rho));
    }
  */
}

// compute W_t(rho) zeta(rho/2)/zeta'(s) for s=1/2+i gamma, t=t0, t0/10, t0/100
void do_rho(arb_t gamma, acb_t residue, arb_t t0, int64_t prec)
{
  static acb_t zeta;
  static arb_t pi2,delta0,delta1,delta2;
  static acb_t rho,ctmp1,z,rho_2;
  static arb_t tmp1,t1,t2,varp,inv_zetap,zeta_abs,inv_s_abs,zeta_bit;
  static bool init=false;
  if(!init)
    {
      init=true;
      acb_init(zeta);acb_init(rho);acb_init(ctmp1);acb_init(z);
      arb_set_d(acb_realref(rho),0.5);
      arb_init(tmp1);arb_init(varp);
      arb_init(inv_s_abs);arb_init(inv_zetap);arb_init(zeta_abs);
      arb_init(t1);
      arb_init(t2);
      arb_init(zeta_bit);
      arb_div_ui(t1,t0,10,prec);
      arb_div_ui(t2,t0,100,prec);
      acb_init(rho_2);
      arb_set_d(acb_realref(rho_2),0.25);
      arb_init(pi2);
      arb_const_pi(pi2,prec);
      arb_mul_2exp_si(pi2,pi2,-1); // pi/2
      arb_init(delta0);arb_init(delta1);arb_init(delta2);
      arb_div(delta0,pi2,t0,prec);
      arb_div(delta1,pi2,t1,prec);
      arb_div(delta2,pi2,t2,prec);
    }

  arb_set(acb_imagref(rho),gamma); // rho

  arb_mul_2exp_si(acb_imagref(rho_2),gamma,-1); // rho/2
  acb_dirichlet_zeta(zeta,rho_2,prec/5); // zeta(rho/2)
  acb_abs(zeta_abs,zeta,prec); // |zeta(rho/2)|
  
  acb_abs(inv_zetap,residue,prec); // 1/|zeta'|
  //max_min_zetap(inv_zetap,gamma,prec);

  arb_mul(zeta_bit,inv_zetap,zeta_abs,prec); // |zeta(rho/2)/zeta'(rho)|
  do_sum4(zeta_bit,rho,prec);
  W_t(varp,rho,delta0,prec); // |W_t_0(rho)|
  //printf("|W_t| = ");arb_printd(varp,20);printf("\n");
  arb_mul(tmp1,varp,zeta_bit,prec); // W(t)|zeta(rho/2)/zeta'(rho)|
  arb_add(sum1,sum1,tmp1,prec);

  arb_sub(tmp1,gamma,t1,prec);
  if(arb_is_positive(tmp1))
    return;


  
  // now with t1 into sum2
  W_t(varp,rho,delta1,prec); // |W_t_0(rho)|
  //printf("|W_t| = ");arb_printd(varp,20);printf("\n");
  arb_mul(tmp1,varp,zeta_bit,prec); // W(t)|zeta(rho/2)/zeta'(rho)|
  arb_add(sum2,sum2,tmp1,prec);

  
  arb_sub(tmp1,gamma,t2,prec);
  if(arb_is_positive(tmp1))
    return;
  // now with t2 into sum2
  W_t(varp,rho,delta2,prec); // |W_t_0(rho)|
  //printf("|W_t| = ");arb_printd(varp,20);printf("\n");
  arb_mul(tmp1,varp,zeta_bit,prec); // W(t)|zeta(rho/2)/zeta'(rho)|
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
      printf("Fatal error in %s: usage %s <residues file> <t0>. Exiting.\n",argv[0],argv[0]);
      exit(0);
    }
  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Fatal error in %s: failed to open residues file %s for input. Exiting.\n",argv[0],argv[1]);
      exit(0);
    }

  double t0=atof(argv[2]);
  if(t0<=0)
    {
      printf("Need t0 > 0. Exiting.\n");
      exit(0);
    }
  arb_t a_t0;
  arb_init(a_t0);
  arb_set_d(a_t0,t0);
  
  arb_t gamma,tmp,gamma0;
  acb_t residue;

  arb_init(sum1);//arb_init(sum1a);arb_init(sum1b);
  arb_init(sum2);//arb_init(sum2a);arb_init(sum2b);
  arb_init(sum3);//arb_init(sum3a);arb_init(sum3b);
  arb_init(sum4);// sum 1/zeta'(rho)
  arb_init(sum5);// sum 1/|rho zeta'(rho)|
  arb_init(sum6);

  arb_init(gamma);arb_init(tmp);arb_init(gamma0);
  acb_init(residue);
  int zero=0;
  while(arb_load_file(gamma,infile)==0)
    {
      arb_sub(tmp,gamma,a_t0,prec);
      if(arb_is_positive(tmp))
	break;
      if(arb_load_file(acb_realref(residue),infile)!=0)
	{
	  printf("Error reading real part of residue. Exiting.\n");
	  return 0;
	}
      if(arb_load_file(acb_imagref(residue),infile)!=0)
	{
	  printf("Error reading imaginary part of residue. Exiting.\n");
	  return 0;
	}
      if(zero==0)
	arb_set(gamma0,gamma);
      zero++;
      do_rho(gamma,residue,a_t0,prec);
      /*
      if(zero==1000)
	break;
      
      if((zero%1000)==0)
	{
	  printf("Sum |zeta(rho/2)/(rho zeta'(rho))| for gamma in [");
	  arb_printd(gamma0,20);
	  printf(",");
	  arb_printd(gamma,20);
	  printf("] is ");
	  arb_printd(sum4,20);printf("\n");
	  arb_dump_file(stdout,gamma);printf("\n");
	  arb_dump_file(stdout,sum4);printf("\n");
	  arb_zero(sum4);
	  arb_set(gamma0,gamma);
	}
      */  
    }
      



  printf("\nWe processed %d zeros.\n",zero);
  printf("sum (t0=%f) ",(double) t0);arb_printd(sum1,20);
  printf("\nsum (t0=%f) ",(double) t0 / 10.0);arb_printd(sum2,20);
  printf("\nsum (t0=%f) ",(double) t0 / 100.0);arb_printd(sum3,20);
  printf("\nSum |zeta(rho/2)/(rho zeta'(rho))| for gamma in [");
  arb_printd(gamma0,20);
  printf(",");
  arb_printd(gamma,20);
  printf("] is ");
  arb_printd(sum4,20);printf("\n");
  arb_dump_file(stdout,sum1);printf("\n");
  arb_dump_file(stdout,sum2);printf("\n");
  arb_dump_file(stdout,sum3);printf("\n");
  arb_dump_file(stdout,sum4);printf("\n");
  arb_dump_file(stdout,gamma0);printf("\n");
  arb_dump_file(stdout,gamma);printf("\n");

  
  return 0;
}


