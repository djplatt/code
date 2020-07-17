/*

File: skewes.c

Created: 24 January 2012

Version: 1.0

Last Modified: 

Dialect: C

Requires: GMP
          MPFR
          MPFI

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 1.0 Initial implementation


Build instructions: gcc -oskewes skewes -O2 -lmpfi

By: DJ Platt
    Bristol University

Copyright 2012.

The author is funded by the UK ESPRC. */

#include "stdio.h"
#include "stdlib.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
#include "../includes/mpfi_c.h"

#define OP_ACC (101)
#include "inttypes.h"

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

//
// 
//
/*
mpfi_t do_rho_tmp1,do_rho_tmp2,do_rho_tmp3,do_rho_tmp4;
void init_do_rho()
{
  mpfi_init(do_rho_tmp1);
  mpfi_init(do_rho_tmp2);
  mpfi_init(do_rho_tmp3);
  mpfi_init(do_rho_tmp4);
}

inline void do_rho(mpfi_t res, mpfi_t gamma, mpfi_t omega, mpfi_t two_alpha)
{
  mpfi_sqr(do_rho_tmp1,gamma); // 1=gam^2
  mpfi_mul(do_rho_tmp2,omega,gamma); // 2=gam*om
  mpfi_sin(do_rho_tmp3,do_rho_tmp2); // 3=sin gam*om
  mpfi_cos(do_rho_tmp4,do_rho_tmp2); // 4=cos
  mpfi_mul(do_rho_tmp2,do_rho_tmp3,gamma); // 2=gam*sin
  mpfi_mul_2ui(do_rho_tmp2,do_rho_tmp2,1); // 2=2gam*sin
  mpfi_add(do_rho_tmp3,do_rho_tmp2,do_rho_tmp4); // 3= cos+2gam sin
  //printf("cos+2 gam sin=");mpfi_out_str(stdout,10,0,do_rho_tmp3);printf("\n");
  mpfi_add_d(do_rho_tmp2,do_rho_tmp1,0.25); // 2=gam^2+1/4
  mpfi_div(do_rho_tmp4,do_rho_tmp3,do_rho_tmp2); // 4=(cos+2gam sin)/(gam^2+1/4)
  mpfi_mul(do_rho_tmp2,do_rho_tmp1,two_alpha); // 1=-gam^2/2al
  //mpfi_neg(do_rho_tmp1,do_rho_tmp2); // 2=-gam^2/2a;
  mpfi_exp(do_rho_tmp1,do_rho_tmp2); // 1=exp(-gam^2/2al)
  mpfi_mul(res,do_rho_tmp1,do_rho_tmp4); // res=(cos+2gam sin)*exp(-gam^2/2al)/(gam^2+1/4)
}
*/

mpfi_t g2,g4,om_et,d1,d2,gam_om,cs,ss,n1,n2;
void init_do_rho(mpfi_t omega, double eta)
{
  mpfi_init(g2);
  mpfi_init(g4);
  mpfi_init(om_et);
  mpfi_set_d(g4,eta);
  mpfi_neg(g2,g4);
  mpfi_put(g4,g2);
  mpfi_add(om_et,g4,omega);
  mpfi_t one;mpfi_init(one);mpfi_set_ui(one,1);
  mpfi_div(om_et,one,om_et); // 1/[om-et,om+et]
  mpfi_clear(one);
  mpfi_init(d1);
  mpfi_init(d2);
  mpfi_init(gam_om);
  mpfi_init(ss);
  mpfi_init(cs);
  mpfi_init(n1);
  mpfi_init(n2);
}

// this version includes the second term in expansion of Li(exp(u*rho))
inline void do_rho(mpfi_t res, mpfi_t gamma, mpfi_t omega, mpfi_t two_alpha, double eta)
{
  mpfi_sqr(g2,gamma);
  mpfi_mul(n1,g2,two_alpha);
  mpfi_exp(n2,n1);
  mpfi_mul_2ui(n2,n2,2); // 4exp(-gam^2/2alpha)
  //mpfi_print_str("4exp(-gam^2/2alpha=",n2);
  mpfi_mul_2ui(g2,g2,2); // 4 gam^2
  mpfi_mul(gam_om,gamma,omega);
  mpfi_cos(cs,gam_om); // cos(gamma*omega)
  mpfi_sin(ss,gam_om);
  mpfi_mul_2ui(d1,ss,1);
  mpfi_mul(ss,d1,gamma); // 2*gamma*sin(gamma*omega)
  mpfi_sub_ui(d1,g2,1);
  mpfi_neg(d2,d1);
  mpfi_mul_2ui(d2,d2,1);
  mpfi_mul(d1,d2,om_et); // 2(1-4gam^2)/[om_et]
  mpfi_add(d1,d1,g2);
  mpfi_add_ui(d1,d1,1);
  mpfi_mul(d2,cs,d1); // cos()*(1+4gam^2+2(1-4gam^2)/[om_et])
  //mpfi_print_str("cos term=",d2);
  mpfi_mul_2ui(d1,om_et,2); // 4/[om_et]
  mpfi_add(d1,d1,g2);
  mpfi_add_ui(d1,d1,1);
  mpfi_mul(cs,d1,ss); // 2gam*sin()(1+4gam^2+8/[om_et])
  //mpfi_print_str("sin term=",cs);
  mpfi_add(d1,d2,cs);
  mpfi_add_ui(g2,g2,1);
  mpfi_sqr(cs,g2);
  mpfi_div(n1,n2,cs);
  mpfi_mul(res,n1,d1);

}

int main(int argc, char **argv)
{
  int i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=7)
    {
      printf("Fatal error in main: usage skewes <prec> <zeros file> <int omega> <omega div> <eta> <alpha>. Exiting.\n");
      exit(0);
    }
  mpfr_set_default_prec(atoi(argv[1]));
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }

  mpfi_t gamma,omega,res,tot,del_t,t,two_alpha;
  mpfi_init(gamma);mpfi_init(omega);mpfi_init(res);mpfi_init(tot);mpfi_init(del_t);mpfi_init(t);




  mpfi_set_ui(omega,atol(argv[3]));
  mpfi_div_ui(omega,omega,atol(argv[4]));
  //mpfi_out_str(stdout,10,0,omega);printf("\n");
  
  double alpha_d=atof(argv[6]);//66.0e11;//5444646098003630.0;

  double eta=atof(argv[5]);
  printf("eta=%30.28e\nalpha=%30.28e\n",eta,alpha_d);
  mpfi_print_str("omega=",omega);

  mpfi_init(two_alpha);
  mpfi_sub_d(two_alpha,omega,eta);
  mpfi_print_str("Lower bound=",two_alpha);
  mpfi_add_d(two_alpha,omega,eta);
  mpfi_print_str("Upper bound=",two_alpha);

  mpfi_set_d(two_alpha,-alpha_d);
  mpfi_mul_2ui(two_alpha,two_alpha,1);
  mpfi_t one;
  mpfi_init(one);
  mpfi_set_ui(one,1);
  mpfi_div(two_alpha,one,two_alpha); // -1/(2*alpha)
  
  init_in_bytes();
  init_do_rho(omega,eta);

  mpfi_set_ui(tot,0);

  long int num_its,it,z;
  double st[2];
  long int zs[2];
  fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(&zs[0],sizeof(long int),1,infile);
      if(st[0]==0.0)
	continue;
      fread(&zs[1],sizeof(long int),1,infile);
      //printf("Processing zero %ld to %ld=%ld in total.\nt0=%f\n",zs[0]+1,zs[1],zs[1]-zs[0],st[0]);
      mpfi_set_d(gamma,st[0]);
      mpfi_set(t,gamma);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile);
          if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_add(t,t,del_t); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  mpfi_add(gamma,t,pm1); // not exact
	  //printf("Doing zero at ");mpfi_out_str(stdout,10,0,gamma);printf("\n");
	  do_rho(res,gamma,omega,two_alpha,eta);
	  mpfi_add(tot,tot,res);
	  //printf("total=");mpfi_out_str(stdout,10,0,tot);printf("\n");exit(0);
	}
    }
  printf("skewes on file %s returned ",argv[2]);
  mpfi_out_str(stdout,10,0,tot);
  printf("\n");
  return(0);
}
