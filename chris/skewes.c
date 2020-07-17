/*

File: skewes.c

Created: 3 May 2016

Version: 1.0

Last Modified: 

Dialect: C

Requires: GMP
          MPFR
          MPFI


V 1.0 Initial implementation


Build instructions: gcc -oskewes skewes -O2 -lmpfi

By: DJ Platt
    Bristol University

Copyright 2016.

The author is funded by the UK ESPRC. 

Compute sum over zeros of zeta in interval arithmetic for Chris Smith

*/

#include "stdio.h"
#include "stdlib.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "time.h"
//#include "../includes/mpfi_c.h"

#define mpfi_print_str(a,b) {printf("%s",a);mpfi_out_str(stdout,10,0,b);printf("\n");}

#define OP_ACC (101)
#include "inttypes.h"

mpfi_t pm1; // +/- 1 in the last place

// initialise some variables
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

// read a 13 byte number from file
// structured 8,4,1
// read as if its exact
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


// Chris Smith (York) four term version
inline void do_rho(mpfi_t res1, mpfi_t res2, mpfi_t res3, mpfi_t res4, mpfi_t gamma, mpfi_t omega, mpfi_t two_alpha)
{
  // set up static variables
  static int init=(1==0);
  static mpfi_t g2,gam_om,tmp,tmp1,tmp2,tmp3,tmp4,gsgo,cgo,den,deninc;
  if(!init)
    {
      init=(1==1);
      mpfi_init(g2);
      mpfi_init(gam_om);
      mpfi_init(tmp);
      mpfi_init(tmp1);
      mpfi_init(tmp2);
      mpfi_init(tmp3);
      mpfi_init(tmp4);
      mpfi_init(gsgo);
      mpfi_init(cgo);
      mpfi_init(den);
      mpfi_init(deninc);
    }

  mpfi_sqr(g2,gamma); // gamma^2
  mpfi_mul(gam_om,gamma,omega); // gamma*omega
  mpfi_sin(tmp,gam_om);
  mpfi_mul(gsgo,tmp,gamma); // gamma*sin(gamma*omega)
  mpfi_cos(cgo,gam_om);
  mpfi_mul_ui(tmp,g2,4);
  mpfi_add_ui(tmp1,tmp,1);
  mpfi_set_ui(tmp,1);
  mpfi_div(den,tmp,tmp1); // 1/(1+4gamma^2)

  mpfi_mul_ui(tmp,cgo,4);
  mpfi_mul_ui(tmp1,gsgo,8);
  mpfi_add(tmp2,tmp,tmp1);
  mpfi_mul(res1,tmp2,den); // 1st term
  //mpfi_print_str("First term = ",res1);

  mpfi_div(deninc,den,omega);
  mpfi_mul(den,den,deninc);
  //mpfi_print_str("den now ",den);

  mpfi_mul_ui(tmp,gsgo,32); // 32gamma*sin(gamma*omega)
  mpfi_mul_si(tmp2,g2,-32); // -32gamma^2
  mpfi_add_ui(tmp3,tmp2,8); //(8-32gamma^2)
  mpfi_mul(tmp2,tmp3,cgo); // (8-32gamma^2)*cos(gamma*omega)
  mpfi_add(tmp3,tmp,tmp2);
  mpfi_mul(res2,tmp3,den);
  //mpfi_print_str("Second term = ",res2);


  mpfi_mul(den,den,deninc);

  mpfi_mul_si(tmp,g2,-192);
  mpfi_add_ui(tmp1,tmp,16);
  mpfi_mul(tmp,tmp1,cgo);
  mpfi_mul_si(tmp1,g2,-128);
  mpfi_add_ui(tmp2,tmp1,96);
  mpfi_mul(tmp1,tmp2,gsgo);
  mpfi_add(tmp2,tmp1,tmp);
  mpfi_mul(res3,tmp2,den); // 3rd term

  mpfi_mul(den,den,deninc);

  mpfi_mul_si(tmp1,g2,-1024);
  mpfi_add_ui(tmp2,tmp1,256);
  mpfi_mul(tmp1,tmp2,gsgo); // (256-1024g^2)gsin(gw)
  mpfi_sqr(tmp2,g2);
  mpfi_mul_ui(tmp3,tmp2,512);
  mpfi_mul_si(tmp2,g2,-768);
  mpfi_add(tmp4,tmp2,tmp3);
  mpfi_add_ui(tmp2,tmp4,32);
  mpfi_mul(tmp4,tmp2,cgo);
  mpfi_add(tmp3,tmp4,tmp1);
  mpfi_mul(res4,tmp3,den);


  mpfi_mul(tmp,g2,two_alpha);
  mpfi_exp(tmp2,tmp);
  mpfi_mul(res1,res1,tmp2);
  mpfi_mul(res2,res2,tmp2);
  mpfi_mul(res3,res3,tmp2);
  mpfi_mul(res4,res4,tmp2);
  /*
  mpfi_print_str("gamma = ",gamma);
  mpfi_print_str("res1 = ",res1);
  mpfi_print_str("res2 = ",res2);
  mpfi_print_str("res3 = ",res3);
  mpfi_print_str("res4 = ",res4);
  mpfi_add(tmp,res1,res2);
  mpfi_add(tmp1,res3,res4);
  mpfi_add(tmp2,tmp,tmp1);
  mpfi_print_str("Returning ",tmp2);exit(0);
  */
}

int main(int argc, char **argv)
{
  int i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=6)
    {
      printf("Fatal error in main: usage skewes <prec> <zeros file> <int omega> <omega div> <alpha>. Exiting.\n");
      exit(0);
    }
  mpfr_set_default_prec(atoi(argv[1]));
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }

  mpfi_t gamma,omega,res1,res2,res3,res4,tot1,tot2,tot3,tot4,del_t,t,two_alpha;
  mpfi_init(gamma);mpfi_init(omega);mpfi_init(del_t);mpfi_init(t);
  mpfi_init(res1);mpfi_init(tot1);
  mpfi_init(res2);mpfi_init(tot2);
  mpfi_init(res3);mpfi_init(tot3);
  mpfi_init(res4);mpfi_init(tot4);



  mpfi_set_ui(omega,atol(argv[3]));
  mpfi_div_ui(omega,omega,atol(argv[4]));
  //mpfi_out_str(stdout,10,0,omega);printf("\n");
  
  double alpha_d=atof(argv[5]);//66.0e11;//5444646098003630.0;

  mpfi_print_str("omega=",omega);

  mpfi_init(two_alpha);

  mpfi_set_d(two_alpha,-alpha_d);
  mpfi_print_str("-alpha = ",two_alpha);
  mpfi_mul_2ui(two_alpha,two_alpha,1);
  mpfi_t one;
  mpfi_init(one);
  mpfi_set_ui(one,1);
  mpfi_div(two_alpha,one,two_alpha); // -1/(2*alpha)
  
  init_in_bytes();

  mpfi_set_ui(tot1,0);
  mpfi_set_ui(tot2,0);
  mpfi_set_ui(tot3,0);
  mpfi_set_ui(tot4,0);

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
      mpfi_set_d(gamma,st[0]);
      mpfi_set(t,gamma);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile); // distance to next gamma
          if(mpfi_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  mpfi_add(t,t,del_t); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  mpfi_add(gamma,t,pm1); // not exact

	  // compute Chris's sum at this gamma
	  do_rho(res1,res2,res3,res4,gamma,omega,two_alpha);
	  // add the results to our running totals
	  mpfi_add(tot1,tot1,res1);
	  mpfi_add(tot2,tot2,res2);
	  mpfi_add(tot3,tot3,res3);
	  mpfi_add(tot4,tot4,res4);
	}
    }
  printf("term by term we got\n1     ");
  mpfi_out_str(stdout,10,0,tot1);
  printf("\n");
  printf("2     ");
  mpfi_out_str(stdout,10,0,tot2);
  printf("\n");
  printf("3     ");
  mpfi_out_str(stdout,10,0,tot3);
  printf("\n");
  printf("4     ");
  mpfi_out_str(stdout,10,0,tot4);
  printf("\n");
  mpfi_add(tot1,tot1,tot2);
  mpfi_add(tot1,tot1,tot3);
  mpfi_add(tot1,tot1,tot4);

  printf("skewes on file %s at omega = %f returned ",(double)atol(argv[3])/atol(argv[4]),argv[2]);
  mpfi_out_str(stdout,10,0,tot1);
  printf("\n");
  return(0);
}
