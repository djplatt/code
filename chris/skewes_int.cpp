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

// Re (Li(exp(u(1/2 + i gam)))+Li(exp(u(1/2-i gam))))/exp(u/2)
void li_star(arb_t res, const arb_t u, const arb_t gam, uint64_t n, int64_t prec)
{
  static bool init=false;
  static arb_t gamu,tmp1,tmp2;
  static acb_t rho,ctmp1,ctmp2,ctmp3,rhou,rhouk;
  if(!init)
    {
      init=true;
      arb_init(gamu);
      arb_init(tmp1);
      arb_init(tmp2);
      acb_init(rho);
      acb_init(ctmp1);
      acb_init(ctmp2);
      acb_init(ctmp3);
      acb_init(rhou);
      acb_init(rhouk);
    }
  uint64_t kfac=1;

  arb_set_d(acb_realref(rho),0.5);
  arb_set(acb_imagref(rho),gam); // rho =1/2+i gam

  acb_mul_arb(ctmp1,rho,u,prec);
  acb_inv(rhou,ctmp1,prec); // rhou=1/rho u
  acb_set(rhouk,rhou); // 1/(rho u)^1
  acb_set(ctmp1,rhou);
  for(uint64_t k=2;k<=n;k++)
    {
      kfac*=k-1;
      acb_mul(rhouk,rhouk,rhou,prec); // 1/(rho u)^k
      acb_mul_ui(ctmp2,rhouk,kfac,prec); // (k-1)!/(rho u)^k
      acb_add(ctmp1,ctmp1,ctmp2,prec);
    }

  arb_zero(acb_realref(ctmp2));
  arb_set(acb_imagref(ctmp2),gamu);
  acb_exp(ctmp3,ctmp2,prec);

  acb_mul(ctmp2,ctmp1,ctmp3,prec);

  arb_mul(gamu,gam,u,prec);
  arb_pow_ui(tmp1,gamu,n+1,prec);
  arb_inv(tmp2,tmp1,prec);
  arb_mul_ui(tmp1,tmp2,kfac*n,prec);
  arb_add_error(acb_realref(ctmp2),tmp1);
  arb_mul_2exp_si(res,acb_realref(ctmp2),1);
}

double alpha_d;
arb_t omega,eta;

// Gaussian Kernel K(y-omega) as per S. T. D. 2.3
void K(arb_t res, const arb_t y, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,min_al_by_2,sqrt_al_by_2pi;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(min_al_by_2);
      arb_set_d(min_al_by_2,-alpha_d/2.0);
      arb_init(sqrt_al_by_2pi);
      arb_const_pi(tmp1,prec);
      arb_div(tmp2,min_al_by_2,tmp1,prec); // -alpha/2Pi
      printf("-al/2 = ");arb_printd(min_al_by_2,20);printf("\n");
      arb_neg(tmp2,tmp2);
      arb_sqrt(sqrt_al_by_2pi,tmp2,prec);
      printf("sqrt(al/2Pi) = ");arb_printd(sqrt_al_by_2pi,20);printf("\n");
    }
  //printf("K called with u = ");arb_printd(y,20);printf("\n");
  arb_sub(tmp2,y,omega,prec);
  //printf("Doing exp with y = ");arb_printd(tmp2,20);printf("\n");
  arb_mul(tmp1,tmp2,tmp2,prec); // y^2
  arb_mul(tmp2,tmp1,min_al_by_2,prec); // -y^2 alpha/2
  arb_exp(tmp1,tmp2,prec);
  arb_mul(res,tmp1,sqrt_al_by_2pi,prec);
  //printf("K returning ");arb_printd(res,20);printf("\n");
}

void one(arb_t res, const arb_t u, int64_t prec)
{
  arb_one(res);
}

arb_t gam;
void li_bit(arb_t res, const arb_t u, int64_t prec)
{
  //printf("Li bit called with u = ");arb_printd(u,20);printf("\n");
  li_star(res,u,gam,10,prec);
  arb_mul(res,res,u,prec);
  //printf("Li bit returning ");arb_printd(res,20);printf("\n");
}


// 
inline void do_rho(arb_t res, arb_t lo, arb_t hi, arb_t maxd, uint64_t n, int64_t prec)
{
  molin_int(res,n,K,li_bit,maxd,lo,hi,prec);
}

int main(int argc, char **argv)
{
  int i;
  int64_t prec=200;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=6)
    {
      printf("Fatal error in main: usage skewes <zeros file> <int omega> <omega div> <alpha> <eta>. Exiting.\n");
      exit(0);
    }
  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }

  arb_t res1,tot1,del_t,t;
  arb_init(gam);arb_init(omega);arb_init(del_t);arb_init(t);
  arb_init(res1);arb_init(tot1);



  arb_set_ui(omega,atol(argv[2]));
  arb_div_ui(omega,omega,atol(argv[3]),prec);
  printf("omega = ");arb_printd(omega,20);printf("\n");
  
  alpha_d=atof(argv[4]);//66.0e11;//5444646098003630.0;
  double eta_d=atof(argv[5]);
  arb_t eta;
  arb_init(eta);
  arb_set_d(eta,eta_d);
  printf("eta = ");arb_printd(eta,20);printf("\n");


  arb_t lo,hi,maxd;
  arb_init(lo);arb_init(hi);arb_init(maxd);
  arb_set_d(maxd,alpha_d*eta_d*eta_d*2);
  arb_exp(maxd,maxd,prec);
  arb_sub(lo,omega,eta,prec);
  arb_add(hi,omega,eta,prec);


  
  init_in_bytes();

  arb_zero(tot1);

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
      arb_set_d(gam,st[0]);
      arb_set(t,gam);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile,prec); // distance to next gam
          if(arb_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  arb_add(t,t,del_t,prec); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  arb_add(gam,t,pm1,prec); // not exact

	  // compute Chris's sum at this gam
	  do_rho(res1,lo,hi,maxd,400,prec);
	  printf("rho on ");arb_printd(gam,20);printf(" returned ");arb_printd(res1,20);printf("\n");
	  // add the results to our running totals
	  arb_add(tot1,tot1,res1,prec);
	}
    }
  printf("term by term we got\n1     ");
  arb_printd(tot1,20);
  printf("\n");
  printf("skewes on file %s at omega = %f returned ",argv[1],(double)atol(argv[2])/atol(argv[3]));
  arb_printd(tot1,prec);
  printf("\n");
  return(0);
}
