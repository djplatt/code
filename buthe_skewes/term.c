#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "acb.h"

#define OP_ACC (101)
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
void next_rho(arb_t del_t, FILE *infile, int64_t prec)
{
  in_bytes(del_t,infile,prec);
}

// error from log 2 term
void e1(arb_t res)
{
}

// error from int du/u^2... term
void e2(arb_t res)
{
}

// error from pi(x^1/n)/n n>=2 terms
void e4(arb_t res)
{
}

// error from sum zeros T_0<gam<T1
void e5(arb_t res)
{
}

// error from sum zeros T_1<gam
void e6(arb_t res)
{
}



// Khat
// expects c^2-x^2<0
void Khat(arb_t res, arb_t x, arb_t csinhc, arb_t c2, int64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
    }
  arb_mul(tmp2,x,x,prec);
  arb_sub(tmp1,c2,tmp2,prec);
  if(!arb_is_positive(tmp1))
    {
      fprintf(stderr,"Error in Khat. x^2-c^2 < 0. Exiting.\n");
      exit(0);
    }
  arb_sqrt(tmp2,tmp1,prec);
  arb_sinh(tmp3,tmp2,prec);
  arb_div(tmp1,tmp3,tmp2,prec);
  arb_mul(res,tmp1,csinhc,prec);
  return;
}
  

// Lemma 3 with a=1, k=1
// does both p=1/2+i gam and 1/2-i gam
void do_term(arb_t res, arb_t gam, arb_t om, arb_t c, arb_t eps, int64_t prec)
{
  static arb_t tmp1,tmp2,tmp3,tmp4,csinhc,c2,err;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp4);
      arb_init(csinhc);
      arb_sinh(tmp1,c,prec);
      arb_div(csinhc,c,tmp1,prec);
      arb_init(c2);
      arb_mul(c2,c,c,prec);
      arb_init(err);
      arb_mul_2exp_si(tmp1,eps,-1);
      arb_exp(tmp2,tmp1,prec); // exp(eps/2)
      arb_sub(tmp1,om,eps,prec);
      arb_div(err,tmp2,tmp1,prec);
    }

  arb_mul(tmp1,gam,eps,prec);
  Khat(tmp2,tmp1,csinhc,c2,prec);
  //printf("Khat returned ");arb_printd(tmp2,50);printf("\n");
  arb_mul(tmp1,om,gam,prec);
  arb_sin_cos(tmp3,tmp4,tmp1,prec);
  arb_mul(tmp1,gam,tmp3,prec);
  arb_mul_2exp_si(tmp1,tmp1,1);
  arb_add(tmp3,tmp1,tmp4,prec);
  arb_mul(tmp1,tmp3,tmp2,prec);
  arb_set_d(tmp2,0.25);
  arb_mul(tmp3,gam,gam,prec);
  arb_add(tmp4,tmp2,tmp3,prec);
  arb_div(res,tmp1,tmp4,prec);
  arb_div(tmp1,err,gam,prec);
  arb_div(tmp2,tmp1,gam,prec);
  arb_add_error(res,tmp2);
  
  return;
}

// omega = <om num>/<om den>
// c/eps=T1 limit to which we know RH
// ac/eps=T0 limit to which we will sum zeros
int main(int argc, char **argv)
{
  printf("Command Line:- ");
  for(uint64_t i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");

  if(argc!=8)
    {
      printf("usage:- %s <zeros file> <om num> <om den> <T0> <T1> <eps> < prec>.\n",argv[0]);
      return 0;
    }

  int64_t prec=atol(argv[7]);
  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    return 0;

  double deps=atof(argv[6]);
  
  
  arb_t res,gam,om,c,eps,a;
  arb_init(res);
  arb_init(gam);
  arb_init(om);
  arb_init(c);
  arb_init(eps);
  arb_init(a);
  arb_set_ui(om,atol(argv[2]));
  arb_div_ui(om,om,atol(argv[3]),prec);
  
  arb_set_d(eps,deps);
  arb_set_d(res,atof(argv[5]));
  arb_mul(c,res,eps,prec); // c=eps*T1
  arb_set_d(gam,atof(argv[4])); // T0
  arb_div(a,gam,res,prec); // a= T0/T1
  printf("Setup:-\n");
  printf("   c = ");arb_printd(c,20);printf("\n");
  printf(" eps = ");arb_printd(eps,20);printf("\n");
  printf("   a = ");arb_printd(a,20);printf("\n");
  printf("  om = ");arb_printd(om,20);printf("\n");

  arb_t t,del_t,tot,pm1;
  arb_init(t);
  arb_init(del_t);
  arb_init(tot);
  arb_init(pm1);
  arb_set_ui(del_t,1);
  arb_mul_2exp_si(del_t,del_t,-OP_ACC-1);
  arb_zero(pm1);
  arb_add_error(pm1,del_t);
  
  int64_t num_its,it,z;
  double st[2];
  int64_t zs[2];
  if(fread(&num_its,sizeof(int64_t),1,infile)!=1)
    return 0;
  //printf("Doing %ld iterations.\n",num_its);
  for(it=0;it<num_its;it++)
    {
      if(fread(st,sizeof(double),2,infile)!=2)
	return 0; // starting/ending t, exact
      if(fread(&zs[0],sizeof(int64_t),1,infile)!=1)
	return 0; // starting zero number
      if(st[0]==0.0)
	continue;
      if(fread(&zs[1],sizeof(int64_t),1,infile)!=1)
	return 0; // ending zero number

      arb_set_d(gam,st[0]);
      arb_set(t,gam);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,infile,prec); // distance to next gamma
          if(arb_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  arb_add(t,t,del_t,prec); // exact, used to avoid accumulation of +/- 2^-(OP_ACC+1)
	  arb_add(gam,t,pm1,prec); // not exact
	  do_term(res,gam,om,c,eps,prec);
	  arb_add(tot,tot,res,prec);
	}
    }

  
  
  arb_printd(tot,30);printf("\n");
  fmpz_t aa,bb,ee;
  fmpz_init(aa);
  fmpz_init(bb);
  fmpz_init(ee);
  arb_get_interval_fmpz_2exp(aa,bb,ee,tot);
  fmpz_print(aa);printf(" ");
  fmpz_print(bb);printf(" ");
  fmpz_print(ee);printf("\n");

  return 0;
}
