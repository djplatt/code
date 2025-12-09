#include "stdlib.h"
#include "stdbool.h"
#include "inttypes.h"
#include "flint/arb.h"
#include "flint/acb.h"
#include "flint/acb_dirichlet.h"

#define MAX_NEWTON (10)

void newton(acb_t res, acb_t zp, const acb_t s, int64_t prec)
{
  static bool init=false;
  static acb_struct *zeta;
  if(!init)
    {
      init=true;
      zeta=(acb_struct *)malloc(sizeof(acb_struct)*3);
      for(uint64_t i=0;i<3;i++)
	acb_init(zeta+i);
    }
  //printf("    ");acb_printd(s,20);printf("\n");
  acb_dirichlet_zeta_jet(zeta,s,0,3,prec); // zeta, zeta', zeta''/2
  acb_set(zp,zeta+1);
  acb_div(zeta,zeta+1,zeta+2,prec);
  acb_mul_2exp_si(zeta,zeta,-1);
  acb_sub(res,s,zeta,prec);
}

bool find_zero(acb_t res, acb_t zzp, const acb_t s0, const arb_t tol, int64_t prec)
{
  static bool init=false;
  static acb_t old,new,zp;
  static arb_t tmp,tmp1;
  if(!init)
    {
      init=true;
      acb_init(old);
      acb_init(new);
      arb_init(tmp);
      arb_init(tmp1);
      acb_init(zp);
    }
  acb_set(old,s0);
  int count=0;
  while(true)
    {
      newton(new,zp,old,prec);
      if(!arb_is_finite(acb_realref(new)))
	return false;
      arb_sub_ui(tmp,acb_realref(new),3,prec);
      if(arb_is_positive(tmp))
	return false; // failed to converge
	
      acb_swap(new,old);
      acb_abs(tmp,zp,prec);
      arb_sub(tmp1,tmp,tol,prec);
      if(arb_is_negative(tmp1))
	break;
      
      arb_get_rad_arb(tmp1,tmp);
      arb_sub(tmp1,tmp1,tol,prec);
      if(arb_is_positive(tmp1)) // error is > tolerance
	return false;
      count++;
      if(count>MAX_NEWTON)
	return false;
    }
  acb_set(res,new);
  acb_set(zzp,zp);
  return true;
}

int main(int argc, char**argv)
{

  printf("Command:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=6)
    {
      printf("Usage:- %s <START> <STEP> <ITS> <tol> <PREC>.\n",argv[0]);
      return(0);
    }
  double start=atof(argv[1]);
  int64_t prec=atol(argv[5]);
  arb_t tol;
  arb_init(tol);
  arb_set_d(tol,atof(argv[4]));
  double step=atof(argv[2]);
  int64_t num_its=atol(argv[3]);


  
  acb_t s,rho;
  acb_struct *zeta;
  arb_t tmp;
  zeta=(acb_struct *)malloc(sizeof(acb_struct)*2);
  acb_init(zeta);acb_init(zeta+1);
  acb_init(s);acb_init(rho);
  arb_init(tmp);
  acb_t zp;
  acb_init(zp);

  double t0=start,t1=start+step;

  for(int64_t it=0;it<num_its;it++)
    {
      for(double sig=0.0;sig<1.5;sig+=0.25)
	{
	  arb_set_d(acb_realref(s),sig);
	  for(double t=t0+0.125;t<t1;t+=0.25)
	    {
	      arb_set_d(acb_imagref(s),t);
	      if(find_zero(rho,zp,s,tol,prec))
		{
		  acb_abs(tmp,zp,prec);
		  printf("Zero (sig0 = %f t0 = %f) near ",sig,t);
		  acb_printd(rho,50);
		  printf("|zeta'(rho) = ");arb_printd(tmp,20);
		  printf("\n");
		}
	    }
	}
      t0+=step;
      t1+=step;
    }

  return 0;
}
