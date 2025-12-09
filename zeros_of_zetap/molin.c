#include "stdlib.h"
#include "stdbool.h"
#include "quad.h"
#include "flint/acb_dirichlet.h"

double BOTTOM,TOP,LEFT,RIGHT;

// Return Im zeta''/zeta' (t+20i)
void f2(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_struct *zeta;
  static acb_t s;
  if(!init)
    {
      init=true;
      zeta=(acb_struct *)malloc(sizeof(acb_struct)*3);
      for(uint64_t i=0;i<3;i++)
	acb_init(zeta+i);
      acb_init(s);
    }
  arb_set_d(acb_imagref(s),TOP);
  arb_set(acb_realref(s),t);
  acb_dirichlet_zeta_jet(zeta,s,0,3,prec); // zeta, zeta', zeta''/2
  acb_div(zeta,zeta+2,zeta+1,prec);
  arb_mul_2exp_si(res,acb_imagref(zeta),1);
}

// Return Im zeta''/zeta' (t+10i)
void f4(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_struct *zeta;
  static acb_t s;
  if(!init)
    {
      init=true;
      zeta=(acb_struct *)malloc(sizeof(acb_struct)*3);
      for(uint64_t i=0;i<3;i++)
	acb_init(zeta+i);
      acb_init(s);
    }
  arb_set_d(acb_imagref(s),BOTTOM);
  arb_set(acb_realref(s),t);
  acb_dirichlet_zeta_jet(zeta,s,0,3,prec); // zeta, zeta', zeta''/2
  acb_div(zeta,zeta+2,zeta+1,prec);
  arb_mul_2exp_si(res,acb_imagref(zeta),1);
}



// Return Re zeta''/zeta' (0+It)
void f3(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_struct *zeta;
  static acb_t s;
  if(!init)
    {
      init=true;
      zeta=(acb_struct *)malloc(sizeof(acb_struct)*3);
      for(uint64_t i=0;i<3;i++)
	acb_init(zeta+i);
      acb_init(s);
    }
  arb_set_d(acb_realref(s),LEFT);
  arb_set(acb_imagref(s),t);
  acb_dirichlet_zeta_jet(zeta,s,0,3,prec); // zeta, zeta', zeta''/2
  acb_div(zeta,zeta+2,zeta+1,prec);
  arb_mul_2exp_si(res,acb_realref(zeta),1);
}
// Return Re zeta''/zeta' (4+It)
void f1(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_struct *zeta;
  static acb_t s;
  if(!init)
    {
      init=true;
      zeta=(acb_struct *)malloc(sizeof(acb_struct)*3);
      for(uint64_t i=0;i<3;i++)
	acb_init(zeta+i);
      acb_init(s);
    }
  arb_set_d(acb_realref(s),RIGHT);
  arb_set(acb_imagref(s),t);
  acb_dirichlet_zeta_jet(zeta,s,0,3,prec); // zeta, zeta', zeta''/2
  acb_div(zeta,zeta+2,zeta+1,prec);
  arb_mul_2exp_si(res,acb_realref(zeta),1);
}

int main(int argc, char **argv)
{
  printf("Command:- %s ",argv[0]);
  for(uint64_t i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=8)
    {
      printf("Usage:- %s <LEFT> <RIGHT> <BOTTOM> <STEP> <ITS> <N> <PREC>.\n",argv[0]);
      return(0);
    }
  LEFT=atof(argv[1]);
  RIGHT=atof(argv[2]);
  int64_t prec=atol(argv[7]);
  int64_t n=atol(argv[6]);
  BOTTOM=atof(argv[3]);
  double STEP=atof(argv[4]);
  TOP=BOTTOM+STEP;
  int64_t its=atol(argv[5]);
  
  arb_t maxd,low,high,res3,res1,res2,res4,res;
  arb_t two_pi;
  arb_init(two_pi);
  arb_const_pi(two_pi,prec);
  arb_mul_2exp_si(two_pi,two_pi,1);
  
  arb_init(maxd);arb_init(low);arb_init(high);
  arb_init(res3);arb_init(res1);arb_init(res2);arb_init(res4);
  arb_init(res);
  bool first=true;
  for(uint64_t it=0;it<its;it++)
    {
      if(first)
	{
	  first=false;
	  arb_set_d(low,LEFT);
	  arb_set_d(high,RIGHT);
	  molin_int(res4,n,f4,maxd,low,high,prec);
	}
      else
	arb_neg(res4,res2);


      
      arb_set_d(low,BOTTOM);
      arb_set_d(high,TOP);
      molin_int(res1,n,f1,maxd,low,high,prec);

      arb_set_d(low,LEFT);
      arb_set_d(high,RIGHT);
      molin_int(res2,n,f2,maxd,low,high,prec);
      arb_neg(res2,res2);

      arb_set_d(low,BOTTOM);
      arb_set_d(high,TOP);
      molin_int(res3,n,f3,maxd,low,high,prec);
      arb_neg(res3,res3);


      printf("int Re zeta''/zeta'(%f+it) t=%f..%f = ",RIGHT,BOTTOM,TOP);arb_printd(res1,20);printf("\n");

      printf("int Im zeta''/zeta'(t+%fi) t=%f..%f = ",TOP,RIGHT,LEFT);arb_printd(res2,20);printf("\n");

      printf("int Re zeta''/zeta'(%f+it) t=%f..%f = ",LEFT,TOP,BOTTOM);arb_printd(res3,20);printf("\n");

      printf("int Im zeta''/zeta'(t+%fi) t=%f..%f = ",BOTTOM,LEFT,RIGHT);arb_printd(res4,20);printf("\n");

      arb_set(res,res1);
      arb_add(res,res,res2,prec);
      arb_add(res,res,res3,prec);
      arb_add(res,res,res4,prec);

      printf("I = ");arb_printd(res,20);printf("\n");
  
      arb_div(res,res,two_pi,prec);
      printf("Found ");arb_printd(res,20);printf(" zeros.\n");

      BOTTOM=TOP;
      TOP+=STEP;
    }
  return 0;
}
