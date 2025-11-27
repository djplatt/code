#include "stdlib.h"
#include "stdbool.h"
#include "inttypes.h"
#include "flint/mag.h"
#include "flint/acb.h"
#include "flint/acb_dirichlet.h"
#include "flint/acb_calc.h"


void f(acb_t res, const acb_t s, int64_t prec)
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
  //acb_inv(res,s,prec);return;
  acb_dirichlet_zeta_jet(zeta,s,0,3,prec); // zeta, zeta', zeta''/2
  acb_div(res,zeta+2,zeta+1,prec);
  acb_mul_2exp_si(res,res,1);
  //printf("f called with ");acb_printd(s,20);printf(" returning ");acb_printd(res,20);printf("\n");
}

int func(acb_ptr out, const acb_t inp, void *param, slong order, slong prec)
{
  if (order > 1) flint_abort();

  
  f(out,inp,prec);
  /*
  if(order==1)
    {
      printf("calling order 1 with s = ");acb_printd(inp,20);printf("\n");
      printf("returning ");acb_printd(out,20);printf("\n");
    }
  */
  return 0;
}
  



int main(int argc, char **argv)
{
  if(argc!=6)
    {
      printf("Usage:- %s <TOP> <BOTTOM> <LEFT> <RIGHT> <prec in bits>.\n",argv[0]);
      return 0;
    }
  int64_t TOP=atol(argv[1]);
  int64_t BOTTOM=atol(argv[2]);
  int64_t LEFT=atol(argv[3]);
  int64_t RIGHT=atol(argv[4]);
  int64_t prec=atol(argv[5]);

  acb_t a,b,res,res1,res2,res3,res4;
  acb_init(a);acb_init(b);
  acb_init(res);acb_init(res1);acb_init(res2);acb_init(res3);acb_init(res4);

  arb_t two_pi;
  arb_init(two_pi);
  arb_const_pi(two_pi,prec);
  arb_mul_2exp_si(two_pi,two_pi,1);

  mag_t tol;
  mag_init(tol);
  mag_set_ui_2exp_si(tol,1,-prec);

  acb_calc_integrate_opt_t options;
  acb_calc_integrate_opt_init(options);
  options->verbose=1;
  options->eval_limit=1000000;
  
  arb_set_ui(acb_realref(a),RIGHT);
  arb_set_ui(acb_realref(b),RIGHT);
  arb_set_si(acb_imagref(a),BOTTOM);
  arb_set_ui(acb_imagref(b),TOP);  
  acb_calc_integrate(res1,func,NULL,a,b,prec,tol,options,prec); 

  printf("I1=");acb_printd(res1,20);printf("\n");
  
  arb_set_si(acb_realref(a),RIGHT);
  arb_set_si(acb_imagref(a),TOP);
  arb_set_si(acb_realref(b),LEFT);
  arb_set_si(acb_imagref(b),TOP);
  acb_calc_integrate(res2,func,NULL,a,b,prec,tol,options,prec);

  printf("I2=");acb_printd(res2,30);printf("\n");

  
  arb_set_si(acb_realref(a),LEFT);
  arb_set_si(acb_imagref(a),TOP);
  arb_set_si(acb_realref(b),LEFT);
  arb_set_si(acb_imagref(b),BOTTOM);
  acb_calc_integrate(res3,func,NULL,a,b,prec,tol,options,prec);


  printf("I3=");acb_printd(res3,30);printf("\n");

  arb_set_si(acb_realref(a),LEFT);
  arb_set_si(acb_imagref(a),BOTTOM);
  arb_set_si(acb_realref(b),RIGHT);
  arb_set_si(acb_imagref(b),BOTTOM);
  acb_calc_integrate(res4,func,NULL,a,b,prec,tol,options,prec);

  printf("I4=");acb_printd(res4,30);printf("\n");

  acb_add(res,res1,res2,prec);
  acb_add(res,res,res3,prec);
  acb_add(res,res,res4,prec);

  
  acb_div_arb(res,res,two_pi,prec);
  acb_mul_onei(res,res);
  acb_neg(res,res);

  printf("Zeros from %lu to %lu = ",TOP,BOTTOM);acb_printd(res,20);printf("\n");
  
  return 0;
}
