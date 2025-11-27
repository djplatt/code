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
  //printf("     In f with s = ");acb_printd(s,20);printf("\n");
  //acb_inv(res,s,prec);return;

  
  acb_dirichlet_zeta_jet(zeta,s,0,3,prec); // zeta/0!, zeta'/1!, zeta''/2!
  acb_div(res,zeta+2,zeta+1,prec);
  acb_mul_2exp_si(res,res,1);
  //printf(" s = ");acb_printd(s,10);printf(" f returning ");acb_printd(res,20);printf("\n");
  /*
  if(!acb_is_finite(res))
    {
      printf("f called with ");acb_printd(s,20);
      printf("\nreturned ");acb_printd(res,20);
      printf("\nzeta' = ");acb_printd(zeta+1,20);
      printf("\nzeta'' = ");acb_printd(zeta+2,20);
      printf("\n");
      flint_abort();
    }
  */
}




// int from a to b in steps

void my_int(acb_t res, const acb_t a, const acb_t b, const int64_t steps, const int64_t prec)
{
  static bool init=false;
  static acb_t tmp1,tmp2,delta,mid,d2;
  if(!init)
    {
      init=false;
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(delta);
      acb_init(mid);
      acb_init(d2);
    }

  acb_sub(tmp1,b,a,prec);
  acb_div_ui(delta,tmp1,steps,prec);
  if(!acb_is_exact(delta))
    {
      printf("(b-a)/steps must be exact. Exiting.\n");
      exit(0);
    }

  
  acb_mul_2exp_si(d2,delta,-1);
  acb_add(mid,a,d2,prec);
  arb_add_error(acb_realref(mid),acb_realref(d2));
  arb_add_error(acb_imagref(mid),acb_imagref(d2));
  acb_zero(res);
  printf("delta = ");acb_printd(delta,30);printf("\n");
  for(int64_t step=0;step<steps;step++)
    {
      f(tmp1,mid,prec); // f[mid +/- delta/2]
      acb_add(res,res,tmp1,prec);
      acb_add(mid,mid,delta,prec);
      if(!acb_is_finite(res))
	return;
    }
  acb_mul(res,res,delta,prec);
  
  //printf("res = ");acb_printd(res,20);printf("\n");
}




int main(int argc, char **argv)
{
  if(argc!=7)
    {
      printf("Usage:- %s <TOP> <BOTTOM> <LEFT> <RIGHT> <steps> <prec in bits>.\n",argv[0]);
      return 0;
    }

  double top,bottom,left,right;
  int64_t steps,prec,steps1;
  top=atof(argv[1]);
  bottom=atof(argv[2]);
  left=atof(argv[3]);
  right=atof(argv[4]);
  steps=atol(argv[5]);
  prec=atol(argv[6]);
  
  acb_t a,b,res,tmp;

  acb_init(a);acb_init(b);acb_init(res);acb_init(tmp);

  
  arb_set_d(acb_imagref(a),bottom);
  arb_set_d(acb_imagref(b),top);
  arb_set_d(acb_realref(a),right);
  arb_set_d(acb_realref(b),right);
  steps1=steps;
  my_int(res,a,b,steps1,prec); // bottom right -> top right
  while(true)
    {
      if(arb_is_finite(acb_realref(res)))
	if(mag_cmp_2exp_si(arb_radref(acb_realref(res)),0)<0) // error < 1
	  break;
      printf("Insufficient convergence with steps = %ld.",steps1);
      steps1*=2;
      printf(" Iterating I1 with steps = %ld.\n",steps1);
      my_int(res,a,b,steps1,prec);
    }
  printf("I1 from ");acb_printd(a,10);printf(" to ");acb_printd(b,10);
  printf(" returned ");acb_printd(res,20);printf("\n");
  
  
  arb_set_d(acb_imagref(a),top);
  arb_set_d(acb_imagref(b),top);
  arb_set_d(acb_realref(a),right);
  arb_set_d(acb_realref(b),left);
  steps1=steps*256;
  my_int(tmp,a,b,steps1,prec); // top right -> top left
  while(true)
    {
      if(acb_is_finite(tmp))
	if(mag_cmp_2exp_si(arb_radref(acb_realref(tmp)),0)<0) // error < 1
	  break;
      printf("Insufficient convergence with steps = %ld.",steps1);
      steps1*=2;
      printf(" Iterating I2 with steps = %ld.\n",steps1);
      my_int(tmp,a,b,steps1,prec);
    }
  printf("I2 from ");acb_printd(a,10);printf(" to ");acb_printd(b,10);
  printf(" returned ");acb_printd(tmp,20);printf("\n");
  acb_add(res,res,tmp,prec);

  arb_set_d(acb_imagref(a),top);
  arb_set_d(acb_imagref(b),bottom);
  arb_set_d(acb_realref(a),left);
  arb_set_d(acb_realref(b),left);
  steps1=steps*1024;
  my_int(tmp,a,b,steps1,prec); // top left -> bottom left
  while(true)
    {
      if(acb_is_finite(tmp))
	if(mag_cmp_2exp_si(arb_radref(acb_realref(tmp)),0)<0) // error < 1
	  break;
      printf("Insufficient convergence with steps = %ld.",steps1);
      steps1*=2;
      printf(" Iterating I3 with steps = %ld.\n",steps1);
      my_int(tmp,a,b,steps1,prec);
    }
  printf("I3 from ");acb_printd(a,10);printf(" to ");acb_printd(b,10);
  printf(" returned ");acb_printd(tmp,20);printf("\n");
  acb_add(res,res,tmp,prec);

  arb_set_d(acb_imagref(a),bottom);
  arb_set_d(acb_imagref(b),bottom);
  arb_set_d(acb_realref(a),left);
  arb_set_d(acb_realref(b),right);
  steps1=steps*64;
  my_int(tmp,a,b,steps1,prec); // bottom left -> bottom right
  while(true)
    {
      if(acb_is_finite(tmp))
	if(mag_cmp_2exp_si(arb_radref(acb_realref(tmp)),0)<0) // error < 1
	  break;
      printf("Insufficient convergence with steps = %ld.",steps1);
      steps1*=2;
      printf(" Iterating I4 with steps = %ld.\n",steps1);
      my_int(tmp,a,b,steps1,prec);
    }
  printf("I4 from ");acb_printd(a,10);printf(" to ");acb_printd(b,10);
  printf(" returned ");acb_printd(tmp,20);printf("\n");
  acb_add(res,res,tmp,prec);

  printf("Total integral = ");acb_printd(res,20);printf("\n");

  arb_t two_pi;
  arb_const_pi(two_pi,prec);
  arb_mul_2exp_si(two_pi,two_pi,1);

  acb_div_arb(tmp,res,two_pi,prec);
  acb_mul_onei(tmp,tmp);
  acb_neg(tmp,tmp);

  printf("I/(2 pi i) = ");acb_printd(tmp,20);printf("\n");
  
  return 0;
}
