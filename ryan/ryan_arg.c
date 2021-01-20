#include "inttypes.h"
#include "stdbool.h"
#include "acb.h"
#include "../quad/quad.h"

void f(acb_t res, acb_t z, int64_t prec)
{
  static bool init=false;
  static acb_t lz,one_z,tmp1,tmp2,tmp3,tmp4;
  if(!init)
    {
      init=true;
      acb_init(lz);
      acb_init(one_z);
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(tmp3);
      acb_init(tmp4);
    }
  //printf("In f with z = ");acb_printd(z,20);printf("\n");
  acb_set_ui(tmp1,1);
  acb_sub(one_z,tmp1,z,prec); // (1-z)
  acb_log(lz,z,prec); // log z
  //printf("log z = ");acb_printd(lz,20);printf("\n");
  acb_mul_2exp_si(tmp1,one_z,2); // 4(1-z)
  acb_mul_ui(tmp2,z,5,prec);
  acb_add_ui(tmp3,tmp2,1,prec); // (1+5z)
  acb_mul(tmp2,tmp3,lz,prec);
  acb_div(tmp3,tmp2,tmp1,prec); // [(1+5z) log z]/4(1-z)
  //printf("(1+5z)logz/(4(1-z)) = ");acb_printd(tmp3,20);printf("\n");
  acb_mul_ui(tmp2,one_z,3,prec); // 3(1-z)
  acb_add_ui(tmp1,z,1,prec);
  acb_mul(tmp4,tmp1,lz,prec); // (1+z) log z
  acb_add(tmp1,tmp4,tmp2,prec); // 3(1-z)+(1+z)log z
  acb_mul(tmp2,z,lz,prec); // z log z
  acb_div(tmp4,tmp2,tmp1,prec);
  acb_neg(tmp4,tmp4);
  acb_log(tmp1,tmp4,prec);
  //printf("log(-zlogz/(.)) = ");acb_printd(tmp1,20);printf("\n");
  acb_add(tmp4,tmp3,tmp1,prec);
  acb_set_d(tmp1,1.5);
  acb_add(res,tmp4,tmp1,prec);
  return;
}

void fd(acb_t res, acb_t z, int64_t prec)
{
  static bool init=false;
  static acb_t lz,one_z,tmp1,tmp2,tmp3,num,den,t1,t2;

  if(!init)
    {
      init=true;
      acb_init(lz);
      acb_init(one_z);
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(tmp3);
      acb_init(num);
      acb_init(den);
      acb_init(t1);
      acb_init(t2);
    }  
  //printf("In f' with z = ");acb_printd(z,20);printf("\n");
  acb_log(lz,z,prec); // log z
  acb_set_ui(tmp1,1); 
  acb_sub(one_z,tmp1,z,prec); // (1-z)
  acb_mul(tmp1,z,lz,prec);
  acb_mul_ui(tmp2,tmp1,6,prec); // 6z log z
  acb_sqr(tmp1,z,prec);
  acb_mul_ui(tmp3,tmp1,5,prec);
  acb_sub(tmp1,tmp2,tmp3,prec); // 6zlog z-5z^2
  acb_mul_2exp_si(tmp2,z,2); // 4z
  acb_add_ui(tmp3,tmp2,1,prec); // 4z+1
  acb_add(num,tmp1,tmp3,prec); // 6zlog z-5z^2+4z+1
  acb_sqr(tmp1,one_z,prec);
  acb_mul(den,tmp1,z,prec);
  acb_mul_2exp_si(den,den,2); // 4z(1-z)^2
  acb_div(t1,num,den,prec); // term 1 

  acb_add_ui(tmp1,lz,3,prec);
  acb_mul(tmp2,tmp1,lz,prec); // log^2 + 3 log
  acb_mul_ui(tmp1,one_z,3,prec); // 3(1-z)
  acb_add(num,tmp1,tmp2,prec); // log^2 + 3 log + 3(1-z)
  acb_add_ui(tmp2,z,1,prec); // z + 1
  acb_mul(tmp3,tmp2,lz,prec); // zlog z + log z
  acb_add(tmp2,tmp1,tmp3,prec); // 3(1-z)+zlog z+log z
  acb_mul(tmp1,tmp2,z,prec);
  acb_mul(den,tmp1,lz,prec); // (3(1-z)+zlog z+log z) z log z
  acb_div(t2,num,den,prec);

  acb_add(res,t1,t2,prec);
  //printf("returning ");acb_printd(z,20);printf("\n");

}

#define LN_R (-4)

// going to integrate from 0 to 2
// result in circle based at 1 of size 2^LN_R
void gam(acb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static arb_t s,c;
  if(!init)
    {
      init=true;
      arb_init(s);
      arb_init(c);
    }
  arb_sin_cos_pi(s,c,t,prec);
  acb_set_arb_arb(res,c,s);
  acb_mul_2exp_si(res,res,LN_R);
  acb_add_ui(res,res,1,prec);
}

void acb_gam(acb_t res, const acb_t t, int64_t prec)
{
  static bool init=false;
  static acb_t s,c;
  if(!init)
    {
      init=true;
      acb_init(s);
      acb_init(c);
    }
  acb_sin_cos_pi(s,c,t,prec);
  acb_mul_onei(s,s);
  acb_add(res,c,s,prec); // cos t + i sin t
  acb_mul_2exp_si(res,res,LN_R);
  acb_add_ui(res,res,1,prec);
}

void contour_fun(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_t gt,dgt,fz,fdz;
  static arb_t pi;
  if(!init)
    {
      init=true;
      acb_init(gt);
      acb_init(dgt);
      acb_init(fz);
      acb_init(fdz);
      arb_init(pi);
      arb_const_pi(pi,prec);
    }
      
  gam(gt,t,prec); // gamma(t)=1+exp(pi i t)/blah
  acb_sub_ui(dgt,gt,1,prec);
  acb_mul_arb(dgt,dgt,pi,prec);
  acb_mul_onei(dgt,dgt); // gamma'(t) = pi i exp(pi i t)/blah
  f(fz,gt,prec);
  fd(fdz,gt,prec);
  acb_div(gt,fdz,fz,prec);
  acb_mul(fz,gt,dgt,prec);
  arb_set(res,acb_imagref(fz));
}

void acb_contour_fun(acb_t res, const acb_t z, int64_t prec)
{
  static bool init=false;
  static acb_t gt,dgt,fz,fdz;
  static arb_t pi;
  if(!init)
    {
      init=true;
      acb_init(gt);
      acb_init(dgt);
      acb_init(fz);
      acb_init(fdz);
      arb_init(pi);
      arb_const_pi(pi,prec);
    }
  //printf("In acb_contour fun with z = ");acb_printd(z,20);printf("\n");
  acb_gam(gt,z,prec); // gamma(t)
  //printf("gam(z) = ");acb_printd(gt,20);printf("\n");
  acb_sub_ui(dgt,gt,1,prec);
  acb_mul_arb(dgt,dgt,pi,prec);
  acb_mul_onei(dgt,dgt); // gamma'(t)
  //printf("gam'(z) = ");acb_printd(dgt,20);printf("\n");
  f(fz,gt,prec);
  fd(fdz,gt,prec);
  acb_div(gt,fdz,fz,prec);
  //printf("f(gam(z)) = ");acb_printd(fz,20);
  //printf("\nf'(gam(z)) = ");acb_printd(fdz,20);
  //printf("\nf'/f(gam(z)) = ");acb_printd(gt,20);printf("\n"); 
  acb_mul(res,gt,dgt,prec);
  //printf("Returning ");acb_printd(res,20);printf("\n");
}

#define PREC (200)

int main()
{
  /*
  acb_t z,fz;
  acb_init(fz);
  acb_init(z);

  arb_set_d(acb_realref(z),1.25);
  arb_set_d(acb_imagref(z),0.25);
  printf("z = ");acb_printd(z,20);printf("\n");
  f(fz,z,PREC);
  printf("f(z)z = ");acb_printd(fz,20);printf("\n");
  fd(fz,z,PREC);
  printf("f'(z) = ");acb_printd(fz,20);printf("\n");
  */
  
  arb_t res,low,high,maxd;
  arb_init(res);
  arb_init(low);arb_init(high);arb_init(maxd);
  arb_set_ui(low,0);
  arb_set_ui(high,2);
  
  //arb_maxd(maxd,acb_contour_fun,low,high,PREC,1ll<<20);
  arb_set_d(maxd,500);
  printf("Max D = ");arb_printd(maxd,20);printf("\n");
  
  molin_int(res,10,contour_fun,maxd,low,high,PREC);
  printf("Int = ");arb_printd(res,20);printf("\n");

  arb_const_pi(low,PREC);
  arb_div(high,res,low,PREC);
  arb_mul_2exp_si(high,high,-1);
  printf("Z-P = ");arb_printd(high,60);printf("\n");
  fmpz_t ZP;
  fmpz_init(ZP);
  if(!arb_get_unique_fmpz(ZP,high))
    printf("Z-P did not bracket a unique integer.\n");
  else
    {
      printf("Z-P = ");
      fmpz_print(ZP);
      printf("\n");
    }
  
  return 0;
}
