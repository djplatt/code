#include "mpfr.h"
#include "gmp.h"
#include "flint/acb.h"
#include "inttypes.h"
#include <stdbool.h>

#define N_A_TMPS (3)
void ak(arb_ptr res, const arb_ptr sinh_kh, int64_t k, const arb_ptr h, uint64_t prec)
{
  static bool ak_setup=false;
  static arb_t a_kh,a_tmp[N_A_TMPS];
  if(!ak_setup)
    {
      arb_init(a_kh);
      for(int64_t i=0;i<N_A_TMPS;i++)
	arb_init(a_tmp[i]);
      ak_setup=true;
    }
  arb_mul_si(a_kh,h,k,prec);
  arb_cosh(a_tmp[0],a_kh,prec);
  arb_mul(a_tmp[1],a_tmp[0],h,prec);
  arb_sinh(sinh_kh,a_kh,prec);
  arb_cosh(a_tmp[0],sinh_kh,prec);
  arb_pow_ui(a_tmp[2],a_tmp[0],2,prec);
  arb_div(res,a_tmp[1],a_tmp[2],prec);
}

void xk(arb_ptr res, const arb_ptr sinh_kh, uint64_t prec)
{
  arb_tanh(res,sinh_kh,prec);
}

// compute maximum f1(z)*f2(z) around the circle |z|=2
void arb_maxd2(arb_ptr res, 
	      void (*f1)(acb_t,const acb_t,int64_t), 
	      void (*f2)(acb_t,const acb_t,int64_t), 
	      const arb_ptr low, const arb_ptr hi,
	      int64_t prec,
	      uint64_t steps)
{
  arb_t theta,s,c,hi_low;
  arb_init(theta);arb_init(s);arb_init(c);arb_init(hi_low);
  arb_sub(s,hi,low,prec);
  arb_mul_2exp_si(hi_low,s,-1); // (hi/low)/2
  acb_t maxd,z,tmp0,tmp1;
  acb_init(maxd);acb_init(z);acb_init(tmp0);acb_init(tmp1);
  mpfr_t ma,mb;
  mpfr_init(ma);mpfr_init(mb);
  double da=0.0,delta=2.0/(double) steps,db=delta;
  mpfr_set_d(ma,0.0,GMP_RNDN);
  mpfr_set_d(mb,delta,GMP_RNDN);
  arb_zero(res);
  while(db<=2.0)
    {
      arb_set_interval_mpfr(theta,ma,mb,prec);
      arb_sin_cos_pi(s,c,theta,prec);
      arb_mul_2exp_si(acb_realref(z),c,1);
      arb_mul_2exp_si(acb_imagref(z),s,1); // z=2exp(i theta)
      acb_mul_arb(tmp0,z,hi_low,prec);
      arb_add(s,acb_realref(tmp0),hi,prec);
      arb_sub(acb_realref(z),s,hi_low,prec);
      arb_swap(acb_imagref(z),acb_imagref(tmp0));
      f1(tmp0,z,prec);
      f2(tmp1,z,prec);
      acb_mul(z,tmp0,tmp1,prec);
      acb_abs(s,z,prec);
      arb_mul(c,s,hi_low,prec);
      arb_union(s,res,c,prec);
      arb_swap(res,s);
      db+=delta;
      mpfr_swap(ma,mb);
      mpfr_set_d(mb,db,GMP_RNDN);
    }
  arb_get_interval_mpfr(ma,mb,res);
  arb_set_interval_mpfr(res,mb,mb,prec);
  arb_clear(theta);arb_clear(s);arb_clear(c);arb_clear(hi_low);
  acb_clear(maxd);acb_clear(z);acb_clear(tmp0);acb_clear(tmp1);
  mpfr_clear(ma);mpfr_clear(mb);
}

void acb_id(acb_t res, const acb_t x, int64_t prec)
{
  acb_set_ui(res,1);
}
// single function version
void arb_maxd(arb_ptr res, 
	      void (*f1)(acb_t,const acb_t,int64_t), 
	      const arb_ptr low, const arb_ptr hi,
	      int64_t prec,
	      uint64_t steps)
{
  arb_maxd2(res,f1,acb_id,low,hi,prec,steps);
}

// compute maxd*exp(4-5n/log(5n))
// per Molin maxd is max f(z) |z|<=2
// use this version when you have computed maxd as a double perhaps by hand
void comp_errord(arb_t err, double maxd, int64_t n, int64_t prec)
{
  arb_t tmp0,tmp1,tmp2;
  arb_init(tmp0);arb_init(tmp1);arb_init(tmp2); 
  arb_set_d(tmp0,5*n); // 5n
  arb_log(tmp1,tmp0,prec); // log 5n
  arb_ui_div(tmp2,5*n,tmp1,prec); // 5n/log(5n)
  arb_sub_ui(tmp1,tmp2,4,prec); // 5n/log(5n)-4
  arb_neg(tmp0,tmp1); // 4-5n/log(5n)
  arb_exp(tmp1,tmp0,prec); // exp (4-5n/log(5n))
  arb_set_d(tmp0,maxd);
  arb_mul(err,tmp0,tmp1,prec);
  arb_clear(tmp0);arb_clear(tmp1);arb_clear(tmp2);
}

// here maxd is an arb_t, probably computed by arb_maxd
void comp_error(arb_t err, const arb_t maxd, int64_t n, int64_t prec)
{
  arb_t tmp0,tmp1,tmp2; 
  arb_init(tmp0);arb_init(tmp1);arb_init(tmp2); 
  arb_set_d(tmp0,5*n); // 5n
  arb_log(tmp1,tmp0,prec); // log 5n
  arb_ui_div(tmp2,5*n,tmp1,prec); // 5n/log(5n)
  arb_sub_ui(tmp1,tmp2,4,prec); // 5n/log(5n)-4
  arb_neg(tmp0,tmp1); // 4-5n/log(5n)
  arb_exp(tmp1,tmp0,prec); // exp (4-5n/log(5n))
  arb_mul(err,maxd,tmp1,prec);
  arb_clear(tmp0);arb_clear(tmp1);arb_clear(tmp2);
}

// integrate f1(t)*f2(t) from lo to hi using 2*n+1 points
// Molin Theorem 1.1
// maxd is max f1(z)*f2(z) for |z|<=2
void molin_int2(arb_ptr res, int64_t n, 
	       void (*f1)(arb_t,const arb_t,int64_t), 
	       void (*f2)(arb_t,const arb_t,int64_t), 
	       const arb_ptr maxd,
	       const arb_ptr low, const arb_t hi,
	       int64_t prec)
{
  arb_t h,tmp0,tmp1,tmp2,tmp3,tmp4,err;
  arb_init(h);arb_init(tmp0);arb_init(tmp1);arb_init(tmp2);arb_init(tmp3);
  arb_init(err);
  comp_error(err,maxd,n,prec);
  //arb_init(tmp4);
  arb_set_d(tmp0,5.0*n);
  arb_log(tmp1,tmp0,prec);
  arb_div_ui(h,tmp1,n,prec);
  arb_t a_k,x_k,sinh_kh,fx_k;
  arb_init(a_k);arb_init(x_k);arb_init(sinh_kh);arb_init(fx_k);
  arb_zero(tmp1);
  arb_sub(tmp2,hi,low,prec);
  arb_mul_2exp_si(tmp3,tmp2,-1); // (b-a)/2
  for(int64_t k=-n;k<=n;k++)
    {
      ak(a_k,sinh_kh,k,h,prec);
      xk(x_k,sinh_kh,prec);
      arb_mul(tmp0,x_k,tmp3,prec); // t*(hi-lo)/2
      arb_add(x_k,tmp0,hi,prec); // t*(hi-lo)/2+hi
      arb_sub(tmp0,x_k,tmp3,prec); // t*(hi-lo)/2+hi-(hi-lo)/2
      f1(fx_k,tmp0,prec);
      arb_mul(res,fx_k,a_k,prec);
      f2(fx_k,tmp0,prec);
      arb_mul(tmp2,res,fx_k,prec);

      arb_add(tmp0,tmp1,tmp2,prec);
      arb_set(tmp1,tmp0);
      //arb_printd(tmp1,20);printf("\n");
    }
  arb_mul(res,tmp1,tmp3,prec); // *(hi-lo)/2
  arb_add_error(res,err);
  arb_clear(h);arb_clear(tmp0);arb_clear(tmp1);arb_clear(tmp2);
  arb_clear(tmp3);arb_clear(err);
  arb_clear(a_k);arb_clear(x_k);arb_clear(sinh_kh);arb_clear(fx_k);
}

void arb_id (arb_t res, const arb_t x, int64_t prec)
{
  arb_set_ui(res,1);
}

// a version with just f1, uses f2(x)->1
void molin_int(arb_ptr res, int64_t n, 
	       void (*f1)(arb_t,const arb_t,int64_t), 
	       const arb_ptr maxd,
	       const arb_ptr low, const arb_t hi,
	       int64_t prec)
{
  molin_int2(res,n,f1,arb_id,maxd,low,hi,prec);
}
