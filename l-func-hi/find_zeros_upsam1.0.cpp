//
// find_zeros1.0.cpp
//
// reads a series of files produced by l-func
// containing equally spaced values of \Lambda_\chi(1/2+it)
// Conjugate characters are presented and handled together.
//
// Use Turing's method to determine how many zeros there should be
// Count the zeros. Upsample until none are missing.
// Then locate the zeros to +/- 2^{-102}
//
// Upsamples on Turing Zone as well. Does not do stat points

#include <iostream>
#include <cstdlib>
#include "characters.h"

using namespace std;

#define mpfr_print_str(str,x) {printf("%s",str);mpfr_out_str(stdout,10,0,x,GMP_RNDN);printf("\n");}
#define PREC (300)
#define T0 ((double) 200.0) // must be >50.0 for Turing Method to work
#define T1 ((double) 210.0) // we use T0..T1 for Turing Zone
#define T2 ((double) 220.0) // we use this region to interpolate on Turing Zone if necessary

#define OP_ACC (101) // Output is correct to +/- 2^(-OP_ACC-1)
#define Ns (150) // how many points each side of t0 to sample 
#define H ((double) 25.0/64.0) // for Gaussian
#define H2 ((double) H*H)
// error from interpolation (aliasing and truncation of sum)
#define inter_err ((double) 2.8e-44) // valid Ns=150 q<=100,000 t1<=210.0 B=A/2=64/5 H=25/64 
#define MAX_Q ((uint64_t) 100000)
#define SAMPLE_RATE ((double) 5.0/128.0)

// defined so that POS&NEG == 0 but no other combination
typedef int sign_t;
#define POS (1)
#define NEG (2)
#define UNK (3)

typedef int dir_t;
#define UP (1)
#define DOWN (2)

mpfi_t im_t,im_1,ln_tmp;
mpfi_t im_err,im_2,im_t0,im_t1;
mpfi_t tm_t0,tm_t1,tm_res,tm_h,tm_tmp,tm_tmp1,mpfi_pi,mpfi_ln_pi;
mpfr_t tm_mpfr;


void turing_setup()
{
  mpfi_init(im_t);
  mpfi_init(im_1);
  mpfi_init(ln_tmp);
  mpfi_init(im_err);
  mpfi_init(im_2);
  mpfi_init(im_t0);
  mpfi_init(im_t1);
  mpfi_init(tm_t0);
  mpfi_init(tm_t1);
  mpfi_init(tm_res);
  mpfi_init(tm_h);
  mpfi_init(tm_tmp);
  mpfi_init(tm_tmp1);
  mpfi_init(mpfi_pi);
  mpfi_init(mpfi_ln_pi);
  mpfr_init(tm_mpfr);
  mpfi_const_pi(mpfi_pi);
  mpfi_log(mpfi_ln_pi,mpfi_pi);
}

inline sign_t sign(mpfr_ptr x)
{
  int s=mpfr_sgn(x);
  if(s<0)
    return(NEG);
  if(s>0)
    return(POS);
  return(UNK);
}

inline sign_t sign(mpfi_ptr x)
{
  if(mpfi_is_neg(x))
    return(NEG);
  if(mpfi_is_pos(x))
    return(POS);
  return(UNK);
}

inline dir_t dir(mpfr_ptr left, mpfr_ptr right)
{
  int cp=mpfr_cmp(left,right);
  if(cp<0)
    return(UP);
  if(cp>0)
    return(DOWN);
  return(UNK);
}

inline dir_t dir(mpfi_ptr left, mpfi_ptr right)
{
  int cp=mpfi_cmp(left,right);
  if(cp<0)
    return(UP);
  if(cp>0)
    return(DOWN);
  return(UNK);
}

mpfr_t mi_pt,mi_mt0,mp_pt,mp_tmp,mp_dt,e_tmp,e_th2;
mpfr_t mpfr_pi,si_tmp,si_tmp1,mi_tmp;
mpfr_t exps[Ns];
mpfr_t csp_l,csp_r,csp_m,csp_tmp,csp_pt,csp_lt,csp_mt,csp_rt,csp_tmp1,csp_pt1;
mpfr_t csp_res1,dz_tmp1,dz_tmp2;
mpfi_t csp_res2,csp_res3;
mpfr_t rfx3,rx3,rfx4,rtmp1,rtmp2,rtmp3,cetmp;
mpfi_t mpfi_inter_err;
mpfi_t mmi_pt,mmi_mt0,mmp_pt,mmp_tmp,mmp_dt,me_tmp,me_th2;
mpfi_t msi_tmp,msi_tmp1,mmi_tmp;
mpfi_t mexps[Ns],fz_tmp1,fz_tmp2,dz_res1,dz_res2,dz_tmp3,dz_tmp4;
double one_click;
mpfi_t lam_click,half_lam_click;
mpfr_t pc_lt,pc_lf,pc_rt,pc_rf;
mpfr_t zc,zd,ze,zfc,zxm,ztmp,ztmp1,ztmp2,zp,zq,zr,zs,zmin1,zmin2,ztol1;
double zbrent_tol,EPS;
#define ZBRENT_EXTRA_TOL (0)

void inter_setup(double one_over_A)
{
  mpfr_init(mi_pt);
  mpfr_init(mi_mt0);
  mpfr_init(mi_tmp);
  mpfr_init(mp_pt);
  mpfr_init(mp_tmp);
  mpfr_init(mp_dt);
  mpfr_init(e_tmp);
  mpfr_init(e_th2);
  mpfr_set_d(e_th2,-0.5,GMP_RNDN);
  mpfr_div_d(e_th2,e_th2,H2,GMP_RNDN);// -1/2H^
  mpfr_init(mpfr_pi);
  mpfr_const_pi(mpfr_pi,GMP_RNDN);
  mpfr_init(si_tmp);
  mpfr_init(si_tmp1);
  for(uint64_t i=0;i<Ns;i++)
    {
      mpfr_init(exps[i]);
      mpfr_mul_d(exps[i],mpfr_pi,(double) i*one_over_A/(-2.0),GMP_RNDN);
    }
  mpfr_init(csp_l);
  mpfr_init(csp_m);
  mpfr_init(csp_r);
  mpfr_init(csp_tmp);
  mpfr_init(csp_pt);
  mpfr_init(csp_lt);
  mpfr_init(csp_mt);
  mpfr_init(csp_rt);
  mpfr_init(csp_tmp1);
  mpfr_init(csp_pt1);
  mpfr_init(csp_res1);
  mpfi_init(csp_res2);
  mpfi_init(csp_res3);
  mpfr_init(rfx3);
  mpfr_init(rx3);
  mpfr_init(rfx4);
  mpfr_init(rtmp1);
  mpfr_init(rtmp2);
  mpfr_init(rtmp3);
  mpfr_init(cetmp);
  mpfi_init(mpfi_inter_err);
  mpfi_init(mmi_pt);
  mpfi_init(mmi_mt0);
  mpfi_init(mmi_tmp);
  mpfi_init(mmp_pt);
  mpfi_init(mmp_tmp);
  mpfi_init(mmp_dt);
  mpfi_init(me_tmp);
  mpfi_init(me_th2);
  mpfi_set_d(me_th2,-0.5);
  mpfi_div_d(me_th2,me_th2,H2);// -1/2H^
  mpfi_init(msi_tmp);
  mpfi_init(msi_tmp1);
  for(uint64_t i=0;i<Ns;i++)
    {
      mpfi_init(mexps[i]);
      mpfi_mul_d(mexps[i],mpfi_pi,(double) i*one_over_A/(-2.0));
    }

  mpfi_set_d(mpfi_inter_err,inter_err);
  mpfi_neg(msi_tmp,mpfi_inter_err);
  mpfi_put(mpfi_inter_err,msi_tmp);
  mpfi_init(fz_tmp1);
  mpfi_init(fz_tmp2);
  mpfr_init(dz_tmp1);
  mpfr_init(dz_tmp2);
  mpfi_init(dz_res1);
  mpfi_init(dz_res2);
  mpfi_init(dz_tmp3);
  mpfi_init(dz_tmp4);
  mpfr_init(pc_lt);
  mpfr_init(pc_lf);
  mpfr_init(pc_rt);
  mpfr_init(pc_rf);
  zbrent_tol=1.0;
  for(uint64_t i=0;i<OP_ACC+ZBRENT_EXTRA_TOL;i++)
    zbrent_tol*=0.5;
  EPS=1.0;
  for(uint64_t i=0;i<PREC;i++)
    EPS*=0.5;
 
  mpfr_init(zp);
  mpfr_init(zq);
  mpfr_init(zr);
  mpfr_init(zs);
  mpfr_init(ztmp);
  mpfr_init(ztmp1);
  mpfr_init(ztmp2);
  mpfr_init(zmin1);
  mpfr_init(zmin2);
  mpfr_init(zc);
  mpfr_init(zd);
  mpfr_init(ze);
  mpfr_init(zfc);
  mpfr_init(zxm);
  mpfr_init(ztol1);
  one_click=one_over_A;
  for(uint64_t i=0;i<OP_ACC;i++)
    one_click*=0.5;
  mpfi_init(lam_click);
  mpfi_init(half_lam_click);
  mpfi_set_ui(half_lam_click,1);
  mpfi_div_d(half_lam_click,half_lam_click,one_over_A);
  mpfi_div_2ui(half_lam_click,half_lam_click,OP_ACC+1);
  mpfi_mul_2ui(lam_click,half_lam_click,1);
}

/*
#define is_exact(f) (f==MPFI_FLAGS_BOTH_ENDPOINTS_EXACT)


void scale_t(mpfr_ptr res, mpfr_ptr t, double one_over_A)
{
  mpfr_mul_d(res,t,one_over_A,GMP_RNDN);
  mpfr_mul_2ui(res,res,OP_ACC,GMP_RNDN);
  mpfr_add_d(res,res,0.5,GMP_RNDN);
  mpfr_floor(res,res);
}


void bracket_t(mpfi_ptr res1, mpfi_ptr res2, mpfr_t btn, mpfr_ptr t_now, double one_over_A)
{
  scale_t(btn,t_now,one_over_A);
  if(!is_exact(mpfi_set_fr(res1,btn)))
    {
      printf("Wierd fatal error in bracket_t. exiting.\n");
      exit(0);
    }
  mpfi_add_d(res2,res1,0.5);
  mpfi_sub_ui(res1,res2,1);
  mpfi_div_2ui(res1,res1,OP_ACC);
  mpfi_div_2ui(res2,res2,OP_ACC);
  mpfi_div_d(res1,res1,one_over_A);
  mpfi_div_d(res2,res2,one_over_A);
}
*/
//
void mpfr_exp_term(mpfr_ptr res, mpfr_ptr dt, uint64_t ptr)
{
  mpfr_mul(res,dt,dt,GMP_RNDN);
  mpfr_mul(e_tmp,res,e_th2,GMP_RNDN);
  mpfr_add(e_tmp,e_tmp,exps[ptr],GMP_RNDN); // subtract pi t/2 to adjust for negative t, exps[0]=0
  mpfr_exp(res,e_tmp,GMP_RNDN);
}

void mpfi_exp_term(mpfi_ptr res, mpfi_ptr dt, uint64_t ptr)
{
  mpfi_sqr(res,dt);
  mpfi_mul(me_tmp,res,me_th2);
  mpfi_add(me_tmp,me_tmp,mexps[ptr]); // subtract pi t/2 to adjust for negative t, exps[0]=0
  mpfi_exp(res,me_tmp);
}

// sin(pi*x)/pi*x
// x neq 0
void mpfr_sinc(mpfr_ptr res, mpfr_ptr x)
{
  //printf("in sinc with x=");mpfr_out_str(stdout,10,0,x,GMP_RNDN);printf("\n");
  mpfr_mul(si_tmp,x,mpfr_pi,GMP_RNDN);
  mpfr_sin(si_tmp1,si_tmp,GMP_RNDN);
  mpfr_div(res,si_tmp1,si_tmp,GMP_RNDN);
  //printf("sinc returning ");mpfr_out_str(stdout,10,0,res,GMP_RNDN);printf("\n");
}

void mpfi_sinc(mpfi_ptr res, mpfi_ptr x)
{
  if(mpfi_contains_zero(x))
    {
      printf("Sinc called with interval containing zero. Exiting.\n");
      exit(0);
    }
  mpfi_mul(msi_tmp,x,mpfi_pi);
  mpfi_sin(msi_tmp1,msi_tmp);
  mpfi_div(res,msi_tmp1,msi_tmp);
}


// compute the contribution from a single interpolation point
// if we are using the conjugate (or copy of real) to sample
// negative points, we need to adjust exp(pi t/4) to exp(-pi t/4)
void mpfr_inter_point(mpfr_ptr res, uint64_t ptr, mpfr_ptr t0_ptr, mpfi_t *lam, double one_over_A, bool lam_bar)
{
  //printf("inter_point called with %lu and ",ptr);mpfr_out_str(stdout,10,0,t0_ptr,GMP_RNDN);printf("\n");
  mpfr_sub_ui(mp_tmp,t0_ptr,ptr,GMP_RNDN);
  mpfr_sinc(res,mp_tmp);
  mpfr_mul_d(mp_dt,mp_tmp,one_over_A,GMP_RNDN);
  if(lam_bar)
    mpfr_exp_term(mp_pt,mp_dt,ptr); // adjust result by exp(pi t/4)
  else
    mpfr_exp_term(mp_pt,mp_dt,0);
  mpfr_mul(mp_tmp,res,mp_pt,GMP_RNDN);
  mpfi_get_fr(mp_dt,lam[ptr]);
  mpfr_mul(res,mp_tmp,mp_dt,GMP_RNDN);
  //printf("inter_point of %lu returning ",ptr);mpfr_out_str(stdout,10,0,res,GMP_RNDN);printf("\n");
}

void mpfi_inter_point(mpfi_ptr res, uint64_t ptr, mpfi_ptr t0_ptr, mpfi_t *lam, double one_over_A, bool lam_bar)
{
  mpfi_sub_ui(mmp_tmp,t0_ptr,ptr);
  mpfi_sinc(res,mmp_tmp);
  mpfi_mul_d(mmp_dt,mmp_tmp,one_over_A);
  // if we are using lam_bar to compute lam(-t)
  // then we need to adjust for exp(Pi t/4) factor
  if(lam_bar)
    mpfi_exp_term(mmp_pt,mmp_dt,ptr); // adjust result by exp(pi t/4)
  else
    mpfi_exp_term(mmp_pt,mmp_dt,0); // no adjustment
  mpfi_mul(mmp_tmp,res,mmp_pt);
  mpfi_mul(res,mmp_tmp,lam[ptr]);
}

//
// non-rigorously compute lam(t) by Nyquist Shannon
// t=t_ptr*one_over_A
// We need lam_bar in case t is close to 0
// For real characters, lam=lam_bar
//
void mpfr_inter(mpfr_ptr res, mpfi_t *lam, mpfi_t *lam_bar, mpfr_ptr t0_ptr, double one_over_A)
{
  //printf("mpfr_inter called with ");mpfr_out_str(stdout,10,0,t0_ptr,GMP_RNDN);printf("\n");
  int64_t start_ptr=mpfr_get_si(t0_ptr,GMP_RNDU);
  mpfr_set_ui(res,0,GMP_RNDN);
  for(uint64_t n=0,tp=start_ptr;n<Ns;n++,tp++)
    {
      mpfr_inter_point(mi_pt,tp,t0_ptr,lam,one_over_A,false);
      mpfr_add(res,res,mi_pt,GMP_RNDN);
    }
  for(int64_t n=0,tp=start_ptr-1;n<Ns;n++,tp--)
    {
      mpfr_inter_point(mi_pt,tp,t0_ptr,lam,one_over_A,false);
      mpfr_add(res,res,mi_pt,GMP_RNDN);
      if(tp==0)
	{
	  mpfr_neg(mi_mt0,t0_ptr,GMP_RNDN);
	  n++;
	  tp=1;
	  while(n<Ns)
	    {
	      mpfr_inter_point(mi_pt,tp,mi_mt0,lam_bar,one_over_A,true);
	      mpfr_add(res,res,mi_pt,GMP_RNDN);
	      n++;tp++;
	    }
	  break;
	}
    }
}
//
// rigorously compute lam(t) by Nyquist Shannon
// t=t_ptr*one_over_A
// We need lam_bar in case t is close to 0
// For real characters, lam=lam_bar
//
void mpfi_inter(mpfi_ptr res, mpfi_t *lam, mpfi_t *lam_bar, mpfi_ptr t0_ptr, double one_over_A)
{
  //mpfi_print_str("In mpfi_inter with ptr=",t0_ptr);
  int64_t start_ptr=ceil(mpfi_get_d(t0_ptr)); // snap to next point on right
  mpfi_set(res,mpfi_inter_err); // interpolation error (aliasing and truncation)
  for(uint64_t n=0,tp=start_ptr;n<Ns;n++,tp++) // add contribution from Ns points to right
    {
      mpfi_inter_point(mmi_pt,tp,t0_ptr,lam,one_over_A,false);
      mpfi_add(res,res,mmi_pt);
    }
  for(int64_t n=0,tp=start_ptr-1;n<Ns;n++,tp--) // now do Ns points to left
    {                                           // if not enough, use lam_bar
      mpfi_inter_point(mmi_pt,tp,t0_ptr,lam,one_over_A,false);
      mpfi_add(res,res,mmi_pt);
      if(tp==0) // not enough points
	{
	  mpfi_neg(mmi_mt0,t0_ptr);
	  n++;
	  tp=1;
	  while(n<Ns) // use lam_bar for the rest, adjusting exp(pi t/4)
	    {
	      mpfi_inter_point(mmi_pt,tp,mmi_mt0,lam_bar,one_over_A,true);
	      mpfi_add(res,res,mmi_pt);
	      n++;tp++;
	    }
	  break;
	}
    }
}

inline bool mpfr_is_pos(mpfr_ptr x) { return(mpfr_sgn(x)>0);} 
inline bool mpfr_is_neg(mpfr_ptr x) { return(mpfr_sgn(x)<0);} 


#define ZBRENT_ITS (30)
void zbrent(mpfr_ptr res, mpfr_ptr za, mpfr_ptr zfa, mpfr_ptr zb, mpfr_ptr zfb, mpfi_t *lam, mpfi_t *lam_bar, double one_over_A)
{
  mpfr_set(zc,zb,GMP_RNDN);
  mpfr_set(zfc,zfb,GMP_RNDN);
  for(uint64_t it=0;it<ZBRENT_ITS;it++)
    {
      //mpfr_print_str("a=",za);
      //mpfr_print_str("b=",zb);
      //mpfr_print_str("c=",zc);
      //mpfr_print_str("f(a)=",zfa);
      //mpfr_print_str("f(b)=",zfb);
      //mpfr_print_str("f(c)=",zfc);
      if((mpfr_is_pos(zfb)&&mpfr_is_pos(zfc))||(mpfr_is_neg(zfb)&&mpfr_is_neg(zfc)))
	{
	  mpfr_set(zc,za,GMP_RNDN);
	  mpfr_set(zfc,zfa,GMP_RNDN);
	  mpfr_sub(zd,zb,za,GMP_RNDN);
	  mpfr_set(ze,zd,GMP_RNDN);
	}
      if(mpfr_cmp_abs(zfc,zfb)<0)
	{
	  mpfr_set(za,zb,GMP_RNDN);
	  mpfr_set(zb,zc,GMP_RNDN);
	  mpfr_set(zc,za,GMP_RNDN);
	  mpfr_set(zfa,zfb,GMP_RNDN);
	  mpfr_set(zfb,zfc,GMP_RNDN);
	  mpfr_set(zfc,zfa,GMP_RNDN);
	}
      mpfr_mul_d(ztmp,zb,2.0*EPS,GMP_RNDN);
      if(mpfr_is_neg(ztmp))
	{
	  mpfr_neg(ztmp1,ztmp,GMP_RNDN);
	  mpfr_add_d(ztol1,ztmp1,0.5*zbrent_tol,GMP_RNDN);
	}
      else
	mpfr_add_d(ztol1,ztmp,0.5*zbrent_tol,GMP_RNDN);
      //mpfr_print_str("tol1=",ztol1);
      mpfr_sub(zxm,zc,zb,GMP_RNDN);
      mpfr_div_2ui(zxm,zxm,1,GMP_RNDN);
      mpfr_abs(ztmp,zxm,GMP_RNDN);
      //mpfr_print_str("abs(xm)=",ztmp);
      if(mpfr_cmp(ztmp,ztol1)<=0)
	{
	  mpfr_set(res,zb,GMP_RNDN);
	  return;
	}
      mpfr_abs(ztmp,ze,GMP_RNDN);
      if((mpfr_cmp(ztmp,ztol1)>0)&&(mpfr_cmp_abs(zfa,zfb)>0))
	{
	  mpfr_div(zs,zfb,zfa,GMP_RNDN);
	  if(mpfr_cmp(za,zc)==0)
	    {
	      mpfr_mul(zp,zxm,zs,GMP_RNDN);
	      mpfr_mul_2ui(zp,zp,1,GMP_RNDN);
	      mpfr_sub_ui(ztmp,zs,1,GMP_RNDN);
	      mpfr_neg(zq,ztmp,GMP_RNDN);
	    }
	  else
	    {
	      mpfr_div(zq,zfa,zfc,GMP_RNDN);
	      mpfr_div(zr,zfb,zfc,GMP_RNDN);
	      mpfr_sub(zp,zb,za,GMP_RNDN); // b-a
	      mpfr_sub_ui(ztmp,zr,1.0,GMP_RNDN); // r-1
	      mpfr_mul(ztmp2,zp,ztmp,GMP_RNDN); // (b-a)*(r-1)
	      mpfr_sub(zp,zq,zr,GMP_RNDN); // (q-r)
	      mpfr_mul(ztmp,zp,zq,GMP_RNDN);
	      mpfr_mul(ztmp1,zxm,ztmp,GMP_RNDN);
	      mpfr_mul_2ui(ztmp1,ztmp1,1,GMP_RNDN); // 2*xm*q*(q-r)
	      mpfr_sub(ztmp,ztmp1,ztmp2,GMP_RNDN);
	      mpfr_mul(zp,zs,ztmp,GMP_RNDN);
	      mpfr_sub_ui(ztmp,zq,1,GMP_RNDN); // q-1
	      mpfr_sub_ui(ztmp1,zr,1,GMP_RNDN); // r-1
	      mpfr_mul(ztmp2,ztmp1,ztmp,GMP_RNDN); // (q-1)(r-1)
	      mpfr_sub_ui(ztmp,zs,1,GMP_RNDN); // (s-1)
	      mpfr_mul(zq,ztmp,ztmp2,GMP_RNDN); // q=(q-1)(r-1)(s-1)
	    }
	  if(mpfr_is_pos(zp))
	    mpfr_neg(zq,zq,GMP_RNDN);
	  else
	    mpfr_neg(zp,zp,GMP_RNDN);
	  mpfr_mul(ztmp,zq,ztol1,GMP_RNDN);
	  mpfr_mul(zmin1,zxm,zq,GMP_RNDN);
	  mpfr_mul_ui(ztmp1,zmin1,3,GMP_RNDN);
	  if(mpfr_is_pos(ztmp))
	    mpfr_sub(zmin1,ztmp1,ztmp,GMP_RNDN);
	  else
	    mpfr_add(zmin1,ztmp1,ztmp,GMP_RNDN);
	  mpfr_mul(ztmp,ze,zq,GMP_RNDN);
	  if(mpfr_is_neg(ztmp))
	    mpfr_neg(zmin2,ztmp,GMP_RNDN);
	  else
	    mpfr_set(zmin2,ztmp,GMP_RNDN);
	  if(mpfr_cmp(zmin1,zmin2)>0)
	    mpfr_swap(zmin1,zmin2);
	  mpfr_div_2ui(zmin1,zmin1,1,GMP_RNDN);
	  if(mpfr_cmp(zp,zmin1)<0)
	    {
	      mpfr_set(ze,zd,GMP_RNDN);
	      mpfr_div(zd,zp,zq,GMP_RNDN);
	    }
	  else
	    {
	      mpfr_set(zd,zxm,GMP_RNDN);
	      mpfr_set(ze,zd,GMP_RNDN);
	    }
	}
      else
	{
	  mpfr_swap(zd,zxm);
	  mpfr_set(ze,zd,GMP_RNDN);
	}
      mpfr_set(za,zb,GMP_RNDN);
      mpfr_set(zfa,zfb,GMP_RNDN);
      mpfr_abs(ztmp,zd,GMP_RNDN);
      //mpfr_print_str("abs(d)=",ztmp);
      if(mpfr_cmp(ztmp,ztol1)>0)
	mpfr_add(zb,zb,zd,GMP_RNDN);
      else
	{
	  if(mpfr_is_pos(zxm))
	    mpfr_add(zb,zb,ztol1,GMP_RNDN);
	  else
	    mpfr_sub(zb,zb,ztol1,GMP_RNDN);
	}
      mpfr_inter(zfb,lam,lam_bar,zb,one_over_A);
    }
  printf("zbrent failed to converge after %lu iterations. Exiting.\n",ZBRENT_ITS);
  exit(0);
}

/*
#define OP_DEL ((double) exp(-(OP_ACC+10)*log(2))) // if f(t) within this
                                                  // then t within OP_ACC/2  

#define MAX_SECANT_ITS (120)

// use secant method to find root
void secant(mpfr_ptr res, mpfr_ptr x1, mpfr_ptr fx1, mpfr_ptr x2, mpfr_ptr fx2, mpfi_t *lam, mpfi_t *lam_bar, double one_over_A)
{
  uint64_t it_count=0;
  while(true)
    {
      mpfr_sub(rtmp1,fx2,fx1,GMP_RNDN);
      mpfr_sub(rtmp2,x1,x2,GMP_RNDN);
      mpfr_mul(res,rtmp2,fx1,GMP_RNDN);
      mpfr_div(rtmp2,res,rtmp1,GMP_RNDN);
      mpfr_add(rx3,x1,rtmp2,GMP_RNDN);
      mpfr_inter(rfx3,lam,lam_bar,rx3,one_over_A);
      if(fabs(mpfr_get_d(rfx3,GMP_RNDN))<OP_DEL)
	{
	  mpfr_set(res,rx3,GMP_RNDN);
	  return;
	}
      if(mpfr_is_pos(rfx3))
	{
	  if(mpfr_is_neg(fx1)) // use x1 and x3
	    {
	      mpfr_swap(x2,rx3);
	      mpfr_swap(fx2,rfx3);
	    }
	  else // use x2 and x3
	    {
	      mpfr_swap(x1,rx3);
	      mpfr_swap(fx1,rfx3);
	    }
	}
      else
	{
	  if(mpfr_is_pos(fx1)) // use x1 and x3
	    {
	      mpfr_swap(x2,rx3);
	      mpfr_swap(fx2,rfx3);
	    }
	  else // use x2 and x3
	    {
	      mpfr_swap(x1,rx3);
	      mpfr_swap(fx1,rfx3);
	    }
	}
      if(++it_count==MAX_SECANT_ITS)
	{
	  printf("Secant method failed to converge. Exiting.\n");
	  exit(0);
	}
    }
}


#define MAX_RIDDER_ITS (30)
//
// non-rigorous root finding using Ridders method
//
void ridders(mpfr_ptr x4, mpfr_ptr x1, mpfr_ptr fx1, mpfr_ptr x2, mpfr_ptr fx2, mpfi_t *lam, mpfi_t *lam_bar, double one_over_A)
{
  for(uint64_t it_count=0;it_count<MAX_RIDDER_ITS;it_count++)
    {
      mpfr_add(rx3,x1,x2,GMP_RNDN);
      mpfr_div_2ui(rx3,rx3,1,GMP_RNDN);
      mpfr_inter(rfx3,lam,lam_bar,rx3,one_over_A);
      mpfr_mul(rtmp1,fx1,fx2,GMP_RNDN);
      mpfr_mul(rtmp2,rfx3,rfx3,GMP_RNDN);
      mpfr_sub(rtmp3,rtmp2,rtmp1,GMP_RNDN);
      mpfr_sqrt(rtmp1,rtmp3,GMP_RNDN);
      mpfr_div(rtmp2,rfx3,rtmp1,GMP_RNDN);
      if(mpfr_cmp(fx1,fx2)>0)
	mpfr_sub(rtmp1,rx3,x1,GMP_RNDN);
      else
	mpfr_sub(rtmp1,x1,rx3,GMP_RNDN);
      mpfr_mul(rtmp3,rtmp1,rtmp2,GMP_RNDN);
      mpfr_add(x4,rtmp3,rx3,GMP_RNDN);
      mpfr_inter(rfx4,lam,lam_bar,x4,one_over_A);
      if(fabs(mpfr_get_d(rfx4,GMP_RNDN))<OP_DEL)
	return;
      if(mpfr_sgn(rfx4)>0) // f(x4)>0
	{
	  if(mpfr_sgn(rfx3)<0) // f(x3)<0, use x4 and x3
	    {
	      mpfr_swap(x1,rx3);
	      mpfr_swap(fx1,rfx3);
	      mpfr_swap(x2,x4);
	      mpfr_swap(fx2,rfx4);
	    }
	  else
	    {
	      if(mpfr_sgn(fx2)<0) // use x4 and x2
		{
		  mpfr_swap(x1,x4);
		  mpfr_swap(fx1,rfx4);
		}
	      else // use x4 and x1
		{
		  mpfr_swap(x2,x4);
		  mpfr_swap(fx2,rfx4);
		}
	    }
	}
      else // f(x4)<0
	{
	  if(mpfr_sgn(rfx3)>0) // f(x3)>0, use x4 and x3
	    {
	      mpfr_swap(x1,rx3);
	      mpfr_swap(fx1,rfx3);
	      mpfr_swap(x2,x4);
	      mpfr_swap(fx2,rfx4);
	    }
	  else
	    {
	      if(mpfr_sgn(fx2)>0) // use x4 and x2
		{
		  mpfr_swap(x1,x4);
		  mpfr_swap(fx1,rfx4);
		}
	      else // use x4 and x1
		{
		  mpfr_swap(x2,x4);
		  mpfr_swap(fx2,rfx4);
		}
	    }
	}
      if(mpfr_cmp(x1,x2)>0)
	{
	  mpfr_swap(x1,x2);
	  mpfr_swap(fx1,fx2);
	}
      it_count++;
    }
  printf("Ridders failed to converge. Exiting.\n");
  exit(0);
}
*/

mpz_t out_z;
mpfr_t last_t,out_tmp;


void init_out_bytes()
{
  mpz_init(out_z);
  mpfr_init(last_t);
  mpfr_set_ui(last_t,0,GMP_RNDN);
  mpfr_init(out_tmp);
}

uint8_t null_rec[]={0,0,0,0,0,0,0,0,0,0,0,0};

void write_null_record(FILE *outfile)
{
  fwrite(null_rec,sizeof(uint8_t),13,outfile);
}  

void out_bytes(mpfr_ptr t, FILE *outfile)
{
  uint64_t a,i;
  uint32_t b;
  uint8_t c;
  //mpfr_print_str("Outputting zero at 1/2+i",t);
  mpfr_sub(out_tmp,t,last_t,GMP_RNDN);
  mpfr_set(last_t,t,GMP_RNDN);
  mpfr_get_z(out_z,out_tmp,GMP_RNDN);
  a=mpz_get_ui(out_z); // just grabs least sig. long uns int (8)
  fwrite(&a,sizeof(uint64_t),1,outfile);
  mpz_fdiv_q_2exp(out_z,out_z,64);
  i=mpz_get_ui(out_z);
  b=i&0xFFFFFFFF;
  fwrite(&b,sizeof(uint32_t),1,outfile);
  i>>=32;
  if(i>255)
    {
      printf("Argument to out_bytes exceeds 13 bytes. Exiting.\n");
      exit(1);
    }
  c=i;
  fwrite(&c,sizeof(uint8_t),1,outfile);
}


// starting with two mpfr indices into lam (not necessarily integral)
// demonstrating a change of sign,
// use Brent's method to non-rigorously isolate the root
// then rigorously confirm it is within +/- 2^-(OP_ACC+1)
// 
void do_zero(mpfr_ptr lt, mpfr_ptr lf, mpfr_ptr rt, mpfr_ptr rf, mpfi_t *lam, mpfi_t *lam_bar, double one_over_A,FILE *outfile)
{
  //mpfr_out_str("In do zero with lt=",lt);
  zbrent(dz_tmp2,lt,lf,rt,rf,lam,lam_bar,one_over_A);
  // dz_tmp2 is mpfr ptr into lam where we think there is a sign change
  // we want to check +/ 2^{-OP_ACC-1} in t either side
  // this is +/- 2^{-OP_ACC-1}/one_over_A in lam ptr
  mpfi_set_fr(dz_res1,dz_tmp2);
  mpfi_add(dz_res2,dz_res1,half_lam_click); // right hand
  mpfi_sub(dz_res1,dz_res2,lam_click); // left hand

  mpfi_inter(dz_tmp3,lam,lam_bar,dz_res1,one_over_A);
  sign_t lsign=sign(dz_tmp3);

  mpfi_inter(dz_tmp4,lam,lam_bar,dz_res2,one_over_A);
  sign_t rsign=sign(dz_tmp4);

  uint64_t shift_count=1;
  while((lsign&rsign)!=0) // not a valid sign change, must have missed by a bit
    {
      printf("Shifting in do_zero %lu.\n",shift_count++);
      printf("Bracketing didn't work.\n");
      mpfi_print_str("Left hand value=",dz_tmp3);
      mpfi_print_str("Right hand value=",dz_tmp4);
      if(((lsign==POS)&&(mpfi_cmp(dz_tmp3,dz_tmp4)>0))|| // need to move right
	 ((lsign==NEG)&&(mpfi_cmp(dz_tmp3,dz_tmp4)<0)))
	    {
	      mpfi_set(dz_tmp3,dz_tmp4);
	      mpfi_set(dz_res1,dz_res2);
	      lsign=rsign;
	      mpfi_add(dz_res2,dz_res1,lam_click);
	      mpfi_inter(dz_tmp4,lam,lam_bar,dz_res2,one_over_A);
	      rsign=sign(dz_tmp4);
	    }
	  else // need to move left
	    {
	      mpfi_set(dz_tmp4,dz_tmp3);
	      mpfi_set(dz_res2,dz_res1);
	      rsign=lsign;
	      mpfi_sub(dz_res1,dz_res2,lam_click);
	      mpfi_inter(dz_tmp3,lam,lam_bar,dz_res1,one_over_A);
	      lsign=sign(dz_tmp3);
	    }
    }

  // here dz_res1 and dz_res2 are intervals bracketing a zero
  //  specifically right of dz_res1 and left of dz_res2 do it.
  // convert them into t values *2^{OP_ACC+1}
  mpfi_mul_d(dz_res1,dz_res1,one_over_A);
  mpfi_mul_2ui(dz_res1,dz_res1,OP_ACC+1);
  if(mpfi_get_right(dz_tmp1,dz_res1)!=0) // should be exact
    {
      printf("Error getting right end point in do_zeros. Exiting.\n");
      exit(0);
    }
  mpfi_mul_d(dz_res2,dz_res2,one_over_A);
  mpfi_mul_2ui(dz_res2,dz_res2,OP_ACC+1);
  if(mpfi_get_left(dz_tmp2,dz_res2)!=0) // should also be exact
    {
      printf("Error getting left end point in do_zeros. Exiting.\n");
      exit(0);
    }

  // ensure end points are with 2 so mid point is within 1*2^{-OP_ACC-1}
  mpfr_add_ui(dz_tmp1,dz_tmp1,2,GMP_RNDN);
  if(mpfr_cmp(dz_tmp1,dz_tmp2)<0)
    {
      printf("End points more than 2^-{OP_ACC-1} apart. Exiting\n");
      exit(0);
    }
  mpfr_sub_ui(dz_tmp1,dz_tmp1,2,GMP_RNDN);
  mpfr_add(dz_tmp1,dz_tmp1,dz_tmp2,GMP_RNDN);
  mpfr_div_2ui(dz_tmp1,dz_tmp1,2,GMP_RNDN); // average and back to OP_ACC
  out_bytes(dz_tmp1,outfile);
  //printf("Zero found near ");
  //mpfr_div_2ui(dz_tmp1,dz_tmp1,OP_ACC,GMP_RNDN);
  //mpfr_out_str(stdout,10,0,dz_tmp1,GMP_RNDN);
  //printf("\n");
}

// resolve a stat point rigorously.
// if resolve only, find the sign change and return (e.g. stat points in
// Turing Zone). If ! resolve_only, isolate the zeros and write them out.
bool do_stat_point(uint64_t i1, uint64_t i2, uint64_t i3, mpfi_t *lam, mpfi_t *lam_bar, sign_t this_sign, double one_over_A, bool resolve_only,FILE *outfile)
{
  bool solve_left=true; // true if solution is between csp_lt and csp_mt
  /*
  if(this_sign==POS)
    printf("pos ");
  else
    printf("neg ");
  printf("%lu %lu %lu\n",i1,i2,i3);
  */
  mpfi_get_fr(csp_l,lam[i1]);
  mpfi_get_fr(csp_m,lam[i2]);
  mpfi_get_fr(csp_r,lam[i3]);
  mpfr_set_ui(csp_lt,i1,GMP_RNDN);
  mpfr_set_ui(csp_mt,i2,GMP_RNDN);
  mpfr_set_ui(csp_rt,i3,GMP_RNDN);
  while(true)
    {
      //mpfr_print_str("in stat pt    ",csp_lt);
      //mpfr_print_str("              ",csp_mt);
      //mpfr_print_str("              ",csp_rt);

      mpfr_add(csp_tmp,csp_lt,csp_mt,GMP_RNDN);
      mpfr_div_2ui(csp_tmp,csp_tmp,1,GMP_RNDN);
      mpfr_inter(csp_pt,lam,lam_bar,csp_tmp,one_over_A);
      if(this_sign==POS)
	{
	  if(mpfr_sgn(csp_pt)<0) // found our sign change
	    break;
	  if(mpfr_cmp(csp_pt,csp_m)<0) // this is lower than mid
	    {
	      //printf("moving left points\n");
	      mpfr_set(csp_r,csp_m,GMP_RNDN);
	      mpfr_set(csp_m,csp_pt,GMP_RNDN);
	      mpfr_set(csp_rt,csp_mt,GMP_RNDN);
	      mpfr_set(csp_mt,csp_tmp,GMP_RNDN);
	      continue;
	    }
	}
      else
	{      
	  if(mpfr_sgn(csp_pt)>0) // sign change
	    break;
	  if(mpfr_cmp(csp_pt,csp_m)>0) // this is higher than mid
	    {
	      mpfr_set(csp_r,csp_m,GMP_RNDN);
	      mpfr_set(csp_m,csp_pt,GMP_RNDN);
	      mpfr_set(csp_rt,csp_mt,GMP_RNDN);
	      mpfr_set(csp_mt,csp_tmp,GMP_RNDN);
	      continue;
	    }
	}
      // ok so the left two points didn't work out, let's try the right
      mpfr_add(csp_tmp1,csp_mt,csp_rt,GMP_RNDN);
      mpfr_div_2ui(csp_tmp1,csp_tmp1,1,GMP_RNDN);
      mpfr_inter(csp_pt1,lam,lam_bar,csp_tmp1,one_over_A);
      if(this_sign==POS)
	{
	  if(mpfr_sgn(csp_pt1)<0)
	    {
	      solve_left=false;
	      break;
	    }
	  if(mpfr_cmp(csp_pt1,csp_m)<0) // this is lower than mid
	    {
	      //printf("moving right points\n");
	      mpfr_set(csp_l,csp_m,GMP_RNDN);
	      mpfr_set(csp_m,csp_pt1,GMP_RNDN);
	      mpfr_set(csp_lt,csp_mt,GMP_RNDN);
	      mpfr_set(csp_mt,csp_tmp1,GMP_RNDN);
	      continue;
	    }
	}
      else
	{      
	  if(mpfr_sgn(csp_pt1)>0)
	    {
	      solve_left=false;
	      break;
	    }
	  if(mpfr_cmp(csp_pt1,csp_m)>0) // this is higher than mid
	    {
	      mpfr_set(csp_l,csp_m,GMP_RNDN);
	      mpfr_set(csp_m,csp_pt1,GMP_RNDN);
	      mpfr_set(csp_lt,csp_mt,GMP_RNDN);
	      mpfr_set(csp_mt,csp_tmp1,GMP_RNDN);
	      continue;
	    }
	}
      mpfr_set(csp_l,csp_pt,GMP_RNDN);
      mpfr_set(csp_r,csp_pt1,GMP_RNDN);
      mpfr_set(csp_lt,csp_tmp,GMP_RNDN);
      mpfr_set(csp_rt,csp_tmp1,GMP_RNDN);
    }
  if(resolve_only) // we only had to check for a pair of zeros and return
    {
      if(solve_left) // solution was between left and mid
	mpfi_set_fr(fz_tmp2,csp_tmp);
      else
	mpfi_set_fr(fz_tmp2,csp_tmp1);
      mpfi_inter(fz_tmp1,lam,lam_bar,fz_tmp2,one_over_A);
      if(this_sign==POS)
	return(mpfi_is_neg(fz_tmp1));
      else      
	return(mpfi_is_pos(fz_tmp1));
    }
  else // now need to locate the zeros to high accuracy
    {
      if(solve_left) // solution was between left and mid
	{
	  //mpfr_print_str("in solve left ",csp_lt);
	  //mpfr_print_str("              ",csp_tmp);
	  //mpfr_print_str("              ",csp_mt);
	  do_zero(csp_lt,csp_l,csp_tmp,csp_pt,lam,lam_bar,one_over_A,outfile);
	  do_zero(csp_tmp,csp_pt,csp_mt,csp_m,lam,lam_bar,one_over_A,outfile);
	}
      else
	{
	  //mpfr_print_str("in solve right ",csp_mt);
	  //mpfr_print_str("               ",csp_tmp1);
	  //mpfr_print_str("               ",csp_rt);
	  do_zero(csp_mt,csp_m,csp_tmp1,csp_pt1,lam,lam_bar,one_over_A,outfile);
	  do_zero(csp_tmp1,csp_pt1,csp_rt,csp_r,lam,lam_bar,one_over_A,outfile);
	}
      return(true);
    }
}

// find local minima above the line and local maxima below it
// both contravene RH and indicate pairs of missed zeros
// non rigorously assume the pair of zeros will be found at some point
uint64_t num_stat_pts(mpfi_t *f_vec, uint64_t start, uint64_t end)
{
  uint64_t result=0;
  for(uint64_t i=start;i<end-2;i++)
    {
      if(mpfi_is_pos(f_vec[i])&&mpfi_is_pos(f_vec[i+1])&&mpfi_is_pos(f_vec[i+2])
	 &&(mpfi_cmp(f_vec[i],f_vec[i+1])>0)&&(mpfi_cmp(f_vec[i+1],f_vec[i+2])<0))
	{
	  //printf("+ve Stat point at %lu %lu %lu\n",i,i+1,i+2);
	  result+=2;
	  continue;
	}
      if(mpfi_is_neg(f_vec[i])&&mpfi_is_neg(f_vec[i+1])&&mpfi_is_neg(f_vec[i+2])
	 &&(mpfi_cmp(f_vec[i],f_vec[i+1])<0)&&(mpfi_cmp(f_vec[i+1],f_vec[i+2])>0))
	{
	  //printf("-ve Stat point at %lu %lu %lu\n",i,i+1,i+2);
	  result+=2;
	}
    }
  return(result);
}

// used to compute 
// \int\limits_{t_0}^{t_0+h} Im loggamma(3/4+it/2) dt
// needed for Turing's method with Odd characters
void im_int2(mpfi_ptr res, mpfi_ptr t)
{
  mpfi_mul_ui(im_t,t,4); // 4t
  mpfi_div_ui(im_t,t,3); // t/s
  mpfi_atan(res,im_t); // atan(t/s)
  //printf("atan(t/s)=");mpfi_printn(res,30);
  mpfi_mul(res,res,t); // t*atan(t/s)
  mpfi_mul_d(res,res,0.25); // t*atan(t/s)*(s-1/2)
  //printf("-t/4atan(4t)=");mpfi_printn(res,30);
  mpfi_sqr(im_t,t); // t^2
  mpfi_mul_d(im_1,im_t,-3.0/4.0); 
  mpfi_add(res,res,im_1); //-3/4*t^2
  mpfi_mul_ui(im_1,im_t,16); // 16*t^2
  mpfi_div_ui(im_1,im_1,9); // t^2/s^2
  mpfi_add_ui(im_1,im_1,1); // 1+t^2/s^2
  mpfi_log(im_1,im_1); // log(1+t^2/s^2)
  mpfi_mul_d(im_1,im_1,15.0/32.0); // (s^2+s/4) log(1+t^2/s^2)
  mpfi_sub(res,res,im_1);
  mpfi_add_d(res,res,-9.0/64.0); // -s^2/4
  mpfi_add_d(im_t,im_t,9.0/16.0); // t^2+9/16
  mpfi_log(im_1,im_t); // log(t^2+s^2)
  mpfi_mul(im_1,im_t,im_1); // (t^2+s^2)log(t^2+s^2)
  mpfi_mul_d(im_t,im_1,0.25); // log(t^2+s^2)/4
  mpfi_add(res,res,im_t);
}

// to compute
// \int\limits_{t_0}^{t_0+h} Im loggamma(1/4+it/2) dt
// Even characters
void im_int1(mpfi_ptr res, mpfi_ptr t)
{
  mpfi_mul_ui(im_t,t,4); // 4t
  mpfi_atan(res,im_t);
  //printf("atan(4t)=");mpfi_printn(res,30);
  mpfi_mul(res,res,t);
  mpfi_mul_d(res,res,-0.25); // -t/4*atan(4t)
  //printf("-t/4atan(4t)=");mpfi_printn(res,30);
  mpfi_sqr(im_t,t); // t^2
  mpfi_mul_d(im_1,im_t,-3.0/4.0);
  mpfi_add(res,res,im_1); //-3/4*t^2
  mpfi_mul_ui(im_1,im_t,16); // 16t^2
  mpfi_add_ui(im_1,im_1,1);
  mpfi_log(im_1,im_1); // log(1+16t^2)
  mpfi_mul_d(im_1,im_1,1.0/32.0);
  mpfi_add(res,res,im_1);
  mpfi_add_d(res,res,-1.0/64.0);
  mpfi_add_d(im_t,im_t,1.0/16.0); // t^2+1/16
  mpfi_log(im_1,im_t);
  mpfi_mul(im_t,im_t,im_1); // (t^2+1/16)log(t^2+1/16)
  mpfi_mul_d(im_t,im_t,0.25);
  mpfi_add(res,res,im_t);
}

// Compute the integral
void im_int(mpfi_ptr res, mpfi_ptr t0, mpfi_ptr t1, bool even_p)
{
  mpfi_sub(im_err,t1,t0);
  mpfi_div_ui(im_err,im_err,4);
  mpfi_div(im_err,im_err,t0);
  mpfi_neg(im_2,im_err);
  mpfi_put(im_err,im_2); // |error| <= (t1-t0)/4t0
  mpfi_div_ui(im_t0,t0,2);
  mpfi_div_ui(im_t1,t1,2);
  if(even_p)
    {
      im_int1(res,im_t1);
      im_int1(im_2,im_t0);
    }
  else
    {
      im_int2(res,im_t1);
      im_int2(im_2,im_t0);
    }
  mpfi_sub(res,res,im_2);
  mpfi_mul_ui(res,res,2);
  mpfi_add(res,res,im_err);
}

/*
double Nleft_int(long int t0_ptr, long int t1_ptr, mpfi_t *f_vec, double delta)
{
  long int res=0;
  long int ptr=t0_ptr,last_ptr;
  sign_t last_sign=sign(f_vec[ptr++]),this_sign;
  while(last_sign==UNK)
    last_sign=sign(f_vec[ptr++]);
  last_ptr=ptr-1;
  for(;ptr<=t1_ptr;)
    {
      this_sign=sign(f_vec[ptr++]);
      if(this_sign==last_sign) // no sign change here, move on
	{
	  last_ptr=ptr-1;
	  continue;
	}
      if(this_sign==UNK) // might be a sign change coming
	continue;
      // definately a sign change
      //printf("Sign change in Nleft\n");
      res-=(last_ptr-t0_ptr+1);
      last_ptr=ptr-1;
      last_sign=this_sign;
    }
  //printf("Nleft_int returning %f\n",res*delta);
  return(res*delta);
}
*/

double Nright_int(uint64_t t0_ptr, uint64_t t1_ptr, mpfi_t *f_vec, mpfi_t *f_vec_bar, double delta)
{
  int64_t res=0;
  uint64_t ptr=t0_ptr;
  sign_t last_sign=sign(f_vec[ptr++]),this_sign;
  while(last_sign==UNK)
    last_sign=sign(f_vec[ptr++]);
  for(;ptr<=t1_ptr;)
    {
      this_sign=sign(f_vec[ptr++]);
      if(this_sign==last_sign) // no sign change here, move on
	  continue;
      if(this_sign==UNK) // might be a sign change coming
	continue;
      // definately a sign change
      //printf("Sign change in Nright\n");
      res+=(t1_ptr-ptr+1);
      last_sign=this_sign;
    }
  // now look for stationary points indicating missed pairs of zeros
  // assume for safety the zeros are missed just before right hand point
  for(int64_t i=t0_ptr;i<t1_ptr-2;i++)
    {
      if(mpfi_is_pos(f_vec[i])&&mpfi_is_pos(f_vec[i+1])&&mpfi_is_pos(f_vec[i+2])
	 &&(mpfi_cmp(f_vec[i],f_vec[i+1])>0)&&(mpfi_cmp(f_vec[i+1],f_vec[i+2])<0))
	{
	  if(do_stat_point(i,i+1,i+2,f_vec,f_vec_bar,POS,delta,true,NULL))
	    res+=(t1_ptr-i-2)*2;
	  else
	    printf("Failed to resolve stat point in Turing Zone. Continuing.\n");
	  continue;
	}
      if(mpfi_is_neg(f_vec[i])&&mpfi_is_neg(f_vec[i+1])&&mpfi_is_neg(f_vec[i+2])
	 &&(mpfi_cmp(f_vec[i],f_vec[i+1])<0)&&(mpfi_cmp(f_vec[i+1],f_vec[i+2])>0))
	{
	  if(do_stat_point(i,i+1,i+2,f_vec,f_vec_bar,NEG,delta,true,NULL))
	    res+=(t1_ptr-i-2)*2;
	  else
	    printf("Failed to resolve stat point in Turing Zone. Continuing.\n");	}
    }
  return(res*delta);
}

double Nright_int(uint64_t t0_ptr, uint64_t t1_ptr, mpfi_t *lam, mpfi_t *lam_bar, double one_over_A, uint64_t upsam)
{
  //printf("In Nright_int\n");
  double res=0.0;
  uint64_t ptr=t0_ptr*upsam;
  sign_t last_sign=UNK;
  while(last_sign==UNK)
    {
      if(ptr%upsam==0) // on a sample point
	{
	  mpfr_set_ui(pc_lt,ptr/upsam,GMP_RNDN);
	  mpfi_get_fr(pc_lf,lam[ptr/upsam]);
	  last_sign=sign(pc_lf);
	}
      else
	{
	  mpfr_set_d(pc_lt,ptr,GMP_RNDN);
	  mpfr_div_ui(pc_lt,pc_lt,upsam,GMP_RNDN);
	  mpfr_inter(pc_lf,lam,lam_bar,pc_lt,one_over_A);
	  last_sign=sign(pc_lf);
	}
      ptr++;
    }
  if((ptr-1)%upsam!=0) // not on a smaple point, check with mpfi
    {
      mpfi_set_ui(dz_tmp4,ptr-1);
      mpfi_div_ui(dz_tmp3,dz_tmp4,upsam);
      mpfi_inter(dz_tmp4,lam,lam_bar,dz_tmp3,one_over_A);
      if(sign(dz_tmp4)!=last_sign)
	{
	  mpfi_print_str("Mpfi=",dz_tmp4);
	  mpfr_print_str("mpfr=",pc_lf);
	  printf("Error upsampling in Turing Zone. Sign mismatch. Exiting.\n");
	  exit(0);
	}
    }

  while(ptr<t1_ptr*upsam)
    {
      sign_t this_sign;
      if(ptr%upsam==0)
	{
	  mpfr_set_ui(pc_rt,ptr/upsam,GMP_RNDN);
	  mpfi_get_fr(pc_rf,lam[ptr/upsam]);
	  this_sign=sign(pc_rf);
	}
      else
	{
	  mpfr_set_d(pc_rt,ptr,GMP_RNDN);
	  mpfr_div_ui(pc_rt,pc_rt,upsam,GMP_RNDN);
	  mpfr_inter(pc_rf,lam,lam_bar,pc_rt,one_over_A);
	  this_sign=sign(pc_rf);
	}

      if((this_sign&last_sign)==0) // a sign change
	{
	  //printf("Sign change detected before %f\n",(double) ptr/upsam*one_over_A);
	  if((ptr-1)%upsam!=0) // not on a sample point, so check with mpfi
	    {
	      mpfi_set_ui(dz_tmp4,ptr-1);
	      mpfi_div_ui(dz_tmp3,dz_tmp4,upsam);
	      mpfi_inter(dz_tmp4,lam,lam_bar,dz_tmp3,one_over_A);
	      if(sign(dz_tmp4)!=last_sign)
		{
		  printf("Error upsampling in Turing Zone. Sign mismatch. Exiting.\n");
		  exit(0);
		}
	    }
	  /*
	  printf("Zero found between\n");
	  mpfr_print_str("",pc_lt);
	  printf("and\n");
	  mpfr_print_str("",pc_rt);
	  printf("value between\n");
	  mpfr_print_str("",pc_lf);
	  printf("and\n");
	  mpfr_print_str("",pc_rf);
	  */
	  res+=((double) t1_ptr*upsam-ptr+1.0)/upsam;
	  last_sign=this_sign;
	  mpfr_set(pc_lt,pc_rt,GMP_RNDN);
	  mpfr_set(pc_lf,pc_rf,GMP_RNDN);
	  ptr++;
	  continue;
	}
      ptr++;
    }
  return(res*one_over_A);
}

// returns upper limit for int_t0^{t0+h} S(t) dt
// t0>50
// Rumely 
void St_int(mpfi_ptr res, mpfi_ptr t, uint64_t q) 
{
  //assert(mpfi_cmp_d(t,168.0*M_PI+0.1)>0);
  //mpfi_print_str("Doing St_int with t=",t);
  mpfi_mul_d(res,t,(double)q/2.0);
  mpfi_div(res,res,mpfi_pi);
  mpfi_log(res,res);
  mpfi_mul_d(res,res,0.1242001); // ensures > 0.1242 and 1.8397 resp
  mpfi_add_d(res,res,1.8397001);
  //mpfi_neg(tm_tmp,res);
  //mpfi_put(res,tm_tmp);
}

// h+(2ht_0+h^2)/4 log(q/pi)
void ln_term(mpfi_ptr res, mpfi_ptr t0, mpfi_ptr t1, mpfi_ptr h1, uint64_t q)
{
  mpfi_sqr(res,h1); // h^2
  mpfi_mul(ln_tmp,h1,t0); //t0*h
  mpfi_mul_ui(ln_tmp,ln_tmp,2); // 2t0*h
  mpfi_add(res,res,ln_tmp); // 2t0*h+h^2
  mpfi_mul_d(res,res,0.25); // (2t0*h+h^2)/4
  mpfi_set_ui(ln_tmp,q);
  mpfi_div(ln_tmp,ln_tmp,mpfi_pi);
  mpfi_log(ln_tmp,ln_tmp);
  mpfi_mul(res,res,ln_tmp); // (2t0*h+h^2)/4*log(q/pi)
  mpfi_add(res,res,h1);
}


// the maximum number of zeros <=a based on Turing region [a,b]
long int turing_max(mpfi_t *f_vec, mpfi_t *f_vec_bar, long int a_ptr, double a, long int b_ptr, double b, bool even_p,double one_over_A, uint64_t q)
{
  //printf("In turing max from %lu to %lu\n",a_ptr,b_ptr);
  mpfi_set_d(tm_t0,a);
  //printf("t0=");mpfi_printn(tm_t0,10);
  mpfi_set_d(tm_t1,b);
  //printf("t1=");mpfi_printn(tm_t1,10);  
  mpfi_sub(tm_h,tm_t1,tm_t0);
  //  
  St_int(tm_res,tm_t1,q);
  //printf("int S(t)=");mpfi_printn(tm_res,10);
  mpfi_sub_d(tm_res,tm_res,Nright_int(a_ptr,b_ptr,f_vec,f_vec_bar,one_over_A));
  //printf("int S(t) - int Nright");mpfi_printn(tm_res,10);
  ln_term(tm_tmp,tm_t0,tm_t1,tm_h,q);
  //printf("int - t log(pi)/2=");mpfi_printn(tm_tmp,10);
  im_int(tm_tmp1,tm_t0,tm_t1,even_p);
  //printf("Im lnGamma term=");mpfi_printn(tm_tmp1,10);
  mpfi_add(tm_tmp,tm_tmp,tm_tmp1);
  mpfi_div(tm_tmp,tm_tmp,mpfi_pi);
  mpfi_add(tm_res,tm_res,tm_tmp);
  mpfi_div(tm_res,tm_res,tm_h);
  //printf("turing_max computed ");mpfi_printn(tm_res,10);
  mpfi_get_right(tm_mpfr,tm_res);
  return(mpfr_get_si(tm_mpfr,GMP_RNDD)); // note round down here
}
// the maximum number of zeros <=a based on Turing region [a,b]
long int turing_max(mpfi_t *f_vec, mpfi_t *f_vec_bar, long int a_ptr, double a, long int b_ptr, double b, bool even_p,double one_over_A, uint64_t q, uint64_t upsam)
{
  //printf("In turing max from %lu to %lu\n",a_ptr,b_ptr);
  mpfi_set_d(tm_t0,a);
  //printf("t0=");mpfi_printn(tm_t0,10);
  mpfi_set_d(tm_t1,b);
  //printf("t1=");mpfi_printn(tm_t1,10);  
  mpfi_sub(tm_h,tm_t1,tm_t0);
  //  
  St_int(tm_res,tm_t1,q);
  //printf("int S(t)=");mpfi_printn(tm_res,10);
  mpfi_sub_d(tm_res,tm_res,Nright_int(a_ptr,b_ptr,f_vec,f_vec_bar,one_over_A,upsam));
  //printf("int S(t) - int Nright");mpfi_printn(tm_res,10);
  ln_term(tm_tmp,tm_t0,tm_t1,tm_h,q);
  //printf("int - t log(pi)/2=");mpfi_printn(tm_tmp,10);
  im_int(tm_tmp1,tm_t0,tm_t1,even_p);
  //printf("Im lnGamma term=");mpfi_printn(tm_tmp1,10);
  mpfi_add(tm_tmp,tm_tmp,tm_tmp1);
  mpfi_div(tm_tmp,tm_tmp,mpfi_pi);
  mpfi_add(tm_res,tm_res,tm_tmp);
  mpfi_div(tm_res,tm_res,tm_h);
  //printf("turing_max computed ");mpfi_printn(tm_res,10);
  mpfi_get_right(tm_mpfr,tm_res);
  return(mpfr_get_si(tm_mpfr,GMP_RNDD)); // note round down here
}

// do a quick scan to locate sign changes for a real character
uint64_t quick_count_real(mpfi_t *lam, uint64_t num_s)
{
  uint64_t ptr=0,count=0;
  sign_t s=sign(lam[ptr++]);
  while(s==UNK)
    s=sign(lam[ptr++]);
  while(ptr<=num_s)
    {
      sign_t this_sign=sign(lam[ptr++]);
      if(this_sign&s)
	continue;
      count++;
      s=this_sign;
    }
  return(count);
}

// ditto for a pair of conjugate characters
uint64_t quick_count_cmplx(mpfi_t *lam, mpfi_t *lam_bar, uint64_t num_s)
{
  return(quick_count_real(lam,num_s)+quick_count_real(lam_bar,num_s));
}

char signc(mpfi_ptr x)
{
  sign_t s=sign(x);
  if(s==UNK) return('?');
  if(s==POS) return('+');
  return('-');
}

// do a rigorous zero count, isolating them all and writing them out
// we return the number of zeros written
uint64_t proper_count(mpfi_t *lam, mpfi_t *lam_bar, uint64_t end, double one_over_A, FILE *outfile)
{
  // start the difference counter from 0
  mpfr_set_ui(last_t,0,GMP_RNDN);
  uint64_t ptr=0,zero_count=0;
  sign_t last_sign=sign(lam[ptr++]);

  while(last_sign==UNK) last_sign=sign(lam[ptr++]);

  uint64_t last_change=ptr-1;

  while(ptr<end)
    {
      sign_t this_sign=sign(lam[ptr]);
      if((this_sign&last_sign)==0)
	{
	  mpfr_set_ui(pc_lt,last_change,GMP_RNDN);
	  mpfr_set_ui(pc_rt,ptr,GMP_RNDN);
	  mpfi_get_fr(pc_lf,lam[last_change]);
	  mpfi_get_fr(pc_rf,lam[ptr]);
	  do_zero(pc_lt,pc_lf,pc_rt,pc_rf,lam,lam_bar,one_over_A,outfile);
	  zero_count++;
	  last_sign=this_sign;
	  last_change=ptr;
	  ptr++;
	  continue;
	}
      if(this_sign!=UNK) // same sign as before
	last_change=ptr;
      // now check stat points
      if(mpfi_is_pos(lam[ptr-1])&&mpfi_is_pos(lam[ptr])&&mpfi_is_pos(lam[ptr+1])
	 &&(mpfi_cmp(lam[ptr-1],lam[ptr])>0)&&(mpfi_cmp(lam[ptr],lam[ptr+1])<0))
	{
	  if(!do_stat_point(ptr-1,ptr,ptr+1,lam,lam_bar,POS,one_over_A,false,outfile))
	    {
	      printf("Failed to resolve stat point in Main Zone. Exiting.\n");
	      exit(0);
	    }
	  zero_count+=2;
	  ptr++;
	  continue;
	}
      if(mpfi_is_neg(lam[ptr-1])&&mpfi_is_neg(lam[ptr])&&mpfi_is_neg(lam[ptr+1])
	 &&(mpfi_cmp(lam[ptr-1],lam[ptr])<0)&&(mpfi_cmp(lam[ptr],lam[ptr+1])>0))
	{
	  if(!do_stat_point(ptr-1,ptr,ptr+1,lam,lam_bar,NEG,one_over_A,false,outfile))
	    {
	      printf("Failed to resolve stat point in Main Zone. Exiting.\n");	
	      exit(0);
	    }
	  zero_count+=2;
	  ptr++;
	  continue;
	}
      ptr++;
    }
  if((sign(lam[end])&last_sign)==0)
    {
      mpfr_set_ui(pc_lt,last_change,GMP_RNDN);
      mpfr_set_ui(pc_rt,ptr,GMP_RNDN);
      mpfi_get_fr(pc_lf,lam[last_change]);
      mpfi_get_fr(pc_rf,lam[end]);
      do_zero(pc_lt,pc_lf,pc_rt,pc_rf,lam,lam,one_over_A,outfile);
      zero_count++;
    }
  return(zero_count);
}
// do a rigorous zero count with upsampling, isolating them all and writing them out
// we return the number of zeros written
uint64_t proper_count(mpfi_t *lam, mpfi_t *lam_bar, uint64_t end, double one_over_A, FILE *outfile, uint64_t upsam)
{
  // start the difference counter from 0
  mpfr_set_ui(last_t,0,GMP_RNDN);
  uint64_t ptr=0,zero_count=0;
  sign_t last_sign=UNK;

  // pc_lt is the ptr to the last change of sign
  // pc_lf is the value of f there
  // pc_rt is the current ptr
  // pc_rf is the value here

  while(last_sign==UNK)
    {
      if(ptr%upsam==0)
	{
	  mpfr_set_ui(pc_lt,ptr/upsam,GMP_RNDN);
	  mpfi_get_fr(pc_lf,lam[ptr/upsam]);
	  last_sign=sign(pc_lf);
	}
      else
	{
	  mpfr_set_d(pc_lt,ptr,GMP_RNDN);
	  mpfr_div_ui(pc_lt,pc_lt,upsam,GMP_RNDN);
	  mpfr_inter(pc_lf,lam,lam_bar,pc_lt,one_over_A);
	  last_sign=sign(pc_lf);
	}
      ptr++;
    }

  while(ptr<end*upsam)
    {

      sign_t this_sign;
      if(ptr%upsam==0)
	{
	  mpfr_set_ui(pc_rt,ptr/upsam,GMP_RNDN);
	  mpfi_get_fr(pc_rf,lam[ptr/upsam]);
	  this_sign=sign(pc_rf);
	}
      else
	{
	  mpfr_set_d(pc_rt,ptr,GMP_RNDN);
	  mpfr_div_ui(pc_rt,pc_rt,upsam,GMP_RNDN);
	  mpfr_inter(pc_rf,lam,lam_bar,pc_rt,one_over_A);
	  this_sign=sign(pc_rf);
	}

      if((this_sign&last_sign)==0) // a sign change
	{
	  /*
	  printf("Zero found between\n");
	  mpfr_print_str("",pc_lt);
	  printf("and\n");
	  mpfr_print_str("",pc_rt);
	  printf("value between\n");
	  mpfr_print_str("",pc_lf);
	  printf("and\n");
	  mpfr_print_str("",pc_rf);
	  */
	  do_zero(pc_lt,pc_lf,pc_rt,pc_rf,lam,lam_bar,one_over_A,outfile);
	  zero_count++;
	  last_sign=this_sign;
	  mpfr_set(pc_lt,pc_rt,GMP_RNDN);
	  mpfr_set(pc_lf,pc_rf,GMP_RNDN);
	  ptr++;
	  continue;
	}
      if(this_sign!=UNK) // same sign as before
	{
	  mpfr_set_ui(pc_lt,ptr,GMP_RNDN);
	  mpfr_div_ui(pc_lt,pc_lt,upsam,GMP_RNDN);
	  mpfr_set(pc_lf,pc_rf,GMP_RNDN);
	}
      ptr++;
    }

  if((sign(lam[end])&last_sign)==0)
    {
      mpfr_set_ui(pc_rt,end,GMP_RNDN);
      mpfi_get_fr(pc_rf,lam[end]);
      do_zero(pc_lt,pc_lf,pc_rt,pc_rf,lam,lam,one_over_A,outfile);
      zero_count++;
    }
  return(zero_count);
}


// find the zeros for a real character
void find_zeros_real(mpfi_t *lam, uint64_t num_s, double del_t, bool even_p, uint64_t q, FILE *outfile,uint64_t index, uint64_t upsam, FILE *zinfile)
{
  uint64_t zz;
  if(fread(&zz,sizeof(uint64_t),1,zinfile)!=1)
    {
      printf("Error reading zeros count from zeros file. Exiting.\n");
      exit(0);
    }
  if(zz!=0) // worked at previous upsampling rate
    {
      fwrite(&zz,sizeof(uint64_t),1,outfile);
      uint64_t a; uint32_t b; uint8_t c;
      for(uint64_t z=0;z<zz;z++) // no null record here
	{
	  if(fread(&a,sizeof(uint64_t),1,zinfile)!=1)
	    {
	      printf("Error reading zero from zeros file. Exiting.\n");
	      exit(0);
	    }
	  if(fwrite(&a,sizeof(uint64_t),1,outfile)!=1)
	    {
	      printf("Error writing zero to out file. Exiting.\n");
	      exit(0);
	    }
	  if(fread(&b,sizeof(uint32_t),1,zinfile)!=1)
	    {
	      printf("Error reading zero from zeros file. Exiting.\n");
	      exit(0);
	    }
	  if(fwrite(&b,sizeof(uint32_t),1,outfile)!=1)
	    {
	      printf("Error writing zero to out file. Exiting.\n");
	      exit(0);
	    }
	  if(fread(&c,sizeof(uint8_t),1,zinfile)!=1)
	    {
	      printf("Error reading zero from zeros file. Exiting.\n");
	      exit(0);
	    }
	  if(fwrite(&c,sizeof(uint8_t),1,outfile)!=1)
	    {
	      printf("Error writing zero to out file. Exiting.\n");
	      exit(0);
	    }
	}
      return;
    }

  uint64_t a=T0/del_t,b=T1/del_t;

  uint64_t turing_mx=turing_max(lam,lam,a,T0,b,T1,even_p,del_t,q);

  // we know the parity of turing_mx, so use it
  sign_t ls=sign(lam[0]);
  uint64_t ptr=1;
  while(ls==UNK)
    ls=sign(lam[ptr++]);
  sign_t ts=sign(lam[a]);
  ptr=a-1;
  while(ts==UNK)
    ts=sign(lam[ptr--]);
  if(ls&ts) // Lambda starts and ends with same sign
    {
      if(turing_mx&1)
	turing_mx--;
    }
  else // different signs
    {
      if(!(turing_mx&1))
	turing_mx--;
    }

  uint64_t qc=quick_count_real(lam,a);
  if(qc!=turing_mx) // need some more zeros, so try stat points
    qc+=num_stat_pts(lam,0,a);
  if(qc==turing_mx) // all zeros accounted for w/o upsampling
    {
      //printf("Writing %lu zeros\n",qc);
      fwrite(&qc,sizeof(uint64_t),1,outfile);
      uint64_t pc=proper_count(lam,lam,a,del_t,outfile);
      if(pc!=qc)
	{
	  printf("Index=%lu pc=%lu qc=%lu\n",index,pc,qc);
	  printf("Wierd error counting zeros. Exiting.\n");
	  exit(0);
	}
      return;
    }
  uint64_t pc=qc;
  if(upsam!=1) // upsampling
    {
      fwrite(&turing_mx,sizeof(uint64_t),1,outfile);
      pc=proper_count(lam,lam,a,del_t,outfile,upsam);
      if(pc==turing_mx)
	return;
    }
  printf("Expecting %lu zeros.\n",turing_mx);
  printf("Found %lu zeros.\n",pc);
  printf("Signs go from %c to %c\n",signc(lam[0]),signc(lam[a]));
  if(sign(lam[0])&sign(lam[a]))
    printf("Should be even.\n");
  else
    printf("Should be odd.\n");
  uint64_t tz=num_stat_pts(lam,a,b);
  if(tz)
    printf("We missed %lu zeros in Turing Zone.\n",tz);
  if(upsam==1)
    {
      qc=0;
      fwrite(&qc,sizeof(uint64_t),1,outfile);      	
    }
  else
    exit(0);
}

void find_zeros_cmplx(mpfi_t *lam, mpfi_t *lam_bar, uint64_t num_s, double del_t,bool even_p, uint64_t q, FILE *outfile, uint64_t index, uint64_t upsam, FILE *zinfile)
{
  //printf("In find_zeros_cmplx.\n");
  uint64_t zz;
  if(fread(&zz,sizeof(uint64_t),1,zinfile)!=1)
    {
      printf("Error reading zeros count from zeros file. Exiting.\n");
      exit(0);
    }
  if(zz!=0) // worked at previous upsampling rate
    {
      fwrite(&zz,sizeof(uint64_t),1,outfile);
      uint64_t a; uint32_t b; uint8_t c;
      for(uint64_t z=0;z<=zz;z++) // include null record between characters
	{
	  if(fread(&a,sizeof(uint64_t),1,zinfile)!=1)
	    {
	      printf("Error reading zero from zeros file. Exiting.\n");
	      exit(0);
	    }
	  if(fwrite(&a,sizeof(uint64_t),1,outfile)!=1)
	    {
	      printf("Error writing zero to out file. Exiting.\n");
	      exit(0);
	    }
	  if(fread(&b,sizeof(uint32_t),1,zinfile)!=1)
	    {
	      printf("Error reading zero from zeros file. Exiting.\n");
	      exit(0);
	    }
	  if(fwrite(&b,sizeof(uint32_t),1,outfile)!=1)
	    {
	      printf("Error writing zero to out file. Exiting.\n");
	      exit(0);
	    }
	  if(fread(&c,sizeof(uint8_t),1,zinfile)!=1)
	    {
	      printf("Error reading zero from zeros file. Exiting.\n");
	      exit(0);
	    }
	  if(fwrite(&c,sizeof(uint8_t),1,outfile)!=1)
	    {
	      printf("Error writing zero to out file. Exiting.\n");
	      exit(0);
	    }
	}
      return;
    }

  uint64_t a=T0/del_t,b=T1/del_t;
  uint64_t turing_mx=turing_max(lam,lam_bar,a,T0,b,T1,even_p,del_t,q,upsam)+
    turing_max(lam_bar,lam,a,T0,b,T1,even_p,del_t,q,upsam);
  //printf("Turing max=%lu\n",turing_mx);
  sign_t ls=sign(lam[a]);
  uint64_t ptr=a-1;
  while(ls==UNK)
    ls=sign(lam[ptr--]);
  sign_t ts=sign(lam_bar[a]);
  ptr=a-1;
  while(ts==UNK)
    ts=sign(lam_bar[ptr--]);
  if(ls&ts) // signs are same, parity even
    {
      if(turing_mx&1)
	turing_mx--;
    }
  else // parity odd
    {
      if(!(turing_mx&1))
	turing_mx--;
    }
  uint64_t qc=quick_count_cmplx(lam,lam_bar,a),pc;
  //printf("quick_count found %lu\n",qc);
  if(qc!=turing_mx)
    qc+=num_stat_pts(lam,0,a);
  //printf("quick_count now %lu\n",qc);
  if(qc!=turing_mx)
    qc+=num_stat_pts(lam_bar,0,a);
  //printf("quick_count now %lu\n",qc);
  
  if(qc==turing_mx)
    {
      //printf("Writing %lu zeros\n",qc);
      fwrite(&qc,sizeof(uint64_t),1,outfile);
      pc=proper_count(lam,lam_bar,a,del_t,outfile);
      //printf("pc=%lu\n",pc);
      write_null_record(outfile);
      pc+=proper_count(lam_bar,lam,a,del_t,outfile);
      //printf("pc=%lu\n",pc);
      if(pc!=qc)
	{
	  printf("Index=%lu pc=%lu qc=%lu\n",index,pc,qc);
	  printf("Wierd error counting zeros. Exiting.\n");
	  exit(0);
	}
      return;
    }
  if(upsam!=1) // upsampling
    {
      fwrite(&turing_mx,sizeof(uint64_t),1,outfile);
      uint64_t pc=proper_count(lam,lam_bar,a,del_t,outfile,upsam);
      write_null_record(outfile);
      pc+=proper_count(lam_bar,lam,a,del_t,outfile,upsam);
      if(pc==turing_mx)
	return;
    }
  printf("Expecting at most %lu zeros.\n",turing_mx);
  printf("Found %lu zeros.\n",pc);
  if((sign(lam[0])==UNK)||(sign(lam_bar[0])==UNK))
    printf("Unkown central point.\n");
  printf("Signs go from %c to %c\n",signc(lam[a]),signc(lam_bar[a]));
  if(sign(lam_bar[a])&sign(lam[a]))
    printf("Should be even.\n");
  else
    printf("Should be odd.\n");
  uint64_t tz=num_stat_pts(lam,a,b)+num_stat_pts(lam_bar,a,b);
  if(tz)
    printf("We missed %lu zeros in Turing Zone.\n",tz);
  if(upsam==1)
    {
      qc=0;
      fwrite(&qc,sizeof(uint64_t),1,outfile);      	
    }
  else // we are upsampling and have written an incorrect zero count so give up
    exit(0);
}

void read_real(mpfi_ptr res,FILE *infile,mpz_ptr ztemp, mpfr_ptr rtemp)
{ 
  if(!mpz_inp_raw(ztemp,infile))
    {
      printf("Error reading left mantissa. Exiting.\n");
      exit(0);
    }
  mpfr_exp_t ex;
  if(fread(&ex,sizeof(mpfr_exp_t),1,infile)!=1)
    {
      printf("Error reading left exponent. Exiting.\n");
      exit(0);
    }
  mpfr_set_z_2exp(rtemp,ztemp,ex,GMP_RNDD);
  mpfi_set_fr(res,rtemp);

  if(!mpz_inp_raw(ztemp,infile))
    {
      printf("Error reading right mantissa. Exiting.\n");
      exit(0);
    }
  if(fread(&ex,sizeof(mpfr_exp_t),1,infile)!=1)
    {
      printf("Error reading right exponent. Exiting.\n");
      exit(0);
    }
  mpfr_set_z_2exp(rtemp,ztemp,ex,GMP_RNDU);
  mpfi_put_fr(res,rtemp);
}

int main(int argc, char **argv)
{
  if(argc!=7)
    {
      printf("Usage:- %s <specfile> <this_proc> <num procs> <zeros infile> <zeros outfile> <upsam rate>.\n",argv[0]);
      exit(0);
    }

  uint64_t this_proc=atol(argv[2]);
  uint64_t num_procs=atol(argv[3]);
  uint64_t upsam=atol(argv[6]);
  this_proc=this_proc%num_procs;

  dft_init(PREC);
  FILE *specfile=fopen(argv[1],"r");

  uint64_t q,zq;
  double *t_starts,del_t;
  if(!specfile)
    {
      printf("Failed to open %s for input. Exiting.\n",argv[1]);
      exit(0);
    }

  uint64_t num_files;
  if(fscanf(specfile,"%lu\n",&num_files)!=1)
    {
      printf("Error reading number of files from spec file. Exiting.\n");
      exit(0);
    }
  FILE *zinfile=fopen(argv[4],"rb");
  if(!zinfile)
    {
      printf("Failed to open %s for input. Exiting.\n",argv[4]);
      exit(0);
    }


  FILE *outfile=fopen(argv[5],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for output. Exiting.\n",argv[5]);
      exit(0);
    }


  t_starts=(double *)malloc(sizeof(double)*(num_files+1));
  FILE **infiles=(FILE **)malloc(sizeof(FILE *)*num_files);
  if(!t_starts||!infiles)
    {
      printf("Error allocating memory for files and/or t_starts. Exiting.\n");
      exit(0);
    }

  for(uint64_t f=0;f<num_files;f++)
    {
      char fname[1024];
      if(!fscanf(specfile,"%s\n",fname))
	{
	  printf("Error reading filename from spec file. Exiting.\n");
	  exit(0);
	}
      infiles[f]=fopen(fname,"rb");
      if(!infiles[f])
	{
	  printf("Error opening file %s for binary input. Exiting.\n",fname);
	  exit(0);
	}
    }

  if(fread(&q,sizeof(uint64_t),1,infiles[0])!=1)
    {
      printf("Error reading q from file 0. Exiting.\n");
      exit(0);
    }
  if(q>MAX_Q)
    {
      printf("Interpolation error computed for q<=%lu. Exiting.\n",MAX_Q);
      exit(0);
    }
  if(fread(&zq,sizeof(uint64_t),1,zinfile)!=1)
    {
      printf("Error reading q from zeros file. Exiting.\n");
      exit(0);
    }
  if(zq!=q)
    {
      printf("Mismatch between q in zeros file and q in l_func files. Exiting.\n");
      exit(0);
    }
  fwrite(&q,sizeof(uint64_t),1,outfile);
  for(uint64_t f=1,qc;f<num_files;f++)
    if((fread(&qc,sizeof(uint64_t),1,infiles[f])!=1)||qc!=q)
    {
      printf("Error reading q from file %lu. Exiting.\n",f);
      exit(0);
    }


  if((fread(t_starts,sizeof(double),1,infiles[0])!=1)||t_starts[0]!=0.0)
    {
      printf("Error reading t_start=0.0 from file 0. Exiting.\n");
      exit(0);
    }
  if(fread(t_starts+1,sizeof(double),1,infiles[0])!=1)
    {
      printf("Error reading t_end. Exiting.\n");
      exit(0);
    }
  if(fread(&del_t,sizeof(double),1,infiles[0])!=1)
    {
      printf("Error reading del_t. Exiting.\n");
      exit(0);
    }
  if(del_t!=SAMPLE_RATE)
    {
      printf("Interpolation error computed for sample rate=%e. Exiting.\n",SAMPLE_RATE);
      exit(0);
    }

  for(uint64_t f=1;f<num_files;f++)
    {
      double ts;
      if((fread(&ts,sizeof(double),1,infiles[f])!=1)||(ts!=t_starts[f]))
	{
	  printf("Error reading t_start from file %lu. Exiting\n",f);
	  exit(0);
	}
      if(fread(t_starts+f+1,sizeof(double),1,infiles[f])!=1)
	{
	  printf("Error reading t_end from file %lu. Exiting.\n",f);
	  exit(0);
	}
      if((fread(&ts,sizeof(double),1,infiles[f])!=1)||(ts!=del_t))
	{
	  printf("Error reading del_t from file %lu. Exiting\n",f);
	  exit(0);
	}
    }

  if(t_starts[num_files]!=T2)
    {
      printf("We need samples to go to %f. Exiting.\n",T2);
      exit(0);
    }

  mpfi_t *lam,*lam_bar,x;
  mpfi_init(x);
  uint64_t num_s=t_starts[num_files]/del_t;
  lam=(mpfi_t *)malloc(sizeof(mpfi_t)*num_s);
  lam_bar=(mpfi_t *)malloc(sizeof(mpfi_t)*num_s);
  if(!lam||!lam_bar)
    {
      printf("Error allocating memory for lam/lam_bar. Exiting.\n");
      exit(0);
    }

  for(uint64_t s=0;s<num_s;s++)
    {
      mpfi_init(lam[s]);
      mpfi_init(lam_bar[s]);
    }

  //
  // Warning!!!!
  // DirichletGroup G(q) appears to trample on memory
  // So I set up all my variables and constants later.
  DirichletGroup G(q);
  turing_setup();
  inter_setup(del_t);
  init_out_bytes();

  //mpfi_print_str("inter error=",mpfi_inter_err);
  mpfr_t rtemp;
  mpfr_init(rtemp);
  mpz_t ztemp;
  mpz_init(ztemp);

  uint64_t index,index_counter=0,zindex;
  while(fread(&index,sizeof(uint64_t),1,infiles[0])==1)
    {
      bool even_p=G.character(index).is_even(); // this is nasty
      uint64_t i1;
      for(uint64_t f=1;f<num_files;f++)
	{
	  if(fread(&i1,sizeof(uint64_t),1,infiles[f])!=1)
	    {
	      printf("Error reading index from file %lu. Exiting.\n",f);
	      exit(0);
	    }
	  if(index!=i1)
	    {
	      printf("Mismatch error reading index from file %lu. Read %lu, expected %lu. Exiting.\n",f,i1,index);
	      exit(0);
	    }
	}
      for(uint64_t f=0,s=0;f<num_files;f++)
	{
	  for(double t=t_starts[f];t<t_starts[f+1];t+=del_t)
	    {
	      read_real(lam[s++],infiles[f],ztemp,rtemp);
	      //printf("%lu %f %f\n",i,t,mpfi_get_d(x));
	    }
	}

      uint64_t inverse=InvMod(index,q);
      if(inverse==index) // ths character is real, so sort it
	{
	  if((index_counter%num_procs)==this_proc)
	    {
	      printf("Doing real index %lu\n",index);
	      if(fread(&zindex,sizeof(uint64_t),1,zinfile)!=1)
		{
		  printf("Error reading index from zeros file. Exiting.\n");
		  exit(0);
		}
	      if(zindex!=index)
		{
		  printf("Incorrect index from zeros file. Exiting.\n");
		  exit(0);
		}
	      fwrite(&index,sizeof(uint64_t),1,outfile);
	      find_zeros_real(lam,num_s,del_t,even_p,q,outfile,index,upsam,zinfile);
	    }
	  index_counter++;
	}
      else // read its conjugate
	{
	  //printf("Looking for complex conjugate %lu\n",inverse);
	  for(uint64_t f=0;f<num_files;f++)
	    {
	      if(fread(&i1,sizeof(uint64_t),1,infiles[f])!=1)
		{
		  printf("Error reading index from file %lu. Exiting.\n",f);
		  exit(0);
		}
	      if(inverse!=i1)
		{
		  printf("Mismatch error reading index from file %lu. Read %lu, expected %lu. Exiting.\n",f,i1,inverse);
		  exit(0);
		}
	    }
	  for(uint64_t f=0,s=0;f<num_files;f++)
	    {
	      for(double t=t_starts[f];t<t_starts[f+1];t+=del_t)
		{
		  read_real(lam_bar[s++],infiles[f],ztemp,rtemp);
		  //printf("%lu %f %f\n",i,t,mpfi_get_d(x));
		}
	    }
	  mpfi_sub(x,lam[0],lam_bar[0]);
	  if(!mpfi_contains_zero(x))
	    {
	      printf("Character and its conjugate do not agree at t=0. Exiting.\n");
	      exit(0);
	    }
	  if((index_counter%num_procs)==this_proc)
	    {
	      printf("Doing complex index %lu inverse %lu.\n",index,inverse);
	      if(fread(&zindex,sizeof(uint64_t),1,zinfile)!=1)
		{
		  printf("Error reading index from zeros file. Exiting.\n");
		  exit(0);
		}
	      if(zindex!=index)
		{
		  printf("Incorrect index from zeros file. Exiting.\n");
		  exit(0);
		}
	      fwrite(&index,sizeof(uint64_t),1,outfile);
	      find_zeros_cmplx(lam,lam_bar,num_s,del_t,even_p,q,outfile,index,upsam,zinfile);
	    }
	  index_counter++;
	}
    }
  printf("%lu character pairs read.\n",index_counter);
  return(0);
}
