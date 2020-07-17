//
// win_zeta1.8.c
//
// Windowed FFT based Zeta-function calculator
// See Booker - Artin, Turing ....
//
// Version 1.0 Initial implementation
//         1.1 Vary A (and N) between g/f
//         1.2 Use Hermitian symmetry to speed up final iDFT
//         1.3 I forget
//         1.4 Much improved Stationary Point Handling
//         1.5 Removed (most of) spurious sin's when interpolating
//         1.6 Added Binary chop for convergence failures
//         1.7 Changed command line parameters
//         1.8 Improved gamma_factor
//             Fixed Turing estimate (add 1 to allow for arbitary
//             choice of epsilon from +/- 1)
//
// Created: 28 September 2010
// Last Modified: 7 February 2011
//
// DJ Platt
// University of Bristol
//

#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "time.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "../includes/mpfi_fft.h"
#include "../includes/fft_defs.h"
#include "win_zeta.h"


#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define bool int
int HI_PREC;

// parameters for t0<=3*10^10
#define T0_MAX (3.0e10) // maximum t0
#define EXTRA_PREC (50) // more than log2(t0/(2Pi)*log(Msqrt(Pi)))
#define T0_MIN (5000.0) // mimimal input t0
#define UPSAM (32) // the rate of upsampling from g to f
#define N ((int) 1<<20) // final FFT length
#define N1 ((int) N/UPSAM) // intermediate FFT length
#define one_over_A ((double) 21.0/4096.0) // output spacing
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define h ((double) 176431.0/2048.0) // Gaussian width for f(t)
#define B ((double) N*one_over_A)
#define M ((int) 103000) // number of terms in F(x) sum

#define K ((int) 44) // number of Taylor terms
#define gtwiderr_d ((double) 3.2e-82)
#define Gtwiderr_d ((double) 2.3e-213)
#define fhatsumerr_d ((double) 1.7e-83)
#define tayerr_d ((double) 1.5e-82)
#define fhattwiderr_d ((double) 1.0e-307) // actually much smaller
#define ftwiderr_d ((double) 8.1e-211)
#define Fmaxerr_d ((double) 1.0e-307) // Max mod of F(N1/(2B)) actually much smaller

// Turing method parameters
#define TURING_WIDTH (42)
#define TURING_LEN ((int) (TURING_WIDTH/one_over_A))

// Upsampling parameters
#define Ns (70) // take this many points either side of t0
#define H ((double) 2089.0/16384.0) // Gaussian width for upsampling
#define INTER_SPACING (5) // Use every 5th sample
#define intererr_d ((double) 5.0e-41)
#define INTER_A ((double) INTER_SPACING*one_over_A)

// Newton Iteration Parameters
#define NEWTON_ITS (8) // Do this many iterations with Newton Raphson to locate roots
//
// s_vec[n-1] contains log(n sqrt(Pi))/2Pi)-u_m n=1..M
// sk_vec[n-1] contains (log(n sqrt(Pi))/2Pi)-u_m)^k
// buck[n] tells us which sk to put s_vec[n] into
// n_vec[n-1] contains 1/sqrt(n)*(nsqrt(Pi))^(-it0) n=1..M
// g_vec contains g(t)*(2*Pi*I*t)^k/k!
// two_pi_t contains (-2*Pi*t)
// G_vec contains (-2*Pi*i*t)^k*g(t)/k!, then G^(k)(x)/k!, then sum m=1..M G^(k)(x+u_m)S^(k)_m/k!
// f_vec contains f^(x), then f(t)
// sqrts[i] contains (i+1)^-0.5
//

mpfi_t A,A1,two_pi_t[N1],s_vec[M],sk_vec[M],exps[N1],sqrts[M],logisqrtpi[M];
mpfi_c_t g_vec[N1],G_vec[N1],f_vec[N],skn_vec[N1],ws_r[N/4],ws1_r[N1/2],ws1_f[N1/2],n_vec[M];
mpfi_c_t gtwiderr,Gtwiderr,c_tmp,fhatsumerr,fhattwiderr,ftwiderr,tayerr,Fmaxerr;
int buck[M]; // which bucket does log(nsqrt(Pi)/2pi) go in.
mpfi_t *exps1,intererr;
mpfr_t mpfr_pi,inter_tmp,inter_tmp1,inter_tmp2,ip_tmp1,ip_tmp2,ip_tmp3,ip_tmp4,sinc_tmp;
mpfr_t mpfr_two_pi_B,newton_df,two_101,offset,last_t;
mpfi_t msinc_tmp,mpfi_two_pi_B,mip_tmp,minter_tmp,minter_tmp1;
mpfr_t znew_guess,zold_guess,ztn,fl,fm,fr,tl,tm,tr,msin,mcos,fmid;
mpfi_t zlow,zhigh,ztmp,misin;

// note (sign&&sign)==0 <=> signs different and known 

void print_usage()
{
  printf("Usage:- win_zeta1.9 prec t_min num_its step outfile stat_pts?.\n");
  printf("finds zeros in [t_min,t_min+step]*n for n=0..num_its-1.\n");
  exit(1);
}


void set_err(mpfi_c_ptr z, double d)
{
  mpfi_c_init(z);
  mpfi_set_d(z->re,d);
  mpfi_neg(z->im,z->re);
  mpfi_put(z->re,z->im);
  mpfi_set(z->im,z->re);
}

// reverse skn_vec
inline int conj_order(int i)
{
  return(N1-i);
}

// on entry s_vec[i] contains log((i+1)Pi)/2Pi
// on exit                    "" - u_m
void calc_buck()
{
  int i;
  double x;
  mpfi_t tmp;
  mpfi_init(tmp);
  //printf("s_vec[0]=");mpfi_printn(s_vec[0],10);
  for(i=0;i<M;i++)
    {
      buck[i]=mpfi_get_d(s_vec[i])*B+0.5;
      mpfi_set_ui(tmp,buck[i]);
      mpfi_div_d(tmp,tmp,B);
      //printf("n=%d goes in bucket %d\n",i+1,buck[i]);
      mpfi_sub(s_vec[i],s_vec[i],tmp);
    }
  //printf("Before shuffle, n=1 term is in bucket %d\n",buck[0]);
  assert(buck[0]>0); // keep conj_order simple
  assert(buck[M-1]<N1/2); // ditto
  for(i=0;i<M;i++)
    buck[i]=conj_order(buck[i]);
  //printf("buck[0]=%d\nbuck[M-1]=%d\n",buck[0],buck[M-1]);
  //printf("n=1 term is in bucket %d\n",buck[0]);
  mpfi_clear(tmp);
}

mpfi_t atan1_tmp1,atan1_tmp2,im_1,im_2,im_t,im_err;
mpfi_t tm_h,tm_res,tm_tmp,tm_tmp1,im_t0,im_t1;
mpfr_t tm_mpfr,sp_t,sp_ft,sp_t1,sp_ft1;
mpfi_t tm_t0,tm_t1,init_tmp,init_tmp1,sqrt_pi,pi_hi,g_tmp1,g_re,g_im,g_abs,half_log_2_pi,m1,m2;
mpfi_c_t g_z,g_z2;
// approximation to atan for large t
// = pi/2-1/t+[-1/3t^2,1/3t^2]
void mpfi_atan1(mpfi_ptr res, mpfi_ptr t)
{
  mpfi_set_ui(res,1);
  mpfi_div(res,res,t);
  mpfi_sqr(atan1_tmp1,res);
  mpfi_mul(atan1_tmp1,atan1_tmp1,res);
  mpfi_div_ui(atan1_tmp1,atan1_tmp1,3);
  mpfi_neg(atan1_tmp2,atan1_tmp1);
  mpfi_put(atan1_tmp1,atan1_tmp2);
  mpfi_sub(res,mpfi_pi_by_2,res);
  mpfi_add(res,res,atan1_tmp1);
}

mpfi_sin_cos_hi(mpfi_ptr sin_res, mpfi_ptr cos_res, mpfi_ptr t)
{
  mpfi_sin(init_tmp,t);
  mpfi_set(sin_res,init_tmp);
  mpfi_cos(init_tmp,t);
  mpfi_set(cos_res,init_tmp);
}

mpfi_t m1,m2;
#define hi_mul(a,b){mpfi_mul(m1,a->re,b->re);mpfi_mul(m2,a->im,b->im);mpfi_sub(m1,m1,m2);mpfi_mul(m2,a->re,b->im);mpfi_set(a->re,m1);mpfi_mul(m1,a->im,b->re);mpfi_add(a->im,m1,m2);}

#define MAX_LNG_TERMS (17)
int LNG_TERMS;

void set_LNG_TERMS(double t0)
{
  double tmin=(t0-B/2)/2.0;
  if(tmin<1500.0)
    {
      printf("LnGamma is probably not accurate enough this low down. Exiting.\n");
      exit(0);
    }
  if(tmin>9.0e8)
      LNG_TERMS=6;
  else if(tmin>5.0e7)
    LNG_TERMS=7;
  else if(tmin>5.0e6)
    LNG_TERMS=8;
  else if(tmin>9.0e9)
    LNG_TERMS=9;
  else if(tmin>3.0e5)
    LNG_TERMS=10;
  else if(tmin>9.0e4)
    LNG_TERMS=11;
  else if(tmin>4.0e4)
    LNG_TERMS=12;
  else if(tmin>1.0e4)
    LNG_TERMS=14;
  else if(tmin>4.0e3)
    LNG_TERMS=16;
  else if(tmin>2.0e3)
    LNG_TERMS=18;
  else LNG_TERMS=MAX_LNG_TERMS;
}

mpfi_t lng_terms[MAX_LNG_TERMS],lng_err;

double ber_nums[MAX_LNG_TERMS]={1,-1,1,-1,5,-691,7,-3617,43867,-174611,854513,-236364091,8553103,-23749461029.0,8615841276005.0,-7709321041217.0,2577687858367};
unsigned int ber_dens[MAX_LNG_TERMS]={6,30,42,30,66,2730,6,510,798,330,138,2730,6,870,14322,510,6};

void init_lng()
{
  int i;
  for(i=0;i<MAX_LNG_TERMS;i++)
    {
      //printf("%20.18e/%ld\n",ber_nums[i],ber_dens[i]);
      mpfi_init2(lng_terms[i],HI_PREC);
      mpfi_set_d(lng_terms[i],ber_nums[i]);
      mpfi_div_ui(lng_terms[i],lng_terms[i],ber_dens[i]);
      mpfi_div_ui(lng_terms[i],lng_terms[i],(i+i+2)*(i+i+1));
    }
  mpfi_init2(lng_err,HI_PREC);
}

// gamma(1/4+i(t+t0)/2)*exp(Pi*(t+t0)/4)*exp(-t^2/(2h^2))
void gamma_factor(mpfi_c_ptr res, mpfi_ptr gaussian, double t)
{
  int i;
  // set g_re to Re (lngamma(1/4+I*t/2)+Pi*t/4-gaussian)
  // set g_im to Im (lngamma(1/4+I*t/2)

  mpfi_set_d(g_tmp1,0.5);
  mpfi_div_d(g_tmp1,g_tmp1,t);
  mpfi_atan(g_tmp1,g_tmp1);
  mpfi_div_2ui(g_im,g_tmp1,2); // 1/4*atan(1/(2t))
  mpfi_mul_d(g_re,g_tmp1,t/2.0); // t/2*atan(1/(2t))
  mpfi_add(g_re,g_re,half_log_2_pi); // +1/2 log(2 Pi) -1/4
  mpfi_set_d(g_abs,t/2.0);
  mpfi_sqr(g_abs,g_abs);
  mpfi_add_d(g_abs,g_abs,1.0/16.0); // 1/16+t^2/4
  mpfi_log(g_tmp1,g_abs);
  mpfi_div_2ui(g_tmp1,g_tmp1,3); // log(~)/8
  mpfi_sub(g_re,g_re,g_tmp1);
  mpfi_mul_d(g_tmp1,g_tmp1,2.0*t); // t/4 log(~)
  mpfi_add(g_im,g_im,g_tmp1);
  mpfi_sub_d(g_im,g_im,t/2.0);
  mpfi_div_2ui(g_tmp1,pi_hi,3);
  mpfi_sub(g_im,g_im,g_tmp1); // -pi/8

  // now do the sum of bernoulli bits
  mpfi_set_d(g_z->re,0.25);
  mpfi_set_d(g_z->im,-t/2.0);
  mpfi_div(g_z->re,g_z->re,g_abs);
  mpfi_div(g_z->im,g_z->im,g_abs); // 1/z
  
  //printf("t=%13.11e\n1/z=",t);
  //mpfi_c_printn(g_z,100);
  mpfi_c_set(g_z2,g_z);
  hi_mul(g_z2,g_z); // 1/z^2
  //printf("1/z^2=");mpfi_c_printn(g_z2,100);

  for(i=0;i<LNG_TERMS-1;i++)
    {
      mpfi_mul(g_tmp1,g_z->re,lng_terms[i]);
      mpfi_add(g_re,g_re,g_tmp1);
      mpfi_mul(g_tmp1,g_z->im,lng_terms[i]);
      mpfi_add(g_im,g_im,g_tmp1);
      hi_mul(g_z,g_z2);
    }
  mpfi_mul(g_tmp1,g_z->re,lng_terms[i]);
  mpfi_add(g_re,g_re,g_tmp1);
  mpfi_mul(g_tmp1,g_z->im,lng_terms[i]);
  mpfi_add(g_im,g_im,g_tmp1);

  mpfi_c_abs(lng_err,g_z);
  mpfi_mul_2ui(lng_err,lng_err,LNG_TERMS);
  mpfi_mul(lng_err,lng_err,lng_terms[i]);
  mpfi_neg(g_tmp1,lng_err);
  mpfi_put(lng_err,g_tmp1);
  mpfi_add(g_re,g_re,lng_err);
  mpfi_add(g_im,g_im,lng_err);

  //printf("Re lngamma+pit/4=");mpfi_printn(g_re,120);
  //printf("Im lngamma+pit/4=");mpfi_printn(g_im,120);
  // divide by gaussian exp(t^2/(2h^2))
  //printf("gaussian=");mpfi_printn(gaussian,120);
  mpfi_sub(g_re,g_re,gaussian);
  //printf("Re lngamma+pit/4-t^2/(2h^2)=");mpfi_printn(g_re,120);

  // exponentiate and return
  mpfi_exp(g_tmp1,g_re); //exp(x)
  mpfi_cos(g_re,g_im); // cos(y)
  mpfi_mul(g_re,g_re,g_tmp1); // re=exp(x)cos(y)
  mpfi_sin(g_im,g_im); // sin(y)
  mpfi_mul(g_im,g_im,g_tmp1); // sin(y)*exp(x)
  mpfi_set(res->re,g_re);
  mpfi_set(res->im,g_im);
  //printf("exp()=");mpfi_c_printn(res,100);
  //exit(0);
}
 
      
// on exit g[n]=g((n-N/2)/A;0)
void init(double t0)
{
  int i;
  double t;
  mpfr_init(mpfr_pi);
  mpfi_get_fr(mpfr_pi,mpfi_pi);
  mpfr_init(inter_tmp);
  mpfr_init(inter_tmp1);
  mpfr_init(inter_tmp2);
  mpfr_init(ip_tmp2);
  mpfr_init(ip_tmp1);
  mpfr_init(ip_tmp2);
  mpfr_init(ip_tmp3);
  mpfr_init(ip_tmp4);
  mpfr_init(sinc_tmp);
  mpfr_init(mpfr_two_pi_B);
  mpfr_div_d(mpfr_two_pi_B,mpfr_pi,INTER_A,GMP_RNDN);
  mpfr_init(newton_df);
  mpfr_init(two_101);
  mpfr_set_ui(two_101,1,GMP_RNDN);
  mpfr_mul_2ui(two_101,two_101,OP_ACC,GMP_RNDN);
  mpfr_init(znew_guess);
  mpfr_init(zold_guess);
  mpfr_init(ztn);
  mpfr_init(tl);
  mpfr_init(tm);
  mpfr_init(tr);
  mpfr_init(fl);
  mpfr_init(fm);
  mpfr_init(fr);
  mpfr_init(sp_t);
  mpfr_init(sp_ft);
  mpfr_init(sp_t1);
  mpfr_init(sp_ft1);
  mpfr_init(last_t);
  mpfr_init(msin);
  mpfr_init(mcos);
  mpfr_init(fmid);
  init_out_bytes();
  // we assume converting mpz_t to long unsigned int
  // gives us 64 bits
  assert(sizeof(long unsigned int)==8);

  mpfi_init(A);
  mpfi_init(A1);
  mpfi_c_init(c_tmp);
  mpfi_init(msinc_tmp);
  mpfi_init(mpfi_two_pi_B);
  mpfi_div_d(mpfi_two_pi_B,mpfi_pi,INTER_A);
  mpfi_init(mip_tmp);
  mpfi_init(minter_tmp);
  mpfi_init(minter_tmp1);
  mpfi_init(zlow);
  mpfi_init(zhigh);
  mpfi_init(ztmp);
  mpfi_init(intererr);
  mpfi_init2(init_tmp,HI_PREC);
  mpfi_init2(init_tmp1,HI_PREC);
  mpfi_init2(sqrt_pi,HI_PREC);
  mpfi_init2(half_log_2_pi,HI_PREC);
  mpfi_const_pi(sqrt_pi);
  mpfi_init2(pi_hi,HI_PREC);
  mpfi_set(pi_hi,sqrt_pi);
  mpfi_set(half_log_2_pi,sqrt_pi);
  mpfi_mul_2ui(half_log_2_pi,half_log_2_pi,1);
  mpfi_log(half_log_2_pi,half_log_2_pi);
  mpfi_div_2ui(half_log_2_pi,half_log_2_pi,1);
  mpfi_sub_d(half_log_2_pi,half_log_2_pi,0.25); // 1/2log(2 Pi)-1/4
  mpfi_sqrt(sqrt_pi,sqrt_pi);
  mpfi_init2(g_tmp1,HI_PREC);
  mpfi_init2(g_abs,HI_PREC);
  mpfi_init2(g_re,HI_PREC);
  mpfi_init2(g_im,HI_PREC);
  mpfi_init2(m1,HI_PREC);
  mpfi_init2(m2,HI_PREC);
  mpfi_init2(g_z->re,HI_PREC);
  mpfi_init2(g_z->im,HI_PREC);
  mpfi_init2(g_z2->re,HI_PREC);
  mpfi_init2(g_z2->im,HI_PREC);
  init_lng();
  mpfi_init(misin);

  set_LNG_TERMS(t0);
  for(i=0,t=-N1/2*one_over_A1;i<N1;i++,t+=one_over_A1)
    {
      //if ((i&4095)==0)
      //printf("t=%e\n",t);
      mpfi_init(two_pi_t[i]);
      mpfi_mul_d(two_pi_t[i],mpfi_2_pi,-t);
      mpfi_c_init(g_vec[i]);
      mpfi_c_init(G_vec[i]);
      mpfi_c_init(f_vec[i]);
      //mpfi_c_init(df_vec[i]);
      mpfi_c_init(skn_vec[i]);
      mpfi_init2(exps[i],HI_PREC);
      
      //mpfi_c_lngamma_hi(g_vec[i],0.25,(t+t0)/2.0);    // lngamma(1/4+(t+t0)/2*I)
      mpfi_set_d(exps[i],t);
      mpfi_div_d(exps[i],exps[i],h);
      mpfi_sqr(exps[i],exps[i]);
      mpfi_div_2ui(exps[i],exps[i],1); //  t^2/(2h^2)
      //mpfi_sub(g_vec[i]->re,g_vec[i]->re,exps[i]);
      //mpfi_mul_d(G_vec[i]->re,mpfi_pi,(t+t0)/4.0); // Pi(t+t0)/4
      //mpfi_add(g_vec[i]->re,g_vec[i]->re,G_vec[i]->re);
      //mpfi_c_exp(g_vec[i],g_vec[i]);
      gamma_factor(g_vec[i],exps[i],t+t0);
      /*
	if((i%5000)==0)
	{
	  printf("t=%13.11e\n",t);
	  mpfi_c_printn(g_vec[i],80);
	  printf("rel errors = %d %d\n",mpfi_rel_error(g_vec[i]->re),mpfi_rel_error(g_vec[i]->im));
	}
      */
    }
  for(i=N1;i<N;i++)
    {
      mpfi_c_init(f_vec[i]);
      //mpfi_c_init(df_vec[i]);
    }
  for(i=0;i<M;i++)
    {
      mpfi_init(s_vec[i]);
      mpfi_init(sk_vec[i]);
      mpfi_init2(logisqrtpi[i],HI_PREC);
      mpfi_c_init(n_vec[i]);
      mpfi_mul_ui(logisqrtpi[i],sqrt_pi,i+1); // n*sqrt(Pi)
      mpfi_log(logisqrtpi[i],logisqrtpi[i]);        // log (n*sqrt(Pi))
      mpfi_set(s_vec[i],logisqrtpi[i]);
      mpfi_mul_d(init_tmp1,logisqrtpi[i],-t0); // -t0*log(n*sqrt(Pi))
      mpfi_sin_cos_hi(n_vec[i]->im,n_vec[i]->re,init_tmp1);
      mpfi_set_ui(A,1);
      mpfi_div_ui(A,A,i+1);
      mpfi_sqrt(A,A);
      mpfi_init(sqrts[i]);
      mpfi_set(sqrts[i],A);
      mpfi_c_mul_i(n_vec[i],n_vec[i],A);    //  / sqrt(n)
      mpfi_div(s_vec[i],s_vec[i],mpfi_2_pi); // log(n*sqrt(Pi))/(2*Pi)
    }

  // set up the omega vectors for FFTs
  initfft(N/2,ws_r); // ws_r[i]=e(i/N)
  for(i=0;i<N1/2;i++)
    {
      mpfi_c_init(ws1_f[i]);
      mpfi_c_init(ws1_r[i]);
      mpfi_c_set(ws1_r[i],ws_r[i*UPSAM/2]); // ws1_r[i]=e(i/N1)
      mpfi_c_conj(ws1_f[i],ws1_r[i]);     // ws1_f[i]=e(-i/N1)
    }

  set_err(gtwiderr,gtwiderr_d);
  set_err(Gtwiderr,Gtwiderr_d);
  set_err(tayerr,tayerr_d);
  set_err(Fmaxerr,Fmaxerr_d);
  mpfi_set_d(intererr,intererr_d);
  mpfi_neg(A,intererr);
  mpfi_put(intererr,A);
  mpfi_set_ui(A,M);
  mpfi_sqrt(A,A);
  mpfi_mul_ui(A,A,2);
  mpfi_sub_ui(A,A,1); // 2sqrt(M)-1
  mpfi_c_mul_i(tayerr,tayerr,A);
  //printf("Total Taylor Error set to ");mpfi_printn(tayerr->re,30);
  set_err(fhatsumerr,fhatsumerr_d);
  set_err(fhattwiderr,fhattwiderr_d);
  set_err(ftwiderr,ftwiderr_d);
  mpfi_set_ui(A,1);
  mpfi_div_d(A,A,one_over_A); // don't use A as a spare variable anymore
  mpfi_set_ui(A1,1);
  mpfi_div_d(A1,A1,one_over_A1);
  calc_buck();
  mpfi_init(atan1_tmp1);
  mpfi_init(atan1_tmp2);
  mpfi_init(im_1);
  mpfi_init(im_2);
  mpfi_init(im_t);
  mpfi_init(im_err);
  mpfi_init(tm_h);
  mpfi_init(tm_res);
  mpfi_init(tm_tmp);
  mpfi_init(tm_tmp1);
  mpfi_init(im_t0);
  mpfi_init(im_t1);
  mpfi_init(tm_t0);
  mpfi_init(tm_t1);
  mpfr_init(tm_mpfr);
}

// run with a new value of t0
// just need to set g_vec and n_vec up
void re_init(double t0)
{
  int i;
  double t;
  //mpfi_init(A);
  set_LNG_TERMS(t0);
  for(i=0,t=-N1/2*one_over_A1;i<N1;i++,t+=one_over_A1)
    {
      //mpfi_set_d(g_vec[i]->re,-t-t0);mpfi_set_ui(g_vec[i]->im,i);
      //mpfi_c_lngamma_hi(g_vec[i],0.25,(t+t0)/2.0);    // lngamma(1/4+(t+t0)/2*I)
      //mpfi_c_exp(g_vec[i],g_vec[i]);
      //mpfi_mul_d(c_tmp->re,mpfi_pi,(t+t0)/4.0);
      //mpfi_exp(c_tmp->re,c_tmp->re);
      //mpfi_c_mul_i(g_vec[i],g_vec[i],c_tmp->re);
      //mpfi_add(g_vec[i]->re,g_vec[i]->re,c_tmp->re);
      //mpfi_c_div_i(g_vec[i],g_vec[i],exps[i]);
      //mpfi_sub(g_vec[i]->re,g_vec[i]->re,exps[i]);
      //mpfi_c_exp(g_vec[i],g_vec[i]);
      gamma_factor(g_vec[i],exps[i],t+t0);
    }
  for(i=0;i<M;i++)
    {
      //mpfi_mul_ui(n_vec[i]->re,mpfi_sqrt_pi,(i+1)); // ((i+1)*sqrt(Pi))
      //mpfi_log(n_vec[i]->re,n_vec[i]->re);
      mpfi_mul_d(init_tmp1,logisqrtpi[i],-t0);
      mpfi_sin_cos_hi(n_vec[i]->im,n_vec[i]->re,init_tmp1);
      //mpfi_set_ui(c_tmp->re,i+1);
      //mpfi_sqrt(c_tmp->re,c_tmp->re);
      mpfi_c_mul_i(n_vec[i],n_vec[i],sqrts[i]);    //  / sqrt(i+1)
    }
}

// on entry g_vec[i]=g(t)*(-2*Pi*t*I)^k/k!=g(t;k)/k!
// on exit G_vec[i]=G^(k)(i/B)/k!
void G_k(int k)
{
  int i;
  double t,max_t;

  for(i=0;i<N1;i++)
    mpfi_c_add(G_vec[i],g_vec[i],gtwiderr); // G=g~

  //simple_idft(0);
  //simple_idft(965);
  //simple_idft(N-1);

  //nifft(G_vec,N,ws_r);
  fft(G_vec,N1,ws1_f);

  for(i=0;i<=N1/2;i++)
    {
      mpfi_c_div_i(G_vec[i],G_vec[i],A1);
      mpfi_c_add(G_vec[i],G_vec[i],Gtwiderr);
      if(i&1)
	mpfi_c_neg(G_vec[i],G_vec[i]);
    }
  for(i=N1/2+1;i<N1;i++)
    mpfi_c_set_ui(G_vec[i],0,0);

  if(k<K-1)
    for(i=0;i<N1;i++)
      {
	mpfi_c_muli(g_vec[i]); // g(t;k)*i/k!
	mpfi_c_mul_i(g_vec[i],g_vec[i],two_pi_t[i]); // g(t;k+1)/k!
	mpfi_c_div_ui(g_vec[i],g_vec[i],k+1);        // g(t;k+1)/(k+1)!
	/*if((i%10000)==0)
	  {
	    printf("error in g_vec in G_k(%d)= %d %d.\n",k,mpfi_rel_error(g_vec[i]->re),mpfi_rel_error(g_vec[i]->im));
	    }*/
      }

}


// on entry sk_vec = (log(nsqrt(Pi))/2Pi-u_m)^k (or undefined if k==0)
// n_vec[i] = (i+1)^(-1/2)*(sqrt(Pi)*(i+1))^-it0
// on exit skn_vec = sum nsqrt(pi)^(-ito)/sqrt(n)*(log(nsqrt(Pi))/2Pi-u_m)^k
//         sk_vec = (log..-u_m)^(k+1)
void make_skn(int k)
{
  int i,pos;
  for(i=0;i<N1;i++)
    mpfi_c_zero(skn_vec[i]);
  //debug;      
  switch(k)
    {
    case 0: // set sk_vec = s_vec, skn_vec=n_vec
      for(i=0;i<M;i++)
	{
	  //printf("buck[i]=%d\n",buck[i]);
	  mpfi_c_add(skn_vec[buck[i]],skn_vec[buck[i]],n_vec[i]);
	  mpfi_set(sk_vec[i],s_vec[i]);
	}
      break;
    case K-1: // no need to increment sk_vec
      for(i=0;i<M;i++)
	{
	  mpfi_c_mul_i(c_tmp,n_vec[i],sk_vec[i]);
	  mpfi_c_add(skn_vec[buck[i]],skn_vec[buck[i]],c_tmp);
	}
      break;
    default:
      for(i=0;i<M;i++)
	{
	  mpfi_c_mul_i(c_tmp,n_vec[i],sk_vec[i]);
	  mpfi_c_add(skn_vec[buck[i]],skn_vec[buck[i]],c_tmp);
	  mpfi_mul(sk_vec[i],sk_vec[i],s_vec[i]);
	}
    }
}

void my_convolve (mpfi_c_t *res, mpfi_c_t *v1, mpfi_c_t *v2, int n, mpfi_c_t *ws_r, mpfi_c_t *ws_f)
{
  int i;
  fft(v1,n,ws_r);
  fft(v2,n,ws_r);
  for(i=0;i<n;i++)
    mpfi_c_mul(res[i],v1[i],v2[i]);
  fft(res,n,ws_f);
  for(i=0;i<n;i++)
    mpfi_c_div_ui(res[i],res[i],n);
}

// f_vec+=G(k)(x+u_m)/k!*S^(k)_m
void do_conv (int k)
{
  int i;
  make_skn(k);    
  //debug;
  if(k==0)
    my_convolve(f_vec,skn_vec,G_vec,N1,ws1_r,ws1_f);
  else
    {
      my_convolve(G_vec,skn_vec,G_vec,N1,ws1_r,ws1_f);
      for(i=0;i<=N1/2;i++)
	mpfi_c_add(f_vec[i],f_vec[i],G_vec[i]);
    }
}

// compute exp(-t^2/2H^2)
void inter_gaussian(mpfr_ptr res, mpfr_ptr t)
{
  //mpfr_printf("inter_gaussian called with t=%.40Re\n",t);
  mpfr_div_d(res,t,H,GMP_RNDN);
  mpfr_sqr(res,res,GMP_RNDN);
  mpfr_div_d(res,res,-2.0,GMP_RNDN);
  mpfr_exp(res,res,GMP_RNDN);
  //mpfr_printf("inter_gaussian returning %.40Re\n",res);
}

int sincp;
// compute sinc(2*B*Pi*t) into ress, cos(2*B*Pi*t) into resc
void inter_sinc_cos(mpfr_ptr ress, mpfr_ptr resc, mpfr_ptr t)
{
  mpfr_mul(sinc_tmp,t,mpfr_two_pi_B,GMP_RNDN);
  if(sincp)
    {
      mpfr_sin_cos(msin,mcos,sinc_tmp,GMP_RNDN);
      sincp=FALSE;
    }
  else
    {
      mpfr_neg(msin,msin,GMP_RNDN);
      mpfr_neg(mcos,mcos,GMP_RNDN);
    }
  mpfr_set(resc,mcos,GMP_RNDN);
  mpfr_div(ress,msin,sinc_tmp,GMP_RNDN);
}

// compute exp(-t^2/2H^2)
void mpfi_inter_gaussian(mpfi_ptr res, mpfi_ptr t)
{
  mpfi_div_d(res,t,H);
  mpfi_sqr(res,res);
  mpfi_div_d(res,res,-2.0);
  mpfi_exp(res,res);
}

// compute sinc(2*B*Pi*t) into ress
void mpfi_inter_sinc(mpfi_ptr ress, mpfi_ptr t)
{
  mpfi_mul(msinc_tmp,t,mpfi_two_pi_B);
  if(sincp)
    {
      mpfi_sin(misin,msinc_tmp);
      sincp=FALSE;
    }
  else
    mpfi_neg(misin,misin);
  mpfi_div(ress,misin,msinc_tmp);
}

void mpfi_inter_point(mpfi_ptr f_res, mpfi_c_t *f_vec, mpfi_ptr t, int f_ptr)
{
  mpfi_inter_gaussian(f_res,t);
  mpfi_inter_sinc(mip_tmp,t);
  mpfi_mul(f_res,f_res,mip_tmp);
  mpfi_mul(f_res,f_res,f_vec[f_ptr]->re);
}

// t is distance from t0 to t implied by f_ptr
// if fd_res is NULL, don't calc differential
void inter_point(mpfr_ptr f_res, mpfr_ptr fd_res, mpfi_c_t *f_vec, mpfr_ptr t, int f_ptr)
{
  mpfi_get_fr(ip_tmp1,f_vec[f_ptr]->re);
  inter_gaussian(ip_tmp2,t);
  inter_sinc_cos(ip_tmp3,ip_tmp4,t);
  // compute f_res
  mpfr_mul(ip_tmp2,ip_tmp1,ip_tmp2,GMP_RNDN); //f*exp
  mpfr_mul(f_res,ip_tmp2,ip_tmp3,GMP_RNDN); //f*exp*sinc
  // compute fd_res
  if(!fd_res)
    return;
  mpfr_div(fd_res,ip_tmp4,t,GMP_RNDN); // cos/t
  mpfr_set_ui(ip_tmp4,1,GMP_RNDN);
  mpfr_div(ip_tmp4,ip_tmp4,t,GMP_RNDN); // 1/t
  mpfr_div_d(ip_tmp1,t,H*H,GMP_RNDN); // t/H^2
  mpfr_add(ip_tmp4,ip_tmp4,ip_tmp1,GMP_RNDN); // 1/t+t/H^2
  mpfr_mul(ip_tmp4,ip_tmp4,ip_tmp3,GMP_RNDN); // sinc*(1/t+t/H^2)
  mpfr_sub(fd_res,fd_res,ip_tmp4,GMP_RNDN);   // cos/t-sinc*(1/t+t/H^2)
  mpfr_mul(fd_res,fd_res,ip_tmp2,GMP_RNDN);   // f*exp*(cos/t-sinc....)
  //
  // WARNING. COMPLETELY ARBITARY NEGATION TO MAKE f'(t) CORRECT SIGN
  //
  mpfr_neg(fd_res,fd_res,GMP_RNDN);
}

void mpfi_inter_t(mpfi_ptr f_res, mpfi_c_t *f_vec, mpfi_ptr t_ptr)
{
  int i,j,nearest_t0;
  //printf("In inter_t\n");
  double t0_d=mpfi_get_d(t_ptr),dt0;
  mpfi_get_left(inter_tmp,t_ptr);
  mpfr_floor(inter_tmp,inter_tmp);
  mpfi_set_fr(minter_tmp,inter_tmp);
  if(mpfi_cmp(minter_tmp,t_ptr)==0) // t_ptr exactly coincides with a lattice pt
    nearest_t0=t0_d+1.5;
  else
    nearest_t0=t0_d+0.5;
  //printf("Nearest t0=%ld\n",nearest_t0);
  mpfi_set_ui(f_res,0);
  sincp=TRUE;
  if(nearest_t0>t0_d) // start from right of t0
    {
      for(i=nearest_t0,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  mpfi_sub_ui(minter_tmp,t_ptr,i);
	  mpfi_mul_d(minter_tmp,minter_tmp,-one_over_A);
	  mpfi_inter_point(minter_tmp1,f_vec,minter_tmp,i);
	  mpfi_add(f_res,f_res,minter_tmp1);
	}
      sincp=TRUE;
      for(i=nearest_t0-INTER_SPACING,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  mpfi_sub_ui(minter_tmp,t_ptr,i);
	  mpfi_mul_d(minter_tmp,minter_tmp,-one_over_A);
	  mpfi_inter_point(minter_tmp1,f_vec,minter_tmp,i);
	  mpfi_add(f_res,f_res,minter_tmp1);
	}
    }
  else // start from left
    {
      for(i=nearest_t0+INTER_SPACING,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  mpfi_sub_ui(minter_tmp,t_ptr,i);
	  mpfi_mul_d(minter_tmp,minter_tmp,-one_over_A);
	  mpfi_inter_point(minter_tmp1,f_vec,minter_tmp,i);
	  mpfi_add(f_res,f_res,minter_tmp1);
	}
      sincp=TRUE;
      for(i=nearest_t0,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  mpfi_sub_ui(minter_tmp,t_ptr,i);
	  mpfi_mul_d(minter_tmp,minter_tmp,-one_over_A);
	  mpfi_inter_point(minter_tmp1,f_vec,minter_tmp,i);
	  mpfi_add(f_res,f_res,minter_tmp1);
	}
    }
  mpfi_add(f_res,f_res,intererr);
}

// t0 points into f_vec
// let t=t0+(t_ptr-N/2)*one_over_A
// returns f(t) in f_res, f'(t) in fd_res
// return f_vec[t0]->re if t_ptr is integral 
void inter_t(mpfr_ptr f_res, mpfr_ptr fd_res, mpfi_c_t *f_vec, mpfr_ptr t_ptr)
{
  int i,j,nearest_t0;
  //printf("In inter_t\n");
  double t0_d=mpfr_get_d(t_ptr,GMP_RNDN),dt0;
  mpfr_floor(inter_tmp,t_ptr);
  if(mpfr_cmp(inter_tmp,t_ptr)==0) // t_ptr exactly coincides with a lattice pt
    nearest_t0=t0_d+1.5;
  else
    nearest_t0=t0_d+0.5;
  //printf("Nearest t0=%ld\n",nearest_t0);
  mpfr_set_ui(f_res,0,GMP_RNDN);
  if (fd_res) mpfr_set_ui(fd_res,0,GMP_RNDN);
  sincp=TRUE;
  if(nearest_t0>t0_d) // start from right of t0
    {
      for(i=nearest_t0,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  mpfr_sub_ui(inter_tmp,t_ptr,i,GMP_RNDN);
	  mpfr_mul_d(inter_tmp,inter_tmp,-one_over_A,GMP_RNDN);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i);
	  mpfr_add(f_res,f_res,inter_tmp1,GMP_RNDN);
	  if (fd_res) mpfr_add(fd_res,fd_res,inter_tmp2,GMP_RNDN);
	}
      sincp=TRUE;
      for(i=nearest_t0-INTER_SPACING,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  mpfr_sub_ui(inter_tmp,t_ptr,i,GMP_RNDN);
	  mpfr_mul_d(inter_tmp,inter_tmp,-one_over_A,GMP_RNDN);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i);
	  mpfr_add(f_res,f_res,inter_tmp1,GMP_RNDN);
	  if (fd_res) mpfr_add(fd_res,fd_res,inter_tmp2,GMP_RNDN);
	}
    }
  else // start from left
    {
      for(i=nearest_t0+INTER_SPACING,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  //printf("Doing right tail with i=%d\n",i);
	  //
	  mpfr_sub_ui(inter_tmp,t_ptr,i,GMP_RNDN);
	  mpfr_mul_d(inter_tmp,inter_tmp,-one_over_A,GMP_RNDN);
	  //printf("Calling ip\n");
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i);
	  //mpfr_printf("ip returned %.40Re\n",inter_tmp1);
	  mpfr_add(f_res,f_res,inter_tmp1,GMP_RNDN);
	  //
	  if (fd_res) mpfr_add(fd_res,fd_res,inter_tmp2,GMP_RNDN);
	  //
	}
      sincp=TRUE;
      for(i=nearest_t0,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  mpfr_sub_ui(inter_tmp,t_ptr,i,GMP_RNDN);
	  mpfr_mul_d(inter_tmp,inter_tmp,-one_over_A,GMP_RNDN);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i);
	  mpfr_add(f_res,f_res,inter_tmp1,GMP_RNDN);
	  if (fd_res) mpfr_add(fd_res,fd_res,inter_tmp2,GMP_RNDN);
	}
    }
}

#define is_exact(f) (f==MPFI_FLAGS_BOTH_ENDPOINTS_EXACT)



void scale_t(mpfr_ptr res, mpfr_ptr t)
{
  mpfr_mul_d(res,t,one_over_A,GMP_RNDN);
  mpfr_mul_2ui(res,res,OP_ACC,GMP_RNDN);
  mpfr_add_d(res,res,0.5,GMP_RNDN);
  mpfr_floor(res,res);
}

void bracket_t(mpfi_ptr res1, mpfi_ptr res2, mpfr_ptr tn, mpfr_ptr t_now)
{
  scale_t(tn,t_now);
  assert(is_exact(mpfi_set_fr(res1,tn)));
  mpfi_add_d(res2,res1,0.5);
  mpfi_sub_d(res1,res1,0.5);
  mpfi_div_2ui(res1,res1,OP_ACC);
  mpfi_div_2ui(res2,res2,OP_ACC);
  mpfi_mul(res1,res1,A);
  mpfi_mul(res2,res2,A);
}



// NB destroys old_guess
// on entry left and right bracket the zero
// on exit new_guess = improved estimate
void newton(mpfr_ptr new_guess, mpfr_ptr old_guess, mpfi_c_t *f_vec, int num_its)
{
  int i;
  //printf("In newton.\n");
  for(i=0;i<num_its;i++)
    {
      //mpfr_printf("old_guess=%.40Re\n",old_guess);
      inter_t(new_guess,newton_df,f_vec,old_guess); // new=f(t) df=f'(t)
      //mpfr_printf("f()      =%.40Re\n",new_guess);
      mpfr_div(newton_df,new_guess,newton_df,GMP_RNDN);  // f/f'
      mpfr_div_d(newton_df,newton_df,one_over_A,GMP_RNDN); // f/f' normalised
      mpfr_sub(old_guess,old_guess,newton_df,GMP_RNDN); // old-f/f'
    }
  mpfr_set(new_guess,old_guess,GMP_RNDN);
}

// start at point between tl and tr (where f has value fl and fr)
void newton2 (mpfr_ptr new_guess, mpfr_ptr old_guess, mpfi_c_t *f_vec, mpfr_ptr tl, mpfr_ptr tr, mpfr_ptr fl, mpfr_ptr fr, int num_its)
{
  int i;
  mpfr_add(old_guess,tl,tr,GMP_RNDN);
  mpfr_div_2ui(old_guess,old_guess,1,GMP_RNDN);
  newton(new_guess,old_guess,f_vec,num_its);
}

int num_stat_pts(mpfi_c_t *f_vec, int start, int end)
{
  int i=start,res=0,c1,c2;
  c1=mpfi_cmp(f_vec[0]->re,f_vec[1]->re);
  while(i<end-2)
    {
      c2=mpfi_cmp(f_vec[i+1]->re,f_vec[i+2]->re);

      if((c1>0)&&(c2<0)&&(sign(f_vec[i+1]->re)==POS))
	res+=2;
      else
	{
	  if((c1<0)&&(c2>0)&&(sign(f_vec[i+1]->re)==NEG))
	    res+=2;
	}
      c1=c2;
      i++;
    }
  return(res);
}

void binary_chop(mpfr_ptr tmid, mpfr_ptr tl, mpfr_ptr tr, mpfr_ptr fl, mpfr_ptr fr, mpfi_c_t *f_vec)
{
  int i,signl=mpfr_sgn(fl),signr=mpfr_sgn(fr),signm;
  //assert(signl!=signr);
  mpfr_printf("in binary chop with tl=%.40Re and tr=%.40Re\n",tl,tr); 
  for(i=0;i<OP_ACC+10;i++) // 10 is random but should be ample
    {
      mpfr_add(tmid,tl,tr,GMP_RNDN);
      mpfr_div_2ui(tmid,tmid,1,GMP_RNDN);
      inter_t(fmid,NULL,f_vec,tmid);
      signm=mpfr_sgn(fmid);
      if(signm!=signl)
	{
	  mpfr_set(tr,tmid,GMP_RNDN);
	  //signr=signm;
	  continue;
	}
      if(signm!=signr)
	{
	  mpfr_set(tl,tmid,GMP_RNDN);
	  //signl=signm;
	  continue;
	}
      printf("Bad binary chop at i=%d\n",i);
    }
  mpfr_add(tmid,tl,tr,GMP_RNDN);
  mpfr_div_2ui(tmid,tmid,1,GMP_RNDN);
  mpfr_printf("Exiting binary chop with res=%.40Re\n",tmid);
}


// check rigorously for zeros without using stat points
int zeros(mpfi_c_t *f_vec, int start, int end, FILE *outfile)
{
  long int i=start,count=0,last_pos;
  sign_t last_sign=sign(f_vec[i]->re),this_sign,res1_sign,res2_sign;
  mpfr_set_ui(last_t,start,GMP_RNDN);
  scale_t(last_t,last_t);
  while(last_sign==UNK) // skip initial section of UNKs, will break if all UNKs!
    last_sign=sign(f_vec[++i]->re);
  last_pos=i;
  for(i++;i<=end;i++)
    {
      this_sign=sign(f_vec[i]->re);
      if(this_sign==last_sign)
	{
	  last_pos=i;
	  continue;
	}
      if(this_sign==UNK)
	continue;
      // A sign change
      last_sign=this_sign;
      mpfi_get_fr(fl,f_vec[last_pos]->re);
      mpfi_get_fr(fr,f_vec[i]->re);
      mpfr_set_ui(tl,last_pos,GMP_RNDN);
      mpfr_set_ui(tr,i,GMP_RNDN);
      newton2(znew_guess,zold_guess,f_vec,tl,tr,fl,fr,NEWTON_ITS);
      if((mpfr_cmp(znew_guess,tl)<0)||(mpfr_cmp(znew_guess,tr)>0)) // newton wandered off soemwhere
	{
	  printf("Newton failed to converge, binary chopping.\n");
	  binary_chop(znew_guess,tl,tr,fl,fr,f_vec);
	}
      bracket_t(zlow,zhigh,ztn,znew_guess);
      mpfi_inter_t(ztmp,f_vec,zlow);
      res1_sign=sign(ztmp);
      mpfi_inter_t(ztmp,f_vec,zhigh);
      res2_sign=sign(ztmp);
      if((res1_sign&res2_sign)!=0) // not a valid sign change
	{
	  printf("Not a valid sign change (non-stationary).\n");
	  printf("right hand value =");mpfi_printn(ztmp,40);
	  mpfr_printf("At location %.40Re\n",zhigh);
	  mpfi_inter_t(ztmp,f_vec,zlow);
	  printf("left hand value  =");mpfi_printn(ztmp,40);
	  mpfr_printf("At location %.40Re\n",zlow);
	  printf("sign change was between %ld and %ld\n",last_pos,i);
	  return(count);
	}
      mpfr_sub(last_t,ztn,last_t,GMP_RNDN);
      if(mpfr_sgn(last_t)<=0)
	{
	  printf("Co-incident or switched zeros.\n");
	  mpfr_printf("last_t=%.40Re\nztn=%.40Re\ntl=%.40Re\ntr=%.40Re\nfl=%.40Re\nfr=%.40Re\n",last_t,ztn,tl,tr,fl,fr);
	  return(count);
	}
      count++;
      out_bytes(last_t,outfile);
      mpfr_set(last_t,ztn,GMP_RNDN);
      last_pos=i;
    }
  return(count);
}

int resolve_stat_point(mpfr_ptr tl, mpfr_ptr ftl, mpfr_ptr tm, mpfr_ptr ftm, mpfr_ptr tr, mpfr_ptr ftr, int left, sign_t this_sign, mpfi_c_t *f_vec)
{
  dir_t dir1,dir2,dir3,dir4;
  printf("In resolve_stat_points at %d %d %d.\n",left,left+1,left+2);
  //mpfr_t t,ft;
  //mpfr_init(t);mpfr_init(ft);
  mpfr_set_ui(tl,left,GMP_RNDN);
  mpfr_set_ui(tm,left+1,GMP_RNDN);
  mpfr_set_ui(tr,left+2,GMP_RNDN);
  mpfi_get_fr(ftl,f_vec[left]->re);
  mpfi_get_fr(ftm,f_vec[left+1]->re);
  mpfi_get_fr(ftr,f_vec[left+2]->re);

  while(TRUE)
    {
      //mpfr_printf("Looping with t:\nleft=%.40Re\nmid=%.40Re\nright=%.40Re\n",tl,tm,tr);
      //mpfr_printf("F(t):\nleft=%.40Re\nmid=%.40Re\nright=%.40Re\n",ftl,ftm,ftr);

      mpfr_add(sp_t,tl,tm,GMP_RNDN);
      mpfr_div_2ui(sp_t,sp_t,1,GMP_RNDN);
      inter_t(sp_ft,NULL,f_vec,sp_t);
      if((mpfr_sign(sp_ft)&this_sign)==0)
	{
	  //mpfr_printf("Sign change found at t=%.40Re\n",sp_t);
	  mpfr_set(tr,tm,GMP_RNDN);
	  mpfr_set(ftr,ftm,GMP_RNDN);
	  mpfr_set(tm,sp_t,GMP_RNDN);
	  mpfr_set(ftm,sp_ft,GMP_RNDN);
	  return(TRUE);
	}
      dir1=mpfr_dir(ftl,sp_ft);
      dir2=mpfr_dir(sp_ft,ftm);
      if(((this_sign==POS)&&(dir1==DOWN)&&(dir2==UP))||((this_sign==NEG)&&(dir1==UP)&&(dir2==DOWN)))
	{
	  mpfr_set(tr,tm,GMP_RNDN);
	  mpfr_set(ftr,ftm,GMP_RNDN);
	  mpfr_set(tm,sp_t,GMP_RNDN);
	  mpfr_set(ftm,sp_ft,GMP_RNDN);
	  continue;
	}

      mpfr_add(sp_t1,tm,tr,GMP_RNDN);
      mpfr_div_2ui(sp_t1,sp_t1,1,GMP_RNDN);
      inter_t(sp_ft1,NULL,f_vec,sp_t1);
      if((mpfr_sign(sp_ft1)&this_sign)==0)
	{
	  //mpfr_printf("Sign change found at t=%.40Re\n",sp_t);
	  mpfr_set(tl,tm,GMP_RNDN);
	  mpfr_set(ftl,ftm,GMP_RNDN);
	  mpfr_set(tm,sp_t1,GMP_RNDN);
	  mpfr_set(ftm,sp_ft1,GMP_RNDN);
	  return(TRUE);
	}
      dir3=mpfr_dir(ftm,sp_ft1);
      dir4=mpfr_dir(sp_ft1,ftr);
      if(((this_sign==POS)&&(dir2==DOWN)&&(dir3==UP))||((this_sign==NEG)&&(dir2==UP)&&(dir3==DOWN)))
	{
	  mpfr_set(tl,sp_t,GMP_RNDN);
	  mpfr_set(ftl,sp_ft,GMP_RNDN);
	  mpfr_set(tr,sp_t1,GMP_RNDN);
	  mpfr_set(ftr,sp_ft1,GMP_RNDN);
	  continue;
	}
      if(((this_sign==POS)&&(dir3==DOWN)&&(dir4==UP))||((this_sign==NEG)&&(dir3==UP)&&(dir4==DOWN)))
	{
	  mpfr_set(tl,tm,GMP_RNDN);
	  mpfr_set(ftl,ftm,GMP_RNDN);
	  mpfr_set(tm,sp_t1,GMP_RNDN);
	  mpfr_set(ftm,sp_ft1,GMP_RNDN);
	  continue;
	}
      printf("Stat point Process failed to converge.\n");
      return(FALSE);
    }
}


int zeros_st(mpfi_c_t *f_vec, int start, int end, FILE *outfile)
{
  long int i=start+1,count=0,last_pos=start;
  dir_t last_dir,this_dir;
  sign_t last_sign,this_sign,res1_sign,res2_sign;
  mpfr_set_ui(last_t,start,GMP_RNDN);
  scale_t(last_t,last_t);  
  this_sign=sign(f_vec[start]->re);
  this_dir=UNK;
  for(;i<=end;i++)
    {
      
      last_sign=this_sign;
      last_dir=this_dir;
      this_sign=sign(f_vec[i]->re);
      this_dir=dir(f_vec[i-1]->re,f_vec[i]->re);
      if((this_sign&last_sign)==0) // valid sign change
	{
	  /*
	  mpfr_set_ui(tl,last_pos,GMP_RNDN);
	  mpfr_set_ui(tr,i,GMP_RNDN);
	  mpfi_get_fr(fl,f_vec[last_pos]->re);
	  mpfi_get_fr(fr,f_vec[i]->re);
	  newton2(znew_guess,zold_guess,f_vec,tl,tr,fl,fr,NEWTON_ITS);
	  if((mpfr_cmp(znew_guess,tl)<0)||(mpfr_cmp(znew_guess,tr)>0)) // newton wandered off soemwhere
	    {
	      printf("Newton failed to converge, binary chopping.\n");
	      binary_chop(znew_guess,tl,tr,fl,fr,f_vec);
	    }
	  //printf("Result of N-R is ");
	  //mpfr_out_str(NULL,10,0,znew_guess,GMP_RNDN);
	  //printf("\n");
	  //printf("Count %d-Sign change at %d-%d.\nBest guess= ",count,last_pos,i);
	  //mpfr_out_str(NULL,10,0,znew_guess,GMP_RNDN);
	  //printf("\n");
	  bracket_t(zlow,zhigh,ztn,znew_guess);
	  mpfi_inter_t(ztmp,f_vec,zlow);
	  res1_sign=sign(ztmp);
	  mpfi_inter_t(ztmp,f_vec,zhigh);
	  res2_sign=sign(ztmp);
	  if((res1_sign&res2_sign)!=0) // not a valid sign change
	    {
	      binary_chop(znew_guess,tl,tr,fl,fr,f_vec);
	      bracket_t(zlow,zhigh,ztn,znew_guess);
	      mpfi_inter_t(ztmp,f_vec,zlow);
	      printf("Left hand f(t)=");mpfi_printn(ztmp,50);
	      res1_sign=sign(ztmp);
	      mpfi_inter_t(ztmp,f_vec,zhigh);
	      mpfr_printf("Right hand f(t)=");mpfi_printn(ztmp,50);
	      res2_sign=sign(ztmp);
	      if((res1_sign&res2_sign)!=0) // still not a valid sign change
		{
		  printf("Failed to resolve sign (non-stat) even after binary chop.\n");
		 
		  return(count);
		}
	    }
	  mpfr_sub(last_t,ztn,last_t,GMP_RNDN);
	  if(mpfr_sgn(last_t)<=0)
	    {
	      printf("Co-incident or switched zeros.\n");
	      mpfr_printf("last_t=%.40Re\nztn=%.40Re\ntl=%.40Re\ntr=%.40Re\nfl=%.40Re\nfr=%.40Re\n",last_t,ztn,tl,tr,fl,fr);
	      return(count);
	    }
	  out_bytes(last_t,outfile);
	  mpfr_set(last_t,ztn,GMP_RNDN);
	  */	  
	  last_pos=i;
	  count++;
	  continue;
	}
      if(this_sign!=UNK)
	last_pos=i;
      if((this_dir&last_dir)==0) // maximum or minimum
	{
	  if(((this_dir==UP)&&(this_sign==POS))||((this_dir==DOWN)&&(this_sign==NEG)))
	    {
	      printf("Stat point found at %ld.\n",i);
	      //printf("f[%d]=",i-2);mpfi_printn(f_vec[i-2]->re,30);
	      //printf("f[%d]=",i-1);mpfi_printn(f_vec[i-1]->re,30);
	      //printf("f[%d]=",i);mpfi_printn(f_vec[i]->re,30);
	      if(!resolve_stat_point(tl,fl,tm,fm,tr,fr,i-2,this_sign,f_vec))
		return(count);
	      /*
	      mpfr_printf("zeros found between %.40Re %.40Re %.40Re\n",tl,tm,tr);
	      mpfr_printf("with values %.40Re %.40Re %.40Re\n",fl,fm,fr);
	      // do the left one
	      newton2(znew_guess,zold_guess,f_vec,tl,tm,fl,fm,NEWTON_ITS);
	      if((mpfr_cmp(znew_guess,tl)<0)||(mpfr_cmp(znew_guess,tr)>0)) // newton wandered off soemwhere
		{
		  printf("Newton failed to converge, binary chopping.\n");
		  binary_chop(znew_guess,tl,tr,fl,fr,f_vec);
		}

	      //mpfr_printf("Newton (left) returned %.40Re\n",znew_guess);
	      bracket_t(zlow,zhigh,ztn,znew_guess);
	      mpfi_inter_t(ztmp,f_vec,zlow);
	      res1_sign=sign(ztmp);
	      mpfi_inter_t(ztmp,f_vec,zhigh);
	      res2_sign=sign(ztmp);
	      if((res1_sign&res2_sign)!=0)
		{
		  binary_chop(znew_guess,tl,tm,fl,fm,f_vec);
		  bracket_t(zlow,zhigh,ztn,znew_guess);
		  mpfi_inter_t(ztmp,f_vec,zlow);
		  res1_sign=sign(ztmp);
		  mpfi_inter_t(ztmp,f_vec,zhigh);
		  res2_sign=sign(ztmp);
		  if((res1_sign&res2_sign)!=0) // still not a valid sign change
		    {
		      printf("Failed to resolve sign (stat) even after binary chop.\n");
		      return(count);
		    }
		}
	      mpfr_sub(last_t,ztn,last_t,GMP_RNDN);
	      if(mpfr_sgn(last_t)<=0)
		{
		  printf("Co-incident or switched zeros.\n");
		  mpfr_printf("last_t=%.40Re\nztn=%.40Re\ntl=%.40Re\ntr=%.40Re\nfl=%.40Re\nfr=%.40Re\n",last_t,ztn,tl,tr,fl,fr);
		  return(count);
		}
	      */
	      count++;
	      //out_bytes(last_t,outfile);
	      //mpfr_out_str(NULL,10,0,last_t,GMP_RNDN);
	      //mpfr_set(last_t,ztn,GMP_RNDN);
	      //mpfr_printf("zero found near %.50Re\n",ztn);
	      /*
	      newton2(znew_guess,zold_guess,f_vec,tm,tr,fm,fr,NEWTON_ITS);
	      if((mpfr_cmp(znew_guess,tl)<0)||(mpfr_cmp(znew_guess,tr)>0)) // newton wandered off soemwhere
		{
		  printf("Newton failed to converge, binary chopping.\n");
		  binary_chop(znew_guess,tl,tr,fl,fr,f_vec);
		}

	      //mpfr_printf("Newton (right) returned %.40Re\n",znew_guess);
	      bracket_t(zlow,zhigh,ztn,znew_guess);
	      mpfi_inter_t(ztmp,f_vec,zlow);
	      res1_sign=sign(ztmp);
	      mpfi_inter_t(ztmp,f_vec,zhigh);
	      res2_sign=sign(ztmp);
	      if((res1_sign&res2_sign)!=0)
		{
		  binary_chop(znew_guess,tm,tr,fm,fr,f_vec);
		  bracket_t(zlow,zhigh,ztn,znew_guess);
		  mpfi_inter_t(ztmp,f_vec,zlow);
		  res1_sign=sign(ztmp);
		  mpfi_inter_t(ztmp,f_vec,zhigh);
		  res2_sign=sign(ztmp);
		  if((res1_sign&res2_sign)!=0) // still not a valid sign change
		    {
		      printf("Failed to resolve sign (stat) even after binary chop.\n");
		      return(count);
		    }
		}
	      
	      mpfr_sub(last_t,ztn,last_t,GMP_RNDN);
	      if(mpfr_sgn(last_t)<=0)
		{
		  printf("Co-incident or switched zeros.\n");
		  return(count);
		}
	      */
	      count++;
	      //out_bytes(last_t,outfile);
	      //mpfr_set(last_t,ztn,GMP_RNDN);
	      //mpfr_printf("zero found near %.50Re\n",ztn);
	    }
	}
    }
  return(count);
}

// int Im loggamma(1/4+it/2)
void im_int1(mpfi_ptr res, mpfi_ptr t)
{
  mpfi_mul_ui(im_t,t,4); // 4t
  mpfi_atan1(res,im_t);
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
  
void im_int(mpfi_ptr res, mpfi_ptr t0, mpfi_ptr t1)
{
  mpfi_sub(im_err,t1,t0);
  mpfi_div_ui(im_err,im_err,4);
  mpfi_div(im_err,im_err,t0);
  mpfi_neg(im_2,im_err);
  mpfi_put(im_err,im_2);
  mpfi_div_ui(im_t0,t0,2);
  mpfi_div_ui(im_t1,t1,2);
  im_int1(res,im_t1);
  im_int1(im_2,im_t0);
  mpfi_sub(res,res,im_2);
  mpfi_mul_ui(res,res,2);
  mpfi_add(res,res,im_err);
}


double Nleft_int(long int t0_ptr, long int t1_ptr, mpfi_c_t *f_vec, double delta)
{
  long int res=0;
  long int ptr=t0_ptr,last_ptr;
  sign_t last_sign=sign(f_vec[ptr++]->re),this_sign;
  while(last_sign==UNK)
    last_sign=sign(f_vec[ptr++]->re);
  last_ptr=ptr-1;
  for(;ptr<=t1_ptr;)
    {
      this_sign=sign(f_vec[ptr++]->re);
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

double Nright_int(long int t0_ptr, long int t1_ptr, mpfi_c_t *f_vec, double delta)
{
  long int res=0;
  long int ptr=t0_ptr;
  sign_t last_sign=sign(f_vec[ptr++]->re),this_sign;
  while(last_sign==UNK)
    last_sign=sign(f_vec[ptr++]->re);
  for(;ptr<=t1_ptr;)
    {
      this_sign=sign(f_vec[ptr++]->re);
      if(this_sign==last_sign) // no sign change here, move on
	  continue;
      if(this_sign==UNK) // might be a sign change coming
	continue;
      // definately a sign change
      //printf("Sign change in Nright\n");
      res+=(t1_ptr-ptr+1);
      last_sign=this_sign;
    }
  //printf("Nright_int returning %f\n",res*delta);
  return(res*delta);
}

// returns int_t0^{t0+h} S(t) dt
// t0=t0_ptr*delta
// t0>168*Pi
// Trudgian Thesis Theorem 5.2.2
void St_int(mpfi_ptr res, mpfi_ptr t) 
{
  //assert(mpfi_cmp_d(t,168.0*M_PI+0.1)>0);
  mpfi_log(res,t);
  mpfi_mul_d(res,res,0.0590001); // ensures > 0.059,2.067 resp
  mpfi_add_d(res,res,2.0670001);
  //mpfi_neg(tm_tmp,res);
  //mpfi_put(res,tm_tmp);
}

// the log term integrated
void ln_term(mpfi_ptr res, mpfi_ptr t0, mpfi_ptr t1, mpfi_ptr h1)
{
  mpfi_add(res,t0,t1);
  mpfi_mul_d(res,res,-0.25);
  mpfi_mul(res,res,mpfi_ln_pi);
  mpfi_mul(res,res,h1);
}

// the maximum number of zeros <=a based on Turing region [a,b]
long int turing_max(mpfi_c_t *f_vec, long int a_ptr, double a, long int b_ptr, double b)
{
  //printf("In turing max\n");
  mpfi_set_d(tm_t0,a);
  //printf("t0=");mpfi_printn(tm_t0,10);
  mpfi_set_d(tm_t1,b);
  //printf("t1=");mpfi_printn(tm_t1,10);  
  mpfi_sub(tm_h,tm_t1,tm_t0);
  //  
  St_int(tm_res,tm_t1);
  //printf("int S(t)=");mpfi_printn(tm_res,10);
  mpfi_sub_d(tm_res,tm_res,Nright_int(a_ptr,b_ptr,f_vec,one_over_A));
  //printf("int S(t) - int Nright");mpfi_printn(tm_res,10);
  ln_term(tm_tmp,tm_t0,tm_t1,tm_h);
  //printf("int - t log(pi)/2=");mpfi_printn(tm_tmp,10);
  im_int(tm_tmp1,tm_t0,tm_t1);
  //printf("Im lnGamma term=");mpfi_printn(tm_tmp1,10);
  mpfi_add(tm_tmp,tm_tmp,tm_tmp1);
  mpfi_div(tm_tmp,tm_tmp,mpfi_pi);
  mpfi_add(tm_res,tm_res,tm_tmp);
  mpfi_div(tm_res,tm_res,tm_h);
  //printf("turing_max computed ");mpfi_printn(tm_res,10);
  mpfi_get_right(tm_mpfr,tm_res);
  return(mpfr_get_si(tm_mpfr,GMP_RNDD)+1); // +1 because of sqrt(omega) =+/-1
}

long int turing_min(mpfi_c_t *f_vec, long int a_ptr, double a, long int b_ptr, double b)
{
  //printf("In Turing Min\n");
  mpfi_set_d(tm_t0,a);
  //printf("t0=");mpfi_printn(tm_t0,10);  
  mpfi_set_d(tm_t1,b);;
  //printf("t1=");mpfi_printn(tm_t1,10); 
  mpfi_sub(tm_h,tm_t1,tm_t0);
  //  
  St_int(tm_res,tm_t1);
  mpfi_neg(tm_res,tm_res);
  //printf("int S(t)=");mpfi_printn(tm_res,10);
  mpfi_sub_d(tm_res,tm_res,Nleft_int(a_ptr,b_ptr,f_vec,one_over_A));
  //printf("int S(t) - int Nright");mpfi_printn(tm_res,10);
  ln_term(tm_tmp,tm_t0,tm_t1,tm_h);
  //printf("int - t log(pi)/2=");mpfi_printn(tm_tmp,10);
  im_int(tm_tmp1,tm_t0,tm_t1);
  //printf("Im lnGamma term=");mpfi_printn(tm_tmp1,10);
  mpfi_add(tm_tmp,tm_tmp,tm_tmp1);
  mpfi_div(tm_tmp,tm_tmp,mpfi_pi);
  mpfi_add(tm_res,tm_res,tm_tmp);
  mpfi_div(tm_res,tm_res,tm_h);
  //printf("turing_max computed ");mpfi_printn(tm_res,10);
  mpfi_get_left(tm_mpfr,tm_res);
  return(mpfr_get_si(tm_mpfr,GMP_RNDU)+1); // +1 because of sqrt(omega) =+/-1
}

// Use Turing's method to estimate and then check number of zeros
long int turing(mpfi_c_t *f_vec, long int a_ptr, double a, long int b_ptr, double b, long int last_max, FILE *outfile, bool stat_pts)
{
  // etimate maximum for N(b)

  long int i;

  long int min_lo,max_hi=turing_max(f_vec,b_ptr,b,b_ptr+TURING_LEN,b+TURING_WIDTH);
  if(last_max==0) // No previous run or previous run failed, so estimate min N(a)
    min_lo=turing_min(f_vec,a_ptr-TURING_LEN,a-TURING_WIDTH,a_ptr,a);
  else // If previous run succeeded, zeros to height N(a) is known
    min_lo=last_max;
  printf("looking for %ld-%ld=%ld zeros\n",max_hi,min_lo,max_hi-min_lo);
  long int num_found,num_exp=max_hi-min_lo;
  fwrite(&a,sizeof(double),1,outfile);
  fwrite(&b,sizeof(double),1,outfile);
  sign_t ls,rs;
  for(i=a_ptr,ls=sign(f_vec[i]->re);ls==UNK;i++);
  for(i=b_ptr,rs=sign(f_vec[i]->re);rs==UNK;i++);
  i=0;
  if(ls==rs)
    {
      if(num_exp&1)
	{
	  printf("Problem in Turing Zone.\n");
	  printf("Missed All/All zeros in region %f to %f.\n",a,b);
	  fwrite(&i,sizeof(int),1,outfile);
	  return(0);
	}
    }
  else
    {
      if(!(num_exp&1))
	{
	  printf("Problem in Turing Zone.\n");
	  printf("Missed All/All zeros in region %f to %f.\n",a,b);
	  fwrite(&i,sizeof(int),1,outfile);
	  return(0);
	}
    }
  fwrite(&min_lo,sizeof(long int),1,outfile); // number zeros <= a
  fwrite(&max_hi,sizeof(long int),1,outfile);  // number zeros <= b
  if(stat_pts)
    num_found=zeros_st(f_vec,a_ptr,b_ptr,outfile);
  else
    num_found=zeros(f_vec,a_ptr,b_ptr,outfile);
  if(num_exp==num_found)
    {
      if(stat_pts)
	printf("All %ld zeros found in region %f to %f using stat points.\n",num_exp,a,b);
      else
	printf("All %ld zeros found in region %f to %f without using stat points.\n",num_exp,a,b);
      return(max_hi); // use this for next iteration
    }
  if(stat_pts)
    {
      printf("Missed %ld/%ld zeros (running with stat points) in region %f to %f.\n",num_exp-num_found,num_exp,a,b);
      return(0);
    }
  // lets see if checking stats will help
  long int stats=num_stat_pts(f_vec,a_ptr,b_ptr);
  if(num_exp==(num_found+stats))
    {
      printf("All %ld zeros will be found (with %ld stats) in region %f to %f.\n",num_exp,stats,a,b);
      return(0);
    }
  printf("Will miss %ld/%ld zeros in region %f to %f even after %ld stats.\n",num_exp-num_found-stats,num_exp,a,b,stats);
  return(0);
}

int main(int argc, char **argv)
{
  long int prec,num_its,it;
  long int i,j,k,int_step,last_max=0;
  double t0,t1,step,t;
  FILE *outfile;
  fpos_t fpos;
  mpfi_t tmp;
  bool first=TRUE,stat_pts;
  mpfi_c_t omega;
  time_t last_time=time(NULL),this_time;
  if(argc!=7)
    print_usage();

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    print_usage();
  //printf("Running at %d bits of precision.\n",prec);
    
  mpfi_c_setup(prec);
  mpfi_init(tmp);

  mpfi_c_init(omega);
  mpfi_div_ui(omega->re,mpfi_2_pi,N);
  mpfi_sin(omega->im,omega->re);
  mpfi_cos(omega->re,omega->re);

  outfile=fopen(argv[5],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[5]);
      exit(1);
    }

  t0=atof(argv[2]);
  num_its=atoi(argv[3]);
  if(num_its<=0)
    print_usage();
  step=atof(argv[4]);
  if(step<=0)
    print_usage();
  t0=t0+step/2.0;
  assert(t0<=T0_MAX);
  assert(t0>=T0_MIN);
  assert(prec>OP_ACC+1+log(T0_MAX+step/2.0)/log(2.0)); // enough bits to make 2^(OP_ACC+1)*t0 exact ?
  HI_PREC=prec+EXTRA_PREC;
  int_step=step/one_over_A;
  assert((int_step&1)==0); // must be even
  assert(step/one_over_A-(double) int_step ==0.0); // must be an exact number of steps
  //printf("Aiming to get %d values per run.\n",int_step+1);
  if(atoi(argv[6]))
    {
      printf("Running With Stationary Points.\n");
      stat_pts=TRUE;
    }
  else
    {
      printf("Running without stat points.\n");
      exit(0);
      stat_pts=FALSE;
    }
  exps1=(mpfi_t *) malloc((int_step/2+Ns*INTER_SPACING)*sizeof(mpfi_t));
  if(!exps1)
    {
      printf("Failed to allocate memory for exps1. Exiting.\n");
      exit(0);
    }
  for(i=0,t=one_over_A;i<int_step/2+Ns*INTER_SPACING;t+=one_over_A,i++)
    {
      mpfi_init(exps1[i]);
      mpfi_set_d(exps1[i],t);
      mpfi_div_d(exps1[i],exps1[i],h);
      mpfi_sqr(exps1[i],exps1[i]);
      mpfi_div_ui(exps1[i],exps1[i],2);
      mpfi_exp(exps1[i],exps1[i]);
    }

  printf("Aiming to do %ld iteration(s) starting at t=%13.11e.\n",num_its,t0-step/2.0);
  fwrite(&num_its,sizeof(long int),1,outfile);
  for(it=0;it<num_its;it++,t0+=step)
    {
      //printf("Running centred at t0=%f.\n",t0);
      if(first)
	{
	  first=FALSE;
	  //printf("Calling init.\n");
	  //system("date");
	  init(t0);
	  //printf("Init finished.\n");
	  //system("date");
	}
      else
	{
	  //printf("Calling re-init.\n");
	  //system("date");
	  re_init(t0);
	  //printf("re-init finished.\n");
	  //system("date");
	}

      this_time=time(NULL);
      printf("Time to (re-)initialise = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;

      for(k=0;k<K;k++)
	{
	  G_k(k);
	  do_conv(k);

	}
      this_time=time(NULL);
      printf("Time to convolve = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;

      //printf("convolutions finished\n");
      //system("date");

      // upsample by "zero" padding
      for(i=N1/2;i<=N/2;i++)
	mpfi_c_set(f_vec[i],Fmaxerr);

      for(i=0;i<=N/2;i++)
	{
	  mpfi_c_add(f_vec[i],f_vec[i],fhatsumerr);
	  mpfi_c_add(f_vec[i],f_vec[i],tayerr);
	  mpfi_c_add(f_vec[i],f_vec[i],fhattwiderr);
	}

      // f is real so use conjugate symmetry
      // could use a real idtf instead
      for(i=N/2+1;i<N;i++)
	mpfi_c_conj(f_vec[i],f_vec[N-i]);

      for(i=1;i<N;i+=2)
	mpfi_c_neg(f_vec[i],f_vec[i]);

      //printf("Final iFFT\n");
      hermidft(f_vec,N/2,ws_r,omega);
      this_time=time(NULL);
      printf("Time for final iFFT = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;

      //printf("Final iFFT finished.\n");
      //system("date");

      // remove the Gaussian from the region of interest
      for(j=0,i=N/2+1;i<=N/2+int_step/2+INTER_SPACING*Ns;j++,i++)
	{
	  mpfi_div_d(f_vec[i]->re,f_vec[i]->re,B);
	  mpfi_add(f_vec[i]->re,f_vec[i]->re,ftwiderr->re);
	  mpfi_mul(f_vec[i]->re,f_vec[i]->re,exps1[j]);
	}
      for(j=0,i=N/2-1;i>=N/2-int_step/2-INTER_SPACING*Ns;j++,i--)
	{
	  mpfi_div_d(f_vec[i]->re,f_vec[i]->re,B);
	  mpfi_add(f_vec[i]->re,f_vec[i]->re,ftwiderr->re);
	  mpfi_mul(f_vec[i]->re,f_vec[i]->re,exps1[j]);
	}
      
      mpfi_div_d(f_vec[N/2]->re,f_vec[N/2]->re,B);
      mpfi_add(f_vec[N/2]->re,f_vec[N/2]->re,ftwiderr->re);
      //for(i=N/2-int_step/2;i<N/2+int_step/2;i+=int_step/10)
      //printf("error in f_vec[%d]=%d\n",i,mpfi_abs_error(f_vec[i]->re));
      if(fgetpos(outfile,&fpos)!=0)
	{
	  printf("Wierd error if fgetpos. Exiting.\n");
	  exit(0);
	}
      last_max=turing(f_vec,N/2-int_step/2,t0-int_step/2*one_over_A,N/2+int_step/2,t0+int_step/2*one_over_A,last_max,outfile,stat_pts);
      if(last_max==0) // something went wrong
	{
	  if(fsetpos(outfile,&fpos)!=0)
	    {
	      printf("Wierd error if fsetpos. Exiting.\n");
	      exit(0);
	    }
	  double zero[2]={0.0,0.0};
	  long int zs=0;
	  fwrite(zero,sizeof(double),2,outfile);
	  fwrite(&zs,sizeof(long int),1,outfile);
	}
	  
      this_time=time(NULL);
      printf("Time to isolate zeros = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;

    }
  fclose(outfile);
}
