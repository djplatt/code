// residue_sums_fft.cpp
// based on arb_windowed/arb_zeta.cpp
// windowed zeta calculator
//
// Created: 20th June 2025
//
// DJ Platt
// University of Bristol
//

#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "mpfr.h"
#include "flint/arb.h"
#include "flint/acb.h"
#include "../../arb_windowed/arb_win_zeta.h"
#include "../../arb_windowed/arb_fft.h"

int HI_PREC;

// parameters for t0<=3*10^10
// see zeta-fft-errors.gp
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
#define interderr_d ((double) 2.9e-39)
#define INTER_A ((double) INTER_SPACING*one_over_A)


// Newton Iteration Parameters
#define NEWTON_ITS (8) // Do this many iterations with Newton Raphson to locate roots

#define step ((double) 2100)
#define int_step ((uint64_t) (step/one_over_A))

arb_t exps1[int_step/2+Ns*INTER_SPACING];

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

void fatal_error(const char* str)
{
  printf("%s. Exiting.\n",str);
  exit(0);
}

arb_t arb_pi,arb_ln_pi,arb_2_pi,A,A1,two_pi_t[N1],s_vec[M],sk_vec[M],exps[N1],sqrts[M],logisqrtpi[M];
acb_t g_vec[N1],G_vec[N1],f_vec[N],skn_vec[N1],ws_r[N/4],ws1_r[N1/2],ws1_f[N1/2],n_vec[M];
acb_t gtwiderr,Gtwiderr,c_tmp,fhatsumerr,fhattwiderr,ftwiderr,tayerr,Fmaxerr;
int buck[M]; // which bucket does log(nsqrt(Pi)/2pi) go in.
arb_t intererr,interderr;
arf_t ip_tmp1,ip_tmp2,ip_tmp3,ip_tmp4;
arf_t two_101,offset;
arb_t msinc_tmp,arb_two_pi_B,mip_tmp,minter_tmp,minter_tmp1;
//arf_t znew_guess,zold_guess,ztn;
arb_t misin;

// note (sign&&sign)==0 <=> signs different and known 

void acb_add_error(acb_t x, acb_t err)
{
  arb_add_error(acb_realref(x),acb_realref(err));
  arb_add_error(acb_imagref(x),acb_imagref(err));
}

void acb_add_error(acb_t x, arb_t err)
{
  arb_add_error(acb_realref(x),err);
  arb_add_error(acb_imagref(x),err);
}

// inits and sets z to [-d,d]+[-d,d]i
void set_err(acb_t z, double d)
{
  acb_init(z);
  arb_set_d(acb_realref(z),d);
  arb_set(acb_imagref(z),acb_realref(z));
}

void set_err(arb_t x, double d)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_init(x);
  arb_set_d(x,d);
}

// reverse skn_vec
inline int conj_order(int i)
{
  return(N1-i);
}

// on entry s_vec[i] contains log((i+1)Pi)/2Pi
// on exit                    "" - u_m
void calc_buck(int64_t prec)
{
  static bool init=false;
  static arb_t tmp,arb_B;
  if(!init)
    {
      arb_init(tmp);
      arb_init(arb_B);
      arb_set_d(arb_B,B);
    }
  int i;
  double x;
  for(i=0;i<M;i++)
    {
      buck[i]=arb_get_d(s_vec[i])*B+0.5;
      arb_set_ui(tmp,buck[i]);
      arb_div(tmp,tmp,arb_B,prec);
      arb_sub(s_vec[i],s_vec[i],tmp,prec);
    }
  if(buck[0]<=0) // keep conj_order simple
    fatal_error("buck[0]<=0.");
  if(buck[M-1]>=N1/2) // ditto
    fatal_error("buck[M-1]>=N/2.");
  for(i=0;i<M;i++)
    buck[i]=conj_order(buck[i]);
}

arb_t atan1_tmp1,atan1_tmp2,im_1,im_2,im_t,im_err;
arb_t tm_h,tm_res,tm_tmp,tm_tmp1,im_t0,im_t1;
//arf_t tm_arf,sp_t,sp_ft,sp_t1,sp_ft1;
arb_t tm_t0,tm_t1,init_tmp1,sqrt_pi,pi_hi,half_log_2_pi;

// gamma(1/4+i(t+t0)/2)*exp(Pi*(t+t0)/4)*exp(-t^2/(2h^2))
void gamma_factor(acb_ptr res, arb_ptr gaussian, double t,int64_t prec)
{
  static bool init=false;
  static acb_t z,tmp;
  static arb_t two_pi,tmp1; // pi/2
  if(!init)
    {
      init=true;
      acb_init(z);
      acb_init(tmp);
      arb_init(tmp1);
      arb_set_d(acb_realref(z),0.25);
      arb_init(two_pi);
      arb_mul_2exp_si(two_pi,arb_pi,-1);
    }
  arb_set_d(acb_imagref(z),t/2.0);
  acb_lgamma(tmp,z,prec); // lng(1/4+it/2)
  arb_mul(tmp1,acb_imagref(z),two_pi,prec); // pi t/4
  arb_add(acb_realref(tmp),acb_realref(tmp),tmp1,prec); // 
  arb_sub(acb_realref(tmp),acb_realref(tmp),gaussian,prec); // 
  acb_exp(res,tmp,prec);
}

// initialise stuff
// on exit g[n]=g((n-N/2)/A;0)
void init(double t0,int64_t prec)
{
  int i;
  double t;
  arf_t arf_tmp;
  arf_init(arf_tmp);
  arb_init(arb_pi);
  arb_const_pi(arb_pi,prec);
  arb_init(arb_ln_pi);
  arb_log(arb_ln_pi,arb_pi,prec);
  arb_init(arb_2_pi);
  arb_mul_2exp_si(arb_2_pi,arb_pi,1);
  arf_init(ip_tmp2);
  arf_init(ip_tmp1);
  arf_init(ip_tmp2);
  arf_init(ip_tmp3);
  arf_init(ip_tmp4);
  arf_init(two_101);
  arf_set_ui(two_101,1);
  arf_mul_2exp_si(two_101,two_101,OP_ACC);
  // we assume converting mpz_t to long unsigned int
  // gives us 64 bits
  if(sizeof(long unsigned int)!=8)
    fatal_error("Expecting long int to be 64 bits.");

  arb_init(A);
  arb_init(A1);
  acb_init(c_tmp);
  arb_init(msinc_tmp);
  arb_init(arb_two_pi_B);
  arb_div_d(arb_two_pi_B,arb_pi,INTER_A,prec);
  arb_init(mip_tmp);
  arb_init(minter_tmp);
  arb_init(minter_tmp1);
  //arb_init(intererr);
  arb_init(init_tmp1);
  arb_init(sqrt_pi);
  arb_init(half_log_2_pi);
  arb_const_pi(sqrt_pi,HI_PREC);
  arb_init(pi_hi);
  arb_set(pi_hi,sqrt_pi);
  arb_set(half_log_2_pi,sqrt_pi);
  arb_mul_2exp_si(half_log_2_pi,half_log_2_pi,1);
  arb_log(half_log_2_pi,half_log_2_pi,HI_PREC);
  arb_mul_2exp_si(half_log_2_pi,half_log_2_pi,-1);
  arb_sub_d(half_log_2_pi,half_log_2_pi,0.25,HI_PREC); // 1/2log(2 Pi)-1/4
  arb_sqrt(sqrt_pi,sqrt_pi,HI_PREC);
  arb_init(misin);

  for(i=0,t=-N1/2*one_over_A1;i<N1;i++,t+=one_over_A1)
    {
      arb_init(two_pi_t[i]);
      arb_mul_d(two_pi_t[i],arb_2_pi,-t,prec);
      acb_init(g_vec[i]);
      acb_init(G_vec[i]);
      acb_init(f_vec[i]);
      acb_init(skn_vec[i]);
      arb_init(exps[i]);
      arb_set_d(exps[i],t);
      arb_div_d(exps[i],exps[i],h,HI_PREC);
      arb_mul(exps[i],exps[i],exps[i],HI_PREC);
      arb_mul_2exp_si(exps[i],exps[i],-1); //  t^2/(2h^2)
      gamma_factor(g_vec[i],exps[i],t+t0,prec);
    }
  for(i=N1;i<N;i++)
    acb_init(f_vec[i]);
  for(i=0;i<M;i++)
    {
      arb_init(s_vec[i]);
      arb_init(sk_vec[i]);
      arb_init(logisqrtpi[i]);
      acb_init(n_vec[i]);
      arb_mul_ui(logisqrtpi[i],sqrt_pi,i+1,HI_PREC); // n*sqrt(Pi)
      arb_log(logisqrtpi[i],logisqrtpi[i],HI_PREC);        // log (n*sqrt(Pi))
      arb_set(s_vec[i],logisqrtpi[i]);
      arb_mul_d(init_tmp1,logisqrtpi[i],-t0,HI_PREC); // -t0*log(n*sqrt(Pi))
      arb_sin_cos(acb_imagref(n_vec[i]),acb_realref(n_vec[i]),init_tmp1,HI_PREC);
      arb_set_ui(A,1);
      arb_div_ui(A,A,i+1,prec);
      arb_init(sqrts[i]);
      arb_sqrt(sqrts[i],A,prec);
      acb_mul_arb(n_vec[i],n_vec[i],sqrts[i],prec);    //  / sqrt(n)
      arb_div(s_vec[i],s_vec[i],arb_2_pi,prec); // log(n*sqrt(Pi))/(2*Pi)
    }

  // set up the omega vectors for FFTs
  initfft(N/2,ws_r,prec); // ws_r[i]=e(i/N)
  for(i=0;i<N1/2;i++)
    {
      acb_init(ws1_f[i]);
      acb_init(ws1_r[i]);
      acb_set(ws1_r[i],ws_r[i*UPSAM/2]); // ws1_r[i]=e(i/N1)
      acb_conj(ws1_f[i],ws1_r[i]);     // ws1_f[i]=e(-i/N1)
    }

  set_err(gtwiderr,gtwiderr_d);
  set_err(Gtwiderr,Gtwiderr_d);
  set_err(tayerr,tayerr_d);
  set_err(Fmaxerr,Fmaxerr_d);
  set_err(intererr,intererr_d);
  set_err(interderr,interderr_d);

  arb_set_ui(A,M);
  arb_sqrt(A,A,prec);
  arb_mul_2exp_si(A,A,1);
  arb_sub_ui(A,A,1,prec); // 2sqrt(M)-1
  acb_mul_arb(tayerr,tayerr,A,prec);
  //printf("Total Taylor Error set to ");mpfi_printn(tayerr->re,30);
  set_err(fhatsumerr,fhatsumerr_d);
  set_err(fhattwiderr,fhattwiderr_d);
  set_err(ftwiderr,ftwiderr_d);
  arb_set_ui(A,1);
  arb_div_d(A,A,one_over_A,prec); // don't use A as a spare variable anymore
  arb_set_ui(A1,1);
  arb_div_d(A1,A1,one_over_A1,prec);
  calc_buck(prec);
  arb_init(atan1_tmp1);
  arb_init(atan1_tmp2);
  arb_init(im_1);
  arb_init(im_2);
  arb_init(im_t);
  arb_init(im_err);
  arb_init(tm_h);
  arb_init(tm_res);
  arb_init(tm_tmp);
  arb_init(tm_tmp1);
  arb_init(im_t0);
  arb_init(im_t1);
  arb_init(tm_t0);
  arb_init(tm_t1);
}

// run with a new value of t0
// just need to set g_vec and n_vec up
// rest was already done by init
void re_init(double t0, int64_t prec)
{
  int i;
  double t;
  for(i=0,t=-N1/2*one_over_A1;i<N1;i++,t+=one_over_A1)
      gamma_factor(g_vec[i],exps[i],t+t0,prec);
  for(i=0;i<M;i++)
    {
      arb_mul_d(init_tmp1,logisqrtpi[i],-t0,HI_PREC);
      arb_sin_cos(acb_imagref(n_vec[i]),acb_realref(n_vec[i]),init_tmp1,HI_PREC);
      acb_mul_arb(n_vec[i],n_vec[i],sqrts[i],prec);    //  / sqrt(i+1)
    }
}

// on entry g_vec[i]=g(t)*(-2*Pi*t*I)^k/k!=g(t;k)/k!
// on exit G_vec[i]=G^(k)(i/B)/k!
void G_k(int k, int64_t prec)
{
  int i;
  double t,max_t;

  for(i=0;i<N1;i++)
    {
      acb_set(G_vec[i],g_vec[i]);
      acb_add_error(G_vec[i],gtwiderr); // G=g~
    }
  
  fft(G_vec,N1,ws1_f,prec);

  for(i=0;i<=N1/2;i++)
    {
      acb_div_arb(G_vec[i],G_vec[i],A1,prec);
      acb_add_error(G_vec[i],Gtwiderr);
      if(i&1)
	acb_neg(G_vec[i],G_vec[i]);
    }

  for(i=N1/2+1;i<N1;i++)
    acb_zero(G_vec[i]);

  if(k<K-1)
    for(i=0;i<N1;i++)
      {
	acb_mul_onei(g_vec[i],g_vec[i]); // g(t;k)*i/k!
	acb_mul_arb(g_vec[i],g_vec[i],two_pi_t[i],prec); // g(t;k+1)/k!
	acb_div_ui(g_vec[i],g_vec[i],k+1,prec);        // g(t;k+1)/(k+1)!
      }
}


// on entry sk_vec = (log(nsqrt(Pi))/2Pi-u_m)^k (or undefined if k==0)
// n_vec[i] = (i+1)^(-1/2)*(sqrt(Pi)*(i+1))^-it0
// on exit skn_vec = sum nsqrt(pi)^(-ito)/sqrt(n)*(log(nsqrt(Pi))/2Pi-u_m)^k
//         sk_vec = (log..-u_m)^(k+1)
void make_skn(int k,int64_t prec)
{
  int i,pos;
  for(i=0;i<N1;i++)
    acb_zero(skn_vec[i]);
  switch(k)
    {
    case 0: // set sk_vec = s_vec, skn_vec=n_vec
      for(i=0;i<M;i++)
	{
	  acb_add(skn_vec[buck[i]],skn_vec[buck[i]],n_vec[i],prec);
	  arb_set(sk_vec[i],s_vec[i]);
	}
      break;
    case K-1: // no need to increment sk_vec
      for(i=0;i<M;i++)
	{
	  acb_mul_arb(c_tmp,n_vec[i],sk_vec[i],prec);
	  acb_add(skn_vec[buck[i]],skn_vec[buck[i]],c_tmp,prec);
	}
      break;
    default:
      for(i=0;i<M;i++)
	{
	  acb_mul_arb(c_tmp,n_vec[i],sk_vec[i],prec);
	  acb_add(skn_vec[buck[i]],skn_vec[buck[i]],c_tmp,prec);
	  arb_mul(sk_vec[i],sk_vec[i],s_vec[i],prec);
	}
    }
}

// convolve v1 with v2 into res
// on exit v1 and v2 contain fft of input
void my_convolve (acb_t *res, acb_t *v1, acb_t *v2, uint64_t n, acb_t *ws_r, acb_t *ws_f,int64_t prec)
{
  int i;
  fft(v1,n,ws_r,prec); 
  fft(v2,n,ws_r,prec);
  for(i=0;i<n;i++)
    acb_mul(res[i],v1[i],v2[i],prec);
  fft(res,n,ws_f,prec); // use ws_f so this is an iFFT
  for(i=0;i<n;i++)
    acb_div_ui(res[i],res[i],n,prec); // normalise the iFFT
}

// f_vec+=G(k)(x+u_m)/k!*S^(k)_m
void do_conv (int k, int64_t prec)
{
  int i;

  make_skn(k,prec);    

  if(k==0) // convolve straight into f_vec
    my_convolve(f_vec,skn_vec,G_vec,N1,ws1_r,ws1_f,prec);
  else // convolve into G_vec, then add into f_vec 
    {
      my_convolve(G_vec,skn_vec,G_vec,N1,ws1_r,ws1_f,prec);
      for(i=0;i<=N1/2;i++)
	acb_add(f_vec[i],f_vec[i],G_vec[i],prec);
    }
}

// compute exp(-t^2/2H^2)
void inter_gaussian(arb_ptr res, arb_ptr t, int64_t prec)
{
  arb_div_d(res,t,H,prec);
  arb_mul(res,res,res,prec);
  arb_div_d(res,res,-2.0,prec);
  arb_exp(res,res,prec);
}

bool sincp;
// compute sinc(2*B*Pi*t) into ress, cos(2*B*Pi*t) into resc
// first time we call this (with sincp true) we compute sin and cos
// from then on (sincp false) just swap signs each time
void inter_sinc_cos(arb_ptr ress, arb_ptr resc, arb_ptr t, int64_t prec)
{
  static bool init=false;
  static arb_t sinc_tmp,msin,mcos;
  if(!init)
    {
      init=true;
      arb_init(sinc_tmp);
      arb_init(msin);
      arb_init(mcos);
    }
  arb_mul(sinc_tmp,t,arb_two_pi_B,prec);
  if(sincp)
    {
      arb_sin_cos(msin,mcos,sinc_tmp,prec);
      sincp=false;
    }
  else
    {
      arb_neg(msin,msin);
      arb_neg(mcos,mcos);
    }
  arb_set(resc,mcos);
  arb_div(ress,msin,sinc_tmp,prec);
}

#define debug (false)

// t is distance from t0 to t implied by f_ptr
// if fd_res is NULL, don't calc differential
void inter_point(arb_ptr f_res, arb_ptr fd_res, acb_t *f_vec, arb_ptr t, int f_ptr, int64_t prec)
{
  static bool init=false;
  static arb_t ip_tmp1,ip_tmp2,ip_tmp3,ip_tmp4;
  if(!init)
    {
      init=true;
      arb_init(ip_tmp1);
      arb_init(ip_tmp2);
      arb_init(ip_tmp3);
      arb_init(ip_tmp4);
    }
  arb_set(ip_tmp1,acb_realref(f_vec[f_ptr]));
  inter_gaussian(ip_tmp2,t,prec);
  inter_sinc_cos(ip_tmp3,ip_tmp4,t,prec);
  if(debug)
    {
      printf("Lambda(t) = ");arb_printd(ip_tmp1,20);
      printf("\nGuassian = ");arb_printd(ip_tmp2,20);
      printf("\nsinc = ");arb_printd(ip_tmp3,20);
    }
  // compute f_res
  arb_mul(ip_tmp2,ip_tmp1,ip_tmp2,prec); //f*exp
  arb_mul(f_res,ip_tmp2,ip_tmp3,prec); //f*exp*sinc
  if(debug)
    {
      printf("\n Inter_point returning ");
      arb_printd(f_res,20);
      printf("\n");
    }
  // compute fd_res
  if(!fd_res)
    return;
  arb_div(fd_res,ip_tmp4,t,prec); // cos/t
  arb_set_ui(ip_tmp4,1);
  arb_div(ip_tmp4,ip_tmp4,t,prec); // 1/t
  arb_div_d(ip_tmp1,t,H*H,prec); // t/H^2
  arb_add(ip_tmp4,ip_tmp4,ip_tmp1,prec); // 1/t+t/H^2
  arb_mul(ip_tmp4,ip_tmp4,ip_tmp3,prec); // sinc*(1/t+t/H^2)
  arb_sub(fd_res,fd_res,ip_tmp4,prec);   // cos/t-sinc*(1/t+t/H^2)
  arb_mul(fd_res,fd_res,ip_tmp2,prec);   // f*exp*(cos/t-sinc....)
  //
  // WARNING. COMPLETELY ARBITARY NEGATION TO MAKE f'(t) CORRECT SIGN
  //
  // arb_neg(fd_res,fd_res);
}

// t0 points into f_vec
// let t=t0+(t_ptr-N/2)*one_over_A
// returns f(t) in f_res, f'(t) in fd_res
// return f_vec[t0]->re if t_ptr is integral 
void inter_t(arb_ptr f_res, arb_ptr fd_res, acb_t *f_vec, arb_ptr t_ptr, int64_t prec)
{
  static bool init=false;
  static arb_t inter_tmp,inter_tmp1,inter_tmp2;
  if(!init)
    {
      init=true;
      arb_init(inter_tmp);
      arb_init(inter_tmp1);
      arb_init(inter_tmp2);
    }
  int i,j,nearest_t0;
  if(debug) printf("In inter_t\n");
  double t0_d=arb_get_d(t_ptr),dt0;
  arb_set_si(inter_tmp,t0_d);
  if(arb_equal(t_ptr,inter_tmp)) // t_ptr is integral so = a lattice pt
    nearest_t0=t0_d+1.5;
  else
    nearest_t0=t0_d+0.5;
  if(debug) printf("Nearest t0=%ld\n",nearest_t0);
  arb_zero(f_res);
  if (fd_res) arb_zero(fd_res);
  sincp=true;
  if(nearest_t0>t0_d) // start from right of t0
    {
      if(debug) printf("Starting from right.\n");
      for(i=nearest_t0,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if(debug)
	    {
	      printf("i=%d\nptr=",i);
	      arb_printd(inter_tmp,20);
	      printf("\nPoint returned ");
	      arb_printd(inter_tmp1,20);
	      printf("\nsum now ");
	      arb_printd(f_res,20);
	      printf("\n");
	    }
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
      sincp=true;
      for(i=nearest_t0-INTER_SPACING,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if(debug)
	    {
	      printf("i=%d\nptr=",i);
	      arb_printd(inter_tmp,20);
	      printf("\nPoint returned ");
	      arb_printd(inter_tmp1,20);
	      printf("\nsum now ");
	      arb_printd(f_res,20);
	      printf("\n");
	    }
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
    }
  else // start from left
    {
      if(debug) printf("Starting from left.\n");
      for(i=nearest_t0+INTER_SPACING,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  //printf("Doing right tail with i=%d\n",i);
	  //
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  //printf("Calling ip\n");
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  //arf_printf("ip returned %.40Re\n",inter_tmp1);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if(debug)
	    {
	      printf("i=%d\nptr=",i);
	      arb_printd(inter_tmp,20);
	      printf("\nPoint returned ");
	      arb_printd(inter_tmp1,20);
	      printf("\nsum now ");
	      arb_printd(f_res,20);
	      printf("\n");
	    }
	  //
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	  //
	}
      sincp=true;
      for(i=nearest_t0,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if(debug)
	    {
	      printf("i=%d\nptr=",i);
	      arb_printd(inter_tmp,20);
	      printf("\nPoint returned ");
	      arb_printd(inter_tmp1,20);
	      printf("\nsum now ");
	      arb_printd(f_res,20);
	      printf("\n");
	    }
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
    }
  arb_add_error(f_res,intererr);
  if (fd_res) arb_add_error(fd_res,interderr);
}

void inter_t(arb_ptr f_res, acb_t *f_vec, arb_ptr t_ptr, int64_t prec)
{
  inter_t(f_res, NULL, f_vec, t_ptr, prec);
}

// convert t into a ptr into f_res
// t0 is at N/2
// so t is at N/2+A(t-t0)
void t_to_ptr(arb_ptr ptr, arb_ptr t, double t0, int64_t prec)
{
  arb_set_d(ptr,t0);
  arb_sub(ptr,ptr,t,prec); // t0-t
  arb_mul(ptr,ptr,A,prec); // A(t0-t)
  arb_neg(ptr,ptr); // A(t-t0)
  arb_add_ui(ptr,ptr,N/2,prec); // N/2+A(t-t0)
}

/*
// NB destroys old_guess
// on entry left and right bracket the zero
// on exit new_guess = improved estimate
void newton(arb_ptr new_guess, arb_ptr old_guess, acb_t *f_vec, int num_its, int64_t prec)
{
  static bool init=false;
  static arb_t newton_df;
  if(!init)
    {
      init=true;
      arb_init(newton_df);
    }
  int i;
  for(i=0;i<num_its;i++)
    {
      inter_t(new_guess,newton_df,f_vec,old_guess,prec); // new=f(t) df=f'(t)
      arb_div(newton_df,new_guess,newton_df,prec);  // f/f'
      arb_div_d(newton_df,newton_df,one_over_A,prec); // f/f' normalised
      arb_sub(old_guess,old_guess,newton_df,prec); // old-f/f'
    }
  arb_set(new_guess,old_guess);
}
*/

// read a 13 byte number from file
// structured 8,4,1
// read as if its exact
void in_bytes_exact(arb_ptr t, FILE *infile, int64_t prec)
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


// read the delta between zeros into del_t
inline void next_rho(arb_t del_t, FILE *infile, int64_t prec)
{
  in_bytes_exact(del_t,infile,prec);
}

void lam_to_zeta(arb_t res, arb_t lam, arb_t gamma, int64_t prec)
{
  static acb_t s_by_2,ctmp;
  static arb_t pi,tmp;
  static bool init=false;
  if(!init)
    {
      init=true;
      acb_init(s_by_2);
      acb_init(ctmp);
      arb_set_d(acb_realref(s_by_2),0.25);
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(tmp);
    }
  arb_mul_2exp_si(acb_imagref(s_by_2),gamma,-1); // s_by_2 = 1/4 +it/2
  acb_lgamma(ctmp,s_by_2,prec);
  arb_mul(tmp,pi,gamma,prec);
  arb_mul_2exp_si(tmp,tmp,-2); // pi t/4
  arb_add(acb_imagref(ctmp),tmp,acb_realref(ctmp),prec);
  arb_exp(tmp,acb_imagref(ctmp),prec);
  arb_div(res,lam,tmp,prec);
}

void varphi(acb_t res, acb_t x, int64_t prec)
{
  static bool init=false;
  static arb_t pi2;
  static acb_t ctmp1,ctmp2;
  
  if(!init)
    {
      init=true;
      arb_init(pi2);
      arb_const_pi(pi2,prec);
      arb_mul_2exp_si(pi2,pi2,-1); // pi/2
      acb_init(ctmp1);acb_init(ctmp2);
    }

  acb_mul_arb(ctmp1,x,pi2,prec);
  acb_cot(ctmp2,ctmp1,prec); // cot(x pi/2)
  acb_mul(res,ctmp2,ctmp1,prec); // x pi/2 cot(x pi/2)
}

arb_t max_zetap,min_zetap,max_gamma,min_gamma;

void max_min_zetap(arb_t zetap, arb_t gamma, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(max_zetap); // init to 0
      arb_init(min_zetap);
      arb_set_ui(min_zetap,UWORD_MAX); // very large
      arb_init(max_gamma);
      arb_init(min_gamma);
    }
  arb_sub(tmp,max_zetap,zetap,prec);
  if(arb_is_negative(tmp))
    {
      arb_set(max_zetap,zetap);
      arb_set(max_gamma,gamma);
      return;
    }
  else
    {
      if(!arb_is_positive(tmp))
	{
	  printf("Warning. zetap = ");arb_printd(zetap,20);
	  printf(" and max_zetap = ");arb_printd(max_zetap,20);
	  printf(" overlap at gamma = ");arb_printd(gamma,20);printf("\n");
	}
    }
  arb_sub(tmp,zetap,min_zetap,prec);
  if(arb_is_negative(tmp))
    {
      arb_set(min_zetap,zetap);
      arb_set(min_gamma,gamma);
    }
  else
    {
      if(!arb_is_positive(tmp))
	{
	  printf("Warning. zetap = ");arb_printd(zetap,20);
	  printf(" and min_zetap = ");arb_printd(min_zetap,20);
	  printf(" overlap at gamma = ");arb_printd(gamma,20);printf("\n");
	}
    }
}  
  
arb_t sum1,sum2,sum3,sum4,sum5,sum6;//,sum1a,sum2a,sum3a,sum1b,sum2b,sum3b;

// 1/2+i gamma = rho
// zetap = zeta'(rho)
void do_rho(arb_t gamma, arb_t zetap, arb_t t0, int64_t prec)
{
  static acb_t s,z,ctmp1;
  static arb_t tmp1,inv_s_abs,inv_zetap,den,varp,t1,t2;
  static bool init=false;
  if(!init)
    {
      init=true;acb_init(s);acb_init(z);acb_init(ctmp1);
      arb_set_d(acb_realref(s),0.5);
      arb_init(tmp1);arb_init(inv_s_abs);arb_init(den);arb_init(varp);
      arb_init(inv_zetap);arb_init(t1);arb_init(t2);
      arb_mul_2exp_si(t1,t0,-1);
      arb_div_ui(t2,t0,10,prec);
      arb_init(sum1);
      arb_init(sum2);
      arb_init(sum3);
      arb_init(sum4);
      arb_init(sum5);
      arb_init(sum6);
    }


  max_min_zetap(zetap,gamma,prec);
  
  arb_set(acb_imagref(s),gamma);
  acb_abs(tmp1,s,prec);
  arb_inv(inv_s_abs,tmp1,prec); // |1/rho|
  
  arb_inv(inv_zetap,zetap,prec); // 1/|zeta'|
  arb_add(sum4,sum4,inv_zetap,prec);
  arb_mul(den,inv_zetap,inv_s_abs,prec); // 1/|zeta'(rho) rho|
  arb_add(sum5,sum5,den,prec);
  arb_mul(tmp1,den,inv_s_abs,prec); // 1/|zeta'(rho) rho^2|
  arb_add(sum6,sum6,tmp1,prec);
  
  acb_sub_ui(ctmp1,s,1,prec); // rho-1
  acb_div_onei(ctmp1,ctmp1); // (rho-1)/i
  acb_div_arb(z,ctmp1,t0,prec); // (rho-1)/(t0 i)
  
  // first with t0 into sum1
  varphi(ctmp1,z,prec);
  acb_abs(varp,ctmp1,prec); // |varphi(z)|
  arb_mul(tmp1,varp,den,prec); // 1/|zeta' rho|
  arb_add(sum1,sum1,tmp1,prec);

  arb_sub(tmp1,gamma,t1,prec);
  if(arb_is_positive(tmp1))
    return;

  // now with t1 into sum2
  acb_mul_2exp_si(z,z,1); // (rho-1)/(t0/2 i)
  varphi(ctmp1,z,prec);
  acb_abs(varp,ctmp1,prec); // |varphi(z)|
  arb_mul(tmp1,varp,den,prec); // 1/|zeta' rho|
  arb_add(sum2,sum2,tmp1,prec);

  arb_sub(tmp1,gamma,t2,prec);
  if(arb_is_positive(tmp1))
    return;
  acb_mul_ui(z,z,5,prec); // (rho-1)/(t0/10 i)
  varphi(ctmp1,z,prec);
  acb_abs(varp,ctmp1,prec); // |varphi(z)|
  arb_mul(tmp1,varp,den,prec); // 1/|zeta' rho|
  arb_add(sum3,sum3,tmp1,prec);

  
}


int main(int argc, char **argv)
{
  int rval;
  long int prec,num_its,it;
  long int i,j,k,last_max=0;
  double t0,t1,t;
  fpos_t fpos;
  arb_t tmp;
  bool first=true;
  acb_t omega;
  printf("Command line:-");
  for(uint64_t i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  fflush(stdout);
  if(argc!=4)
    {
      printf("Usage:- %s <zeros file> <T0> <prec>\nExiting.\n",argv[0]);
      return 0;
    }

  FILE *zfile=fopen(argv[1],"rb");
  if(zfile==NULL)
    {
      printf("Failed to open %s for binary read. Exiting.\n",argv[1]);
      return 0;
    }

  double T0=atof(argv[2]);
  prec=atoi(argv[3]);

  // the three sums are to T0, T0/2, T0/10
  arb_t arb_T0,arb_T1,arb_T2;
  arb_init(arb_T0);
  arb_set_d(arb_T0,T0);
  

  if(fread(&num_its,sizeof(long int),1,zfile)!=1)
    {
      printf("Error reading number of iterations from %s. Exiting.\n",argv[1]);
      return 0;
    }

  double st[2];
  long int zs[2];
  arb_t gamma,a_t,del_t,res,resd,ptr,pm1;
  arb_init(gamma);arb_init(a_t);arb_init(del_t);
  arb_init(res);arb_init(ptr);arb_init(pm1);
  arb_init(resd);
  for(i=0,t=one_over_A;i<int_step/2+Ns*INTER_SPACING;t+=one_over_A,i++)
    {
      arb_init(exps1[i]);
      arb_set_d(exps1[i],t);
      arb_div_d(exps1[i],exps1[i],h,prec);
      arb_mul(exps1[i],exps1[i],exps1[i],prec);
      arb_mul_2exp_si(exps1[i],exps1[i],-1);
      arb_exp(exps1[i],exps1[i],prec);
    }

  bool done=false;
  for(long int it=0;it<num_its;it++)
    {
      if(done)
	break;
      rval=fread(st,sizeof(double),2,zfile); // starting/ending t, exact
      rval=fread(zs,sizeof(long int),1,zfile); // starting zero number
      if(st[0]==0.0)
        {
          printf("Iteration %lu empty.\n",it);
          continue;
        }

      
      rval=fread(zs+1,sizeof(long int),1,zfile); // ending zero number
      //printf("processing zeros %ld to %ld inclusive\n",zs[0]+1,zs[1]);
      //printf("doing t from %f to %f zeros from %ld to %ld\n",st[0],st[1],zs[0],zs[1]);fflush(stdout);
      t0=st[0];
      if(t0<T0_MIN)
	fatal_error("t0 outside limits.");
      double this_step=st[1]-st[0];
      t0=t0+this_step/2.0;
      if(t0>T0_MAX)
	fatal_error("t0 outside limits.");
	


      HI_PREC=prec+EXTRA_PREC;
      
      //printf("Running centred at t0=%f.\n",t0);fflush(stdout);
      if(first)
	{
	  first=false;
	  init(t0,prec);
	  arb_init(tmp);
	  acb_init(omega);
	  arb_div_ui(tmp,arb_2_pi,N,prec);
	  arb_sin_cos(acb_imagref(omega),acb_realref(omega),tmp,prec);
	  arb_clear(tmp);
	}
      else
	{
	  re_init(t0,prec);
	}

      // do K convolutions
      for(k=0;k<K;k++)
	{
	  G_k(k,prec);
	  do_conv(k,prec);

	}

      // upsample by "zero" padding
      for(i=N1/2;i<=N/2;i++)
	{
	  acb_zero(f_vec[i]);
	  acb_add_error(f_vec[i],Fmaxerr);
	}

      for(i=0;i<=N/2;i++)
	{
	  acb_add_error(f_vec[i],fhatsumerr);
	  acb_add_error(f_vec[i],tayerr);
	  acb_add_error(f_vec[i],fhattwiderr);
	}

      // f is real so use conjugate symmetry
      // could use a real idtf instead
      for(i=N/2+1;i<N;i++)
	acb_conj(f_vec[i],f_vec[N-i]);

      for(i=1;i<N;i+=2)
	acb_neg(f_vec[i],f_vec[i]);

      //Final iFFT
      hermidft(f_vec,N/2,ws_r,omega,prec);

      // remove the Gaussian from the region of interest
      for(j=0,i=N/2+1;i<=N/2+int_step/2+INTER_SPACING*Ns;j++,i++)
	{
	  arb_div_d(acb_realref(f_vec[i]),acb_realref(f_vec[i]),B,prec);
	  arb_add_error(acb_realref(f_vec[i]),acb_realref(ftwiderr));
	  arb_mul(acb_realref(f_vec[i]),acb_realref(f_vec[i]),exps1[j],prec);
	}
      for(j=0,i=N/2-1;i>=N/2-int_step/2-INTER_SPACING*Ns;j++,i--)
	{
	  arb_div_d(acb_realref(f_vec[i]),acb_realref(f_vec[i]),B,prec);
	  arb_add_error(acb_realref(f_vec[i]),acb_realref(ftwiderr));
	  arb_mul(acb_realref(f_vec[i]),acb_realref(f_vec[i]),exps1[j],prec);
	}
      
      arb_div_d(acb_realref(f_vec[N/2]),acb_realref(f_vec[N/2]),B,prec);
      arb_add_error(acb_realref(f_vec[N/2]),acb_realref(ftwiderr));

      // f_vec now contains Lam(t) = Pi^(-it/2) Gamma(1/4+it/2) zeta(1/2+it)
      // with t =n/A+t0, n=-N/2..N/2-1
      // so f_vec[N/2] contains Lam(t0)

      arb_set_d(a_t,st[0]);
      for(long int z=zs[0]+1;z<=zs[1];z++)
	{
	  next_rho(del_t,zfile,prec);
	  if(arb_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      return 0;
	    }
	  arb_add(a_t,a_t,del_t,prec); // exact
	  if(!arb_is_exact(a_t))
	    {
	      printf("Insufficient precision to hold zeros to +/-2^-%u. Exiting.\n",OP_ACC+1);
	      printf("\nzero centred at ");arb_printd(a_t,20);printf("\n");
	      return 0;
	    }
	  arb_set(gamma,a_t);
	  arb_add_error_2exp_si(gamma,-1-OP_ACC);
	  arb_sub(res,gamma,arb_T0,prec);
	  if(arb_is_positive(res))
	    {
	      done=true;
	      break;
	    }
	  t_to_ptr(ptr,gamma,t0,prec);
	  inter_t(res,resd,f_vec,ptr,prec);
	  if(!arb_contains_zero(res))
	    {
	      printf("Interpolation returned non-zero. Exiting.\n");
	      printf("Zero %ld at ",z);arb_printd(gamma,70);
	      printf("\npointer set to ");arb_printd(ptr,20);
	      printf("\nLambda = ");arb_printd(resd,20);printf("\n");
	      return 0;
	    }
	  lam_to_zeta(res,resd,gamma,prec);
	  arb_abs(res,res);
	  do_rho(gamma,res,arb_T0,prec);
	  //printf("\nzeta'(t) ");arb_printd(res,20);printf("\n");
	}

    }

  fclose(zfile);

  printf("sum (t0=%f) ",(double) T0);arb_printd(sum1,20);
  printf("\nsum (t0=%f) ",(double) T0 / 2.0);arb_printd(sum2,20);
  printf("\nsum (t0=%f) ",(double) T0 / 10.0);arb_printd(sum3,20);
  printf("\nsum 1/|zeta'| ");arb_printd(sum4,20);
  printf("\nsum 1/|rho zeta'| ");arb_printd(sum5,20);
  printf("\nsum 1/|rho^2 zeta'| ");arb_printd(sum6,20);
  printf("\nMax zeta' = ");arb_printd(max_zetap,20);
  printf(" seen at ");arb_printd(max_gamma,20);printf("\n");
  printf("Min zeta' = ");arb_printd(min_zetap,20);
  printf(" seen at ");arb_printd(min_gamma,20);printf("\n");
  printf("\n");

  arb_dump_file(stdout,sum1);printf("\n");
  arb_dump_file(stdout,sum2);printf("\n");
  arb_dump_file(stdout,sum3);printf("\n");
  arb_dump_file(stdout,sum4);printf("\n");
  arb_dump_file(stdout,sum5);printf("\n");
  arb_dump_file(stdout,sum6);printf("\n");
  
  
  return 0;
}
