// arb_zeta.cpp
// windowed zeta calculator
// Looking for extremes of S(t)
//
// based on win_zeta1.10.c
//
// Created: September 2024
// Copyright DJ Platt 2024
//
// DJ Platt
// University of Bristol
//

#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "time.h"
#include "mpfr.h"
#include "flint/arb.h"
#include "flint/acb.h"
#include "arb_win_zeta.h"
#include "arb_fft.h"

int HI_PREC;
#include "parameters.h"
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
arb_t *exps1,intererr;
arf_t ip_tmp1,ip_tmp2,ip_tmp3,ip_tmp4;
arf_t two_101,offset;
arb_t msinc_tmp,arb_two_pi_B,mip_tmp,minter_tmp,minter_tmp1;
//arf_t znew_guess,zold_guess,ztn;
arb_t misin;

// note (sign&&sign)==0 <=> signs different and known 

void print_usage(const char* cmd)
{
  printf("Usage:- %s prec t_min num_its step?.\n",cmd);
  printf("Verifies RH in [t_min+step*n,t_min+step*(n+1)] for n=0..num_its-1.\n");
  exit(1);
}

// inits and sets z to [-d,d]+[-d,d]i
void set_err(acb_t z, double d)
{
  acb_init(z);
  arb_set_d(acb_realref(z),d);
  arb_add_error(acb_imagref(z),acb_realref(z));
  arb_set(acb_realref(z),acb_imagref(z));
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
  arb_set_d(tmp,d);
  arb_init(x);
  arb_add_error(x,tmp);
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
  //printf("s_vec[0]=");arb_printd(s_vec[0],10);printf("\n");
  for(i=0;i<M;i++)
    {
      buck[i]=arb_get_d(s_vec[i])*B+0.5;
      arb_set_ui(tmp,buck[i]);
      arb_div(tmp,tmp,arb_B,prec);
      //printf("n=%d goes in bucket %d\n",i+1,buck[i]);
      arb_sub(s_vec[i],s_vec[i],tmp,prec);
    }
  //printf("Before shuffle, n=1 term is in bucket %d\n",buck[0]);
  if(buck[0]<=0) // keep conj_order simple
    fatal_error("buck[0]<=0.");
  printf("Final bucket used is at %d.\n",buck[M-1]);
  if(buck[M-1]>=N1/2) // ditto
    fatal_error("buck[M-1]>=N/2.");
  for(i=0;i<M;i++)
    buck[i]=conj_order(buck[i]);
  //printf("buck[0]=%d\nbuck[M-1]=%d\n",buck[0],buck[M-1]);
  //printf("n=1 term is in bucket %d\n",buck[0]);
  arb_clear(tmp);
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
  //printf("gamma_factor(%f) returning ",t);acb_printd(res,20);printf("\n");
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
      //printf("s_vec[%ld] = ",i);arb_printd(s_vec[i],20);printf("\n");
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

  // the Taylor error from the Pari script already includes
  // 2sqrt(J)-1
  /*
  arb_set_ui(A,M);
  arb_sqrt(A,A,prec);
  arb_mul_2exp_si(A,A,1);
  arb_sub_ui(A,A,1,prec); // 2sqrt(M)-1
  acb_mul_arb(tayerr,tayerr,A,prec);
  */
  printf("Total Taylor Error set to ");acb_printd(tayerr,30);printf("\n");
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
    acb_add(G_vec[i],g_vec[i],gtwiderr,prec); // G=g~
  //printf("k=%lu, last g~ = ",k);acb_printd(G_vec[N1-1],20);printf("\n");

  fft(G_vec,N1,ws1_f,prec);

  for(i=0;i<=N1/2;i++)
    {
      acb_div_arb(G_vec[i],G_vec[i],A1,prec);
      acb_add(G_vec[i],G_vec[i],Gtwiderr,prec);
      if(i&1)
	acb_neg(G_vec[i],G_vec[i]);
    }
  //printf("k=%lu, N1/2'th G = ",k);acb_printd(G_vec[N1/2],20);printf("\n");

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
    {
      my_convolve(f_vec,skn_vec,G_vec,N1,ws1_r,ws1_f,prec);
      //printf("k=0 0'th convolution = ");acb_printd(f_vec[0],20);printf("\n");
    }
  else // convolve into G_vec, then add into f_vec 
    {
      my_convolve(G_vec,skn_vec,G_vec,N1,ws1_r,ws1_f,prec);
      //printf("k=%lu 0'th convolution = ",k);acb_printd(G_vec[0],20);printf("\n");
      for(i=0;i<=N1/2;i++)
	acb_add(f_vec[i],f_vec[i],G_vec[i],prec);
    }
}

// compute exp(-t^2/2H^2)
void inter_gaussian(arb_ptr res, arb_ptr t, int64_t prec)
{
  //arf_printf("inter_gaussian called with t=%.40Re\n",t);
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
  // compute f_res
  arb_mul(ip_tmp2,ip_tmp1,ip_tmp2,prec); //f*exp
  arb_mul(f_res,ip_tmp2,ip_tmp3,prec); //f*exp*sinc
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
  arb_neg(fd_res,fd_res);
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
  //printf("In inter_t\n");
  double t0_d=arb_get_d(t_ptr),dt0;
  arb_set_si(inter_tmp,t0_d);
  if(arb_equal(t_ptr,inter_tmp)) // t_ptr is integral so = a lattice pt
    nearest_t0=t0_d+1.5;
  else
    nearest_t0=t0_d+0.5;
  //printf("Nearest t0=%ld\n",nearest_t0);
  arb_zero(f_res);
  if (fd_res) arb_zero(fd_res);
  sincp=true;
  if(nearest_t0>t0_d) // start from right of t0
    {
      for(i=nearest_t0,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
      sincp=true;
      for(i=nearest_t0-INTER_SPACING,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
    }
  else // start from left
    {
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
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
    }
}

void inter_t(arb_ptr f_res, acb_t *f_vec, arb_ptr t_ptr, int64_t prec)
{
  inter_t(f_res, NULL, f_vec, t_ptr, prec);
}

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


// we have found a stationary pt between left, left+1 and left+2
bool resolve_stat_point(arb_ptr tl, arb_ptr ftl, arb_ptr tm, arb_ptr ftm, arb_ptr tr, arb_ptr ftr, int left, sign_t this_sign, acb_t *f_vec, int64_t prec)
{
  static bool init=false;
  static arb_t sp_t,sp_ft,sp_t1,sp_ft1;  
  if(!init)
    {
      init=true;
      arb_init(sp_t);
      arb_init(sp_ft);
      arb_init(sp_t1);
      arb_init(sp_ft1);
    }
  dir_t dir1,dir2,dir3,dir4;
  //printf("In resolve_stat_points at %d %d %d.\n",left,left+1,left+2);
  //arf_t t,ft;
  //arf_init(t);arf_init(ft);
  arb_set_ui(tl,left);
  arb_set_ui(tm,left+1);
  arb_set_ui(tr,left+2);
  arb_set(ftl,acb_realref(f_vec[left]));
  arb_set(ftm,acb_realref(f_vec[left+1]));
  arb_set(ftr,acb_realref(f_vec[left+2]));

  while(true)
    {
      //arf_printf("Looping with t:\nleft=%.40Re\nmid=%.40Re\nright=%.40Re\n",tl,tm,tr);
      //arf_printf("F(t):\nleft=%.40Re\nmid=%.40Re\nright=%.40Re\n",ftl,ftm,ftr);

      arb_add(sp_t,tl,tm,prec);
      arb_mul_2exp_si(sp_t,sp_t,-1);
      inter_t(sp_ft,NULL,f_vec,sp_t,prec);
      if((sign(sp_ft)&this_sign)==0)
	{
	  //arf_printf("Sign change found at t=%.40Re\n",sp_t);
	  arb_set(tr,tm);
	  arb_set(ftr,ftm);
	  arb_set(tm,sp_t);
	  arb_set(ftm,sp_ft);
	  return(true);
	}
      dir1=dir(ftl,sp_ft);
      dir2=dir(sp_ft,ftm);
      if(((this_sign==POS)&&(dir1==DOWN)&&(dir2==UP))||((this_sign==NEG)&&(dir1==UP)&&(dir2==DOWN)))
	{
	  arb_set(tr,tm);
	  arb_set(ftr,ftm);
	  arb_set(tm,sp_t);
	  arb_set(ftm,sp_ft);
	  continue;
	}

      arb_add(sp_t1,tm,tr,prec);
      arb_mul_2exp_si(sp_t1,sp_t1,-1);
      inter_t(sp_ft1,NULL,f_vec,sp_t1,prec);
      if((sign(sp_ft1)&this_sign)==0)
	{
	  //arf_printf("Sign change found at t=%.40Re\n",sp_t);
	  arb_set(tl,tm);
	  arb_set(ftl,ftm);
	  arb_set(tm,sp_t1);
	  arb_set(ftm,sp_ft1);
	  return(true);
	}
      dir3=dir(ftm,sp_ft1);
      dir4=dir(sp_ft1,ftr);
      if(((this_sign==POS)&&(dir2==DOWN)&&(dir3==UP))||((this_sign==NEG)&&(dir2==UP)&&(dir3==DOWN)))
	{
	  arb_set(tl,sp_t);
	  arb_set(ftl,sp_ft);
	  arb_set(tr,sp_t1);
	  arb_set(ftr,sp_ft1);
	  continue;
	}
      if(((this_sign==POS)&&(dir3==DOWN)&&(dir4==UP))||((this_sign==NEG)&&(dir3==UP)&&(dir4==DOWN)))
	{
	  arb_set(tl,tm);
	  arb_set(ftl,ftm);
	  arb_set(tm,sp_t1);
	  arb_set(ftm,sp_ft1);
	  continue;
	}
      printf("Stat point Process failed to converge.\n");
      return(false);
    }
}

arb_t S_of_t_plus,S_of_t_minus;
void St(arb_t t, int64_t NNt,int64_t prec)
{
  static bool init=false;
  static arb_t t2,tmp1,tmp2,Nt;
  static acb_t z,lng;
  if(!init)
    {
      init=true;
      arb_init(S_of_t_plus);
      arb_init(S_of_t_minus);
      arb_init(t2);
      arb_init(tmp1);
      arb_init(tmp2);
      acb_init(z);
      acb_init(lng);
      arb_init(Nt);
      arb_set_d(acb_realref(z),0.25);
    }
  arb_mul_2exp_si(t2,t,-1);
  arb_set(acb_imagref(z),t2);
  acb_lgamma(lng,z,prec);
  arb_mul(tmp1,t2,arb_ln_pi,prec);
  arb_sub(tmp2,acb_imagref(lng),tmp1,prec);
  arb_div(tmp1,tmp2,arb_pi,prec);
  arb_set_si(Nt,NNt);
  arb_sub(tmp2,Nt,tmp1,prec);
  printf("S+(t) = ");arb_printd(tmp2,20);printf("\n");
  arb_sub_ui(tmp1,tmp2,1,prec);
  printf("S-(t) = ");arb_printd(tmp1,20);printf("\n");
  
}

// given left and right such that start+[left,right]
// destroys left
void arb_st(arb_t left, arb_t right, double dstart, long int zero, int64_t prec)
{
  static bool init=false;
  static arb_t arb_one_over_A,astart;
  if(!init)
    {
      init=true;
      arb_init(arb_one_over_A);
      arb_set_d(arb_one_over_A,one_over_A);
      arb_init(astart);
    }
  //printf("In arb_st with:\n   left = ");
  //arb_printd(left,30);printf("\n  right = ");
  //arb_printd(right,30);printf("\n  start = ");
  //arb_printd(start,30);
  //printf("\n    1/A = %f\n",one_over_A);
  
  arb_set_d(astart,dstart);
  arb_union(left,left,right,prec);
  arb_mul(left,left,arb_one_over_A,prec);
  arb_add(left,left,astart,prec);
  // left now brackets zero number <zero>
  St(left,zero,prec);
  printf("Zero number %ld at ",zero);arb_printd(left,30);printf("\n");
  //exit(0);
}


void int_st(int left, int right, double dstart, long int zero, int64_t prec)
{
  static bool init=false;
  static arb_t aleft,aright;
  if(!init)
    {
      init=true;
      arb_init(aleft);
      arb_init(aright);
    }
  arb_set_si(aleft,left);
  arb_set_si(aright,right);
  arb_st(aleft,aright,dstart,zero,prec);
}


int zeros_st(acb_t *f_vec, int start, double d_start, int end, long int zeros_to_date, int64_t prec)
{
  static bool init=false;
  static arb_t tl,fl,tm,fm,tr,fr,arb_start;
  if(!init)
    {
      init=true;
      arb_init(tl);
      arb_init(fl);
      arb_init(tm);
      arb_init(fm);
      arb_init(tr);
      arb_init(fr);
      arb_init(arb_start);
    }
  printf("In zeros_st with dstart=%f\n",d_start);
  
  long int i=start+1,count=0;
  dir_t last_dir,this_dir;
  sign_t last_sign,this_sign,res1_sign,res2_sign;
  last_sign=sign(f_vec[start]);
  last_dir=UNK;
  if(last_sign==UNK)
    last_sign=sign(f_vec[i++]);
  if(last_sign==UNK)
    {
      printf("Two UNKs in zeros_st.\n");
      return(count);
    }
  for(;i<=end;i++)
    {
      this_sign=sign(f_vec[i]);
      this_dir=dir(f_vec[i-1],f_vec[i]);
      if((this_sign&last_sign)==0) // valid sign change
	{
	  int_st(i-start-1,i-start,d_start,zeros_to_date+count+1,prec);
	  //printf("Zero %ld found at [%f,%f]\n",zeros_to_date+count+1,d_start+(i-1-start)*one_over_A,d_start+(i-start)*one_over_A);
	  count++;
	  last_sign=this_sign;
	  last_dir=this_dir;
	  continue;
	}
      if((this_dir&last_dir)==0) // maximum or minimum
	{
	  if(((this_dir==UP)&&(this_sign==POS))||((this_dir==DOWN)&&(this_sign==NEG)))
	    {
	      //printf("Stat point found at %ld.\n",i);
	      if(!resolve_stat_point(tl,fl,tm,fm,tr,fr,i-2,this_sign,f_vec,prec))
		return(count);
	      arb_sub_si(tl,tl,start,prec);
	      arb_sub_si(tm,tm,start,prec);
	      arb_sub_si(tr,tr,start,prec);
	      printf("Stat Point:\n");
	      arb_st(tl,tm,d_start,zeros_to_date+count+1,prec);
	      arb_st(tm,tr,d_start,zeros_to_date+count+2,prec);
	      //printf("looking between i=[%ld,%ld]\n",i-2,i);
	      //printf("Zero %ld found at ",zeros_to_date+count+1);
	      //arb_union(tl,tl,tm,prec);
	      //arb_printd(tl,10);
	      //printf("\nZero %ld found at ",zeros_to_date+count+2);
	      //arb_union(tm,tm,tr,prec);
	      //arb_printd(tm,10);
	      //printf("\n");
	      //printf("Zeros %ld and %ld found at [%f,%f]\n",zeros_to_date+count+1,zeros_to_date+count+2,d_start+(i-start-2)*one_over_A,d_start+(i-start)*one_over_A);
	      last_dir=this_dir;
	      count++;
	      count++;
	    }
	  last_dir=this_dir;
	}
    }
  return(count);
}

// int Im loggamma(1/4+it/2)
void im_int1(arb_t res, arb_t t, int64_t prec)
{
  arb_mul_2exp_si(im_t,t,2); // 4t
  arb_atan(res,im_t,prec);
  //printf("atan(4t)=");mpfi_printn(res,30);
  arb_mul(res,res,t,prec);
  arb_mul_d(res,res,-0.25,prec); // -t/4*atan(4t)
  //printf("-t/4atan(4t)=");mpfi_printn(res,30);
  arb_mul(im_t,t,t,prec); // t^2
  arb_mul_d(im_1,im_t,-3.0/4.0,prec);
  arb_add(res,res,im_1,prec); //-3/4*t^2
  arb_mul_ui(im_1,im_t,16,prec); // 16t^2
  arb_add_ui(im_1,im_1,1,prec);
  arb_log(im_1,im_1,prec); // log(1+16t^2)
  arb_mul_d(im_1,im_1,1.0/32.0,prec);
  arb_add(res,res,im_1,prec);
  arb_add_d(res,res,-1.0/64.0,prec);
  arb_add_d(im_t,im_t,1.0/16.0,prec); // t^2+1/16
  arb_log(im_1,im_t,prec);
  arb_mul(im_t,im_t,im_1,prec); // (t^2+1/16)log(t^2+1/16)
  arb_mul_d(im_t,im_t,0.25,prec);
  arb_add(res,res,im_t,prec);
}
  
void im_int(arb_t res, arb_t t0, arb_t t1, int64_t prec)
{
  arb_sub(im_err,t1,t0,prec);
  arb_div_ui(im_err,im_err,4,prec);
  arb_div(im_err,im_err,t0,prec);
  //mpfi_neg(im_2,im_err);
  //mpfi_put(im_err,im_2);
  arb_div_ui(im_t0,t0,2,prec);
  arb_div_ui(im_t1,t1,2,prec);
  im_int1(res,im_t1,prec);
  im_int1(im_2,im_t0,prec);
  arb_sub(res,res,im_2,prec);
  arb_mul_ui(res,res,2,prec);
  arb_add_error(res,im_err);
}

bool stat_pt(acb_ptr l, acb_ptr m, acb_ptr r)
{
  sign_t sm=sign(m);
  if(sm==UNK)
    return false;
  if(sm==POS)
    return arb_gt(acb_realref(l),acb_realref(m))&&arb_gt(acb_realref(r),acb_realref(m));
  else
    return arb_gt(acb_realref(m),acb_realref(l))&&arb_gt(acb_realref(m),acb_realref(r));
}

double Nleft_int(long int t0_ptr, long int t1_ptr, acb_t *f_vec, double delta, int64_t prec)
{
  static bool init=false;
  static arb_t t1,t2,t3,t4,t5,t6;
  if(!init)
    {
      init=true;
      arb_init(t1);
      arb_init(t2);
      arb_init(t3);
      arb_init(t4);
      arb_init(t5);
      arb_init(t6);
    }
  long int res=0;
  long int ptr=t0_ptr,last_ptr;
  sign_t last_sign=sign(f_vec[ptr++]),this_sign;
  //printf("In Nleft_int with %ld %ld\n",t0_ptr,t1_ptr);
  while(last_sign==UNK)
    last_sign=sign(f_vec[ptr++]);
  last_ptr=ptr-1;
  for(;ptr<=t1_ptr;)
    {
      this_sign=sign(f_vec[ptr++]);
      if(this_sign==last_sign) // no sign change here, move on
	{
	  if(stat_pt(f_vec[ptr-3],f_vec[ptr-2],f_vec[ptr-1]))
	    {
	      arb_printd(acb_realref(f_vec[ptr-3]),10);printf(" ");
	      arb_printd(acb_realref(f_vec[ptr-2]),10);printf(" ");
	      arb_printd(acb_realref(f_vec[ptr-1]),10);printf("\n");
	      
	      if(resolve_stat_point(t1,t2,t3,t4,t5,t6,ptr-3,sign(acb_realref(f_vec[ptr-2])),f_vec,prec))
		{printf("Resolved.\n");//res-=2*(ptr-3-t0_ptr+1);
		}
	    }
	  last_ptr=ptr-1;
	  continue;
	}
      if(this_sign==UNK) // might be a sign change coming
	continue;
      // definately a sign change
      //printf("Sign change in Nleft at %ld %ld\n",last_ptr,ptr-1);
      res-=(last_ptr-t0_ptr+1);
      //printf("res=%10.8e\n",res*delta);
      last_ptr=ptr-1;
      last_sign=this_sign;
    }
  //printf("Nleft_int returning %f\n",res*delta);
  return(res*delta);
}

double Nright_int(long int t0_ptr, long int t1_ptr, acb_t *f_vec, double delta,int64_t prec)
{
  static bool init=false;
  static arb_t t1,t2,t3,t4,t5,t6;
  if(!init)
    {
      init=true;
      arb_init(t1);
      arb_init(t2);
      arb_init(t3);
      arb_init(t4);
      arb_init(t5);
      arb_init(t6);
    }
  //printf("In Nright_int with %ld %ld\n",t0_ptr,t1_ptr);
  //for(uint64_t i=t0_ptr;i<=t1_ptr;i++)
  //{printf("%lu ",i);arb_printd(acb_realref(f_vec[i]),20);printf("\n");}
  long int res=0;
  long int ptr=t0_ptr;
  sign_t last_sign=sign(f_vec[ptr++]),this_sign;
  while(last_sign==UNK)
    last_sign=sign(f_vec[ptr++]);
  for(;ptr<=t1_ptr;)
    {
      this_sign=sign(f_vec[ptr++]);
      if(this_sign==last_sign) // no sign change here, check for stat pt
	{
	  if(stat_pt(f_vec[ptr-3],f_vec[ptr-2],f_vec[ptr-1]))
	    {
	      arb_printd(acb_realref(f_vec[ptr-3]),10);printf(" ");
	      arb_printd(acb_realref(f_vec[ptr-2]),10);printf(" ");
	      arb_printd(acb_realref(f_vec[ptr-1]),10);printf("\n");
 
	      if(resolve_stat_point(t1,t2,t3,t4,t5,t6,ptr-3,sign(acb_realref(f_vec[ptr-2])),f_vec,prec))
		{printf("Resolved.\n");res+=2*(t1_ptr-ptr);
		  //printf("res=%10.8e\n",res*delta);
		}
	    }
	}
      else
	{
	  if(this_sign!=UNK) // definately a sign change
	    {
	      //printf("Sign change in Nright at %ld %ld\n",ptr-2,ptr-1);
	      res+=(t1_ptr-ptr+1);
	      //printf("res=%10.8e\n",res*delta);
	      last_sign=this_sign;
	    }
	}
    }
  
  //printf("Nright_int returning %f\n",res*delta);
  return(res*delta);
}

// returns int_t0^{t0+h} S(t) dt
// t0=t0_ptr*delta
// t0>168*Pi
// Trudgian Thesis Theorem 5.2.2
void St_int(arb_t res, arb_t t, int64_t prec) 
{
  arb_log(res,t,prec);
  arb_mul_ui(res,res,59,prec);
  arb_div_ui(res,res,1000,prec);
  arb_set_ui(tm_tmp,2067);
  arb_div_ui(tm_tmp,tm_tmp,1000,prec);
  arb_add(res,res,tm_tmp,prec);
}

// the log term integrated
void ln_term(arb_t res, arb_t t0, arb_t t1, arb_t h1, int64_t prec)
{
  arb_add(res,t0,t1,prec);
  arb_mul_d(res,res,-0.25,prec);
  arb_mul(res,res,arb_ln_pi,prec);
  arb_mul(res,res,h1,prec);
}

// the maximum number of zeros <=a based on Turing region [a,b]
long int turing_max(acb_t *f_vec, long int a_ptr, double a, long int b_ptr, double b, int64_t prec)
{
  static bool init=false;
  static fmpz_t fz;
  if(!init)
    {
      init=true;
      fmpz_init(fz);
    }
  //printf("In turing max\n");
  arb_set_d(tm_t0,a);
  //printf("t0=");mpfi_printn(tm_t0,10);
  arb_set_d(tm_t1,b);
  //printf("t1=");mpfi_printn(tm_t1,10);  
  arb_sub(tm_h,tm_t1,tm_t0,prec);
  //  
  St_int(tm_res,tm_t1,prec);
  //printf("int S(t)=");arb_printd(tm_res,10);printf("\n");
  arb_sub_d(tm_res,tm_res,Nright_int(a_ptr,b_ptr,f_vec,one_over_A,prec),prec);
  //printf("\nint S(t) - int Nright=");arb_printd(tm_res,10);
  ln_term(tm_tmp,tm_t0,tm_t1,tm_h,prec);
  //printf("\nint - t log(pi)/2=");arb_printd(tm_tmp,10);
  im_int(tm_tmp1,tm_t0,tm_t1,prec);
  //printf("\nIm lnGamma term=");arb_printd(tm_tmp1,10);
  arb_add(tm_tmp,tm_tmp,tm_tmp1,prec);
  arb_div(tm_tmp,tm_tmp,arb_pi,prec);
  arb_add(tm_res,tm_res,tm_tmp,prec);
  arb_div(tm_tmp,tm_res,tm_h,prec);
  //printf("\nturing_max computed ");arb_printd(tm_tmp,20);printf("\n");
  arb_floor(tm_res,tm_tmp,prec);
  if(!arb_get_unique_fmpz(fz,tm_res))
    {
      printf("Turing max did not bracket an integer. Exiting.\n");
      exit(0);
    }
  return(fmpz_get_si(fz)+1); // +1 because of sqrt(omega) =+/-1
}

long int turing_min(acb_t *f_vec, long int a_ptr, double a, long int b_ptr, double b, int64_t prec)
{
  static bool init=false;
  static fmpz_t fz;
  if(!init)
    {
      init=true;
      fmpz_init(fz);
    }
  //printf("In Turing Min\n");
  arb_set_d(tm_t0,a);
  //printf("t0=");arb_printd(tm_t0,10);  
  arb_set_d(tm_t1,b);;
  //printf("\nt1=");arb_printd(tm_t1,10); 
  arb_sub(tm_h,tm_t1,tm_t0,prec);
  //  
  St_int(tm_res,tm_t1,prec);
  arb_neg(tm_res,tm_res);
  //printf("\nint S(t)=");arb_printd(tm_res,10);
  arb_sub_d(tm_res,tm_res,Nleft_int(a_ptr,b_ptr,f_vec,one_over_A,prec),prec);
  //printf("\nint S(t) - int Nright");arb_printd(tm_res,10);
  ln_term(tm_tmp,tm_t0,tm_t1,tm_h,prec);
  //printf("\nint - t log(pi)/2=");arb_printd(tm_tmp,10);
  im_int(tm_tmp1,tm_t0,tm_t1,prec);
  //printf("\nIm lnGamma term=");arb_printd(tm_tmp1,10);
  arb_add(tm_tmp,tm_tmp,tm_tmp1,prec);
  arb_div(tm_tmp,tm_tmp,arb_pi,prec);
  arb_add(tm_res,tm_res,tm_tmp,prec);
  arb_div(tm_res,tm_res,tm_h,prec);
  //printf("\nturing_min computed ");arb_printd(tm_res,20);printf("\n");
  arb_ceil(tm_res,tm_res,prec);
  if(!arb_get_unique_fmpz(fz,tm_res))
    {
      printf("Turing min did not bracket an integer. Exiting.\n");
      exit(0);
    }
  return(fmpz_get_si(fz)+1); // +1 because of sqrt(omega) =+/-1 
}

// Use Turing's method to estimate and then check number of zeros
long int turing(acb_t *f_vec, long int a_ptr, double a, long int b_ptr, double b, long int last_max, int64_t prec)
{
  // etimate maximum for N(b)

  long int i;

  long int min_lo,max_hi=turing_max(f_vec,b_ptr,b,b_ptr+TURING_LEN,b+TURING_WIDTH,prec);
  if(last_max==0) // No previous run or previous run failed, so estimate min N(a)
    min_lo=turing_min(f_vec,a_ptr-TURING_LEN,a-TURING_WIDTH,a_ptr,a,prec);
  else // If previous run succeeded, zeros to height N(a) is known
    min_lo=last_max;
  printf("looking for %ld-%ld=%ld zeros\n",max_hi,min_lo,max_hi-min_lo);
  long int num_found,num_exp=max_hi-min_lo;
  sign_t ls,rs;
  /*
  for(uint64_t i=0;i<100;i++)
    {
      acb_printd(f_vec[a_ptr+i],20);
      printf("\n");
    }
  */
  ls=sign(f_vec[a_ptr]);
  if(ls==UNK)
    {
      ls=sign(f_vec[a_ptr+1]);
      if(ls==UNK)
	{
	  printf("Missed All/All zeros in region %f to %f.\n",a,b);
	  return(0);
	}
    }
  rs=sign(f_vec[b_ptr]);
  if(rs==UNK)
    {
      rs=sign(f_vec[b_ptr+1]);
      if(rs==UNK)
	{
	  printf("Missed All/All zeros in region %f to %f.\n",a,b);
	  return(0);
	}
    }

  //for(i=a_ptr,ls=sign(f_vec[i]);ls==UNK;i++);
  //for(i=b_ptr,rs=sign(f_vec[i]);rs==UNK;i++);
  i=0;

  if(ls==rs) // same sign at both ends, count should be even
    {
      if(num_exp&1)
	{
	  printf("Problem in Turing Zone.\n");
	  printf("Missed All/All zeros in region %f to %f.\n",a,b);
	  return(0);
	}
    }
  else // count should be odd
    {
      if(!(num_exp&1)) 
	{
	  printf("Problem in Turing Zone.\n");
	  printf("Missed All/All zeros in region %f to %f.\n",a,b);
	  return(0);
	}
    }

  num_found=zeros_st(f_vec,a_ptr,a,b_ptr,min_lo,prec); // go find zeros, using stat pts as well

  if(num_exp==num_found)
    {
	printf("All %ld zeros found in region %f to %f using stat points.\n",num_exp,a,b);
	return(max_hi); // use this for next iteration
    }
  printf("Missed %ld/%ld zeros (running with stat points) in region %f to %f.\n",num_exp-num_found,num_exp,a,b);
  return(0);
}

int main(int argc, char **argv)
{
  long int prec,num_its,it;
  long int i,j,k,int_step,last_max=0;
  double t0,t1,step,t;
  fpos_t fpos;
  arb_t tmp;
  bool first=true;
  acb_t omega;
  time_t last_time=time(NULL),this_time;
  if(argc!=5)
    print_usage(argv[0]);

  prec=atoi(argv[1]);


  t0=atof(argv[2]);
  num_its=atoi(argv[3]);
  if(num_its<=0)
    print_usage(argv[0]);
  step=atof(argv[4]);
  if(step<=0)
    print_usage(argv[0]);
  if(t0<T0_MIN)
    fatal_error("t0 outside limits.");
  if((t0+step*num_its)>T0_MAX)
    fatal_error("t0 outside limits.");

  t0=t0+step/2.0;

  HI_PREC=prec+EXTRA_PREC;
  int_step=step/one_over_A;
  if((int_step&1)!=0) // must be even
    fatal_error("int_step must be even.");
  if(step/one_over_A-(double) int_step !=0.0) // must be an exact number of steps
    fatal_error("step must be integral.");

  exps1=(arb_t *) malloc((TURING_LEN+int_step/2+Ns*INTER_SPACING)*sizeof(arb_t));
  if(!exps1)
    {
      printf("Failed to allocate memory for exps1. Exiting.\n");
      exit(0);
    }

  for(i=0,t=one_over_A;i<TURING_LEN+int_step/2+Ns*INTER_SPACING;t+=one_over_A,i++)
    {
      arb_init(exps1[i]);
      arb_set_d(exps1[i],t);
      arb_div_d(exps1[i],exps1[i],h,prec);
      arb_mul(exps1[i],exps1[i],exps1[i],prec);
      arb_mul_2exp_si(exps1[i],exps1[i],-1);
      arb_exp(exps1[i],exps1[i],prec);
    }

#ifdef UPSAM_HIGH
  printf("Using high precision upsampling parameters. error=%e\n",intererr_d);
#endif
  printf("Aiming to do %ld iteration(s) starting at t=%13.11e.\n",num_its,t0-step/2.0);

  for(it=0;it<num_its;it++,t0+=step)
    {
      //printf("Running centred at t0=%f.\n",t0);
      if(first)
	{
	  first=false;
	  //printf("Calling init.\n");
	  //system("date");
	  init(t0,prec);
	  arb_init(tmp);
	  acb_init(omega);
	  arb_div_ui(tmp,arb_2_pi,N,prec);
	  arb_sin_cos(acb_imagref(omega),acb_realref(omega),tmp,prec);
	  arb_clear(tmp);
	  //acb_printd(omega,20);printf("\n");exit(0);
	  //printf("Init finished.\n");
	  //system("date");
	}
      else
	{
	  //printf("Calling re-init.\n");
	  //system("date");
	  re_init(t0,prec);
	  //printf("re-init finished.\n");
	  //system("date");
	}

      this_time=time(NULL);
      printf("Time to (re-)initialise = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;

      for(k=0;k<K;k++)
	{
	  G_k(k,prec);
	  do_conv(k,prec);

	}
      this_time=time(NULL);
      printf("Time to convolve = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;



      //printf("convolutions finished\n");
      //system("date");

      // upsample by "zero" padding
      for(i=N1/2;i<=N/2;i++)
	acb_set(f_vec[i],Fmaxerr);

      for(i=0;i<=N/2;i++)
	{
	  acb_add(f_vec[i],f_vec[i],fhatsumerr,prec);
	  acb_add(f_vec[i],f_vec[i],tayerr,prec);
	  acb_add(f_vec[i],f_vec[i],fhattwiderr,prec);
	}

      // f is real so use conjugate symmetry
      // could use a real idtf instead
      for(i=N/2+1;i<N;i++)
	acb_conj(f_vec[i],f_vec[N-i]);

      for(i=1;i<N;i+=2)
	acb_neg(f_vec[i],f_vec[i]);

      hermidft(f_vec,N/2,ws_r,omega,prec);

      /*
      printf("Post iFFT:-\n");
      for(i=0;i<50;i++)
	{
	  acb_printd(f_vec[i+N/2],20);
	  printf("\n");
	}
      */

      this_time=time(NULL);
      printf("Time for final iFFT = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;


      //printf("Final iFFT finished.\n");
      //system("date");

      // remove the Gaussian from the region of interest
      for(j=0,i=N/2+1;i<=N/2+TURING_LEN+int_step/2+INTER_SPACING*Ns;j++,i++)
	{
	  //printf(" before f=");acb_printd(f_vec[i],20);printf("\n");
	  arb_div_d(acb_realref(f_vec[i]),acb_realref(f_vec[i]),B,prec);
	  //printf("f=");acb_printd(f_vec[i],20);printf("\n");
	  arb_add(acb_realref(f_vec[i]),acb_realref(f_vec[i]),acb_realref(ftwiderr),prec);
	  //printf("f=");acb_printd(f_vec[i],20);printf("\n");
	  arb_mul(acb_realref(f_vec[i]),acb_realref(f_vec[i]),exps1[j],prec);
	  //printf("after f(%8.6e)=",(double)j*one_over_A);acb_printd(f_vec[i],20);printf("\n");
	}
      for(j=0,i=N/2-1;i>=N/2-TURING_LEN-int_step/2-INTER_SPACING*Ns;j++,i--)
	{
	  arb_div_d(acb_realref(f_vec[i]),acb_realref(f_vec[i]),B,prec);
	  arb_add(acb_realref(f_vec[i]),acb_realref(f_vec[i]),acb_realref(ftwiderr),prec);
	  arb_mul(acb_realref(f_vec[i]),acb_realref(f_vec[i]),exps1[j],prec);
	}
      
      arb_div_d(acb_realref(f_vec[N/2]),acb_realref(f_vec[N/2]),B,prec);
      arb_add(acb_realref(f_vec[N/2]),acb_realref(f_vec[N/2]),acb_realref(ftwiderr),prec);

      /*      
      printf("After normalisation.\n");
      
      for(i=0;i<10;i++)
	{
	  printf("t=%20.18e ",t0+(double)i*one_over_A);
	  arb_printd(acb_realref(f_vec[N/2+i]),20);
	  printf("\n");
	}

      for(i=0;i<10;i++)
	{
	  printf("t=%20.18e ",t0+(double)(i+int_step/2)*one_over_A);
	  arb_printd(acb_realref(f_vec[N/2+i+int_step/2]),20);
	  printf("\n");
	}
      */

      

      last_max=turing(f_vec,N/2-int_step/2,t0-int_step/2*one_over_A,N/2+int_step/2,t0+int_step/2*one_over_A,last_max,prec);
      this_time=time(NULL);
      printf("Time to isolate zeros = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;

    }
}
