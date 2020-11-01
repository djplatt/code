// arb_zeta.c
// windowed zeta calculator
//
// based on win_zeta1.10.c
//
// Created: 27th January 2017
//
// DJ Platt
// University of Bristol
//

#include "stdlib.h"
#include "stdio.h"
#include "stdint.h"
#include "stdbool.h"
#include "time.h"
#include "mpfr.h"
#include "arb.h"
#include "acb.h"

#include "parameters.h"
#include "arb_win_zeta.h"
#include "arb_fft.h"
#include "turing.h"

int HI_PREC;
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

arb_t A,A1,two_pi_t[N1],s_vec[M],sk_vec[M],exps[N1],sqrts[M],logisqrtpi[M];
acb_t g_vec[N1],G_vec[N1],f_vec[N],skn_vec[N1],ws_r[N/4],ws1_r[N1/2],ws1_f[N1/2],n_vec[M];
acb_t gtwiderr,Gtwiderr,c_tmp,fhatsumerr,fhattwiderr,ftwiderr,tayerr,Fmaxerr;
int buck[M]; // which bucket does log(nsqrt(Pi)/2pi) go in.
arb_t *exps1,intererr;
arf_t ip_tmp1,ip_tmp2,ip_tmp3,ip_tmp4;
arf_t offset;
arb_t msinc_tmp,mip_tmp,minter_tmp,minter_tmp1;
//arf_t znew_guess,zold_guess,ztn;
arb_t misin;

// note (sign&&sign)==0 <=> signs different and known 


// inits and sets z to [-d,d]+[-d,d]i
void acb_set_err(acb_ptr z, double d)
{
  acb_init(z);
  arb_set_d(acb_realref(z),d);
  arb_add_error(acb_imagref(z),acb_realref(z));
  arb_set(acb_realref(z),acb_imagref(z));
}

// inits and sets x to [-d,d]
void arb_set_err(arb_t x, double d)
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
  arb_t tmp,arb_B;
  int i;
  double x;
  arb_init(tmp);
  arb_init(arb_B);
  arb_set_d(arb_B,B);
  for(i=0;i<M;i++)
    {
      buck[i]=arb_get_d(s_vec[i])*B+0.5;
      arb_set_ui(tmp,buck[i]);
      arb_div(tmp,tmp,arb_B,prec);
      arb_sub(s_vec[i],s_vec[i],tmp,prec);
    }
  if(buck[0]<=0) // keep conj_order simple
    fatal_error("buck[0]<=0.");
  printf("Final bucket used is at %d.\n",buck[M-1]);
  if(buck[M-1]>=N1/2) // ditto
    fatal_error("buck[M-1]>=N/2.");
  for(i=0;i<M;i++)
    buck[i]=conj_order(buck[i]);
  arb_clear(arb_B);
  arb_clear(tmp);
}

arb_t atan1_tmp1,atan1_tmp2;//,im_1,im_2,im_t,im_err;
//arb_t tm_h,tm_res,tm_tmp,tm_tmp1,im_t0,im_t1;
//arf_t tm_arf,sp_t,sp_ft,sp_t1,sp_ft1;
//arb_t tm_t0,tm_t1,init_tmp1,sqrt_pi,pi_hi,half_log_2_pi;

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
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,-1); // pi/2
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
  arb_t init_tmp1,sqrt_pi,half_log_2_pi,pi_hi,arb_2_pi;
  arb_init(init_tmp1);
  arb_init(sqrt_pi);
  arb_init(half_log_2_pi);
  arb_init(pi_hi);
  arb_init(arb_2_pi);
  arb_const_pi(arb_2_pi,prec);
  arb_mul_2exp_si(arb_2_pi,arb_2_pi,1);
  int i;
  double t;
  arf_t arf_tmp;
  arf_init(arf_tmp);
  //arb_init(arb_pi);
  //arb_const_pi(arb_pi,prec);
  //arb_init(arb_ln_pi);
  //arb_log(arb_ln_pi,arb_pi,prec);
  //arb_init(arb_2_pi);
  //arb_mul_2exp_si(arb_2_pi,arb_pi,1);
  arf_init(ip_tmp2);
  arf_init(ip_tmp1);
  arf_init(ip_tmp2);
  arf_init(ip_tmp3);
  arf_init(ip_tmp4);

  //arf_init(two_101);
  //arf_set_ui(two_101,1);
  //arf_mul_2exp_si(two_101,two_101,OP_ACC);
  // we assume converting mpz_t to long unsigned int
  // gives us 64 bits
  if(sizeof(long unsigned int)!=8)
    fatal_error("Expecting long int to be 64 bits.");

  arb_init(A);
  arb_init(A1);
  acb_init(c_tmp);
  arb_init(msinc_tmp);
  arb_init(mip_tmp);
  arb_init(minter_tmp);
  arb_init(minter_tmp1);
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

  acb_set_err(gtwiderr,gtwiderr_d);
  acb_set_err(Gtwiderr,Gtwiderr_d);
  acb_set_err(tayerr,tayerr_d);
  acb_set_err(Fmaxerr,Fmaxerr_d);
  arb_set_err(intererr,intererr_d);

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
  acb_set_err(fhatsumerr,fhatsumerr_d);
  acb_set_err(fhattwiderr,fhattwiderr_d);
  acb_set_err(ftwiderr,ftwiderr_d);
  arb_set_ui(A,1);
  arb_div_d(A,A,one_over_A,prec); // don't use A as a spare variable anymore
  arb_set_ui(A1,1);
  arb_div_d(A1,A1,one_over_A1,prec);
  calc_buck(prec);
  arb_init(atan1_tmp1);
  arb_init(atan1_tmp2);

  arb_clear(init_tmp1);
  arb_clear(sqrt_pi);
  arb_clear(half_log_2_pi);
  arb_clear(pi_hi);
  arb_clear(arb_2_pi);
}

// run with a new value of t0
// just need to set g_vec and n_vec up
// rest was already done by init
void re_init(double t0, int64_t prec)
{
  static bool init=false;
  static arb_t init_tmp1,arb_2_pi;
  if(!init)
    {
      init=true;
      arb_init(init_tmp1);
      arb_init(arb_2_pi);
      arb_const_pi(arb_2_pi,prec);
      arb_mul_2exp_si(arb_2_pi,arb_2_pi,1);
    }
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

  // step (2)
  for(i=0;i<N1;i++)
    acb_add(G_vec[i],g_vec[i],gtwiderr,prec); // G=g~
  // step (3)
  fft(G_vec,N1,ws1_f,prec);
  for(i=0;i<=N1/2;i++)
    {
      acb_div_arb(G_vec[i],G_vec[i],A1,prec);
      acb_add(G_vec[i],G_vec[i],Gtwiderr,prec); // step (4)
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
void my_convolve (acb_t *res, acb_t *v1, acb_t *v2, uint64_t n, acb_t *ws_r, acb_t *ws_f,int64_t prec)
{
  int i;
  fft(v1,n,ws_r,prec);
  printf("first fft of convolve done.\n");fflush(stdout);
  acb_t foo;
  acb_init(foo);
  for(i=0;i<n;i++)
    acb_add(foo,foo,v2[i],prec);
  printf("dummy add done.\n");fflush(stdout);
  fft(v2,n,ws_r,prec);
  printf("second fft of convolve done.\n");fflush(stdout);
  for(i=0;i<n;i++)
    acb_mul(res[i],v1[i],v2[i],prec);
  fft(res,n,ws_f,prec); // use ws_f so this is an iFFT
  printf("reverse fft of convolve done.\n");fflush(stdout);
  for(i=0;i<n;i++)
    acb_div_ui(res[i],res[i],n,prec); // normalise the iFFT
}

// f_vec+=G(k)(x+u_m)/k!*S^(k)_m
void do_conv (int k, int64_t prec)
{
  int i;

  make_skn(k,prec);    
  printf("make_skn done.\n");fflush(stdout);
  if(k==0) // convolve straight into f_vec
    {
      my_convolve(f_vec,skn_vec,G_vec,N1,ws1_r,ws1_f,prec);
      printf("k=0 0'th convolution = ");acb_printd(f_vec[0],20);printf("\n");fflush(stdout);
    }
  else // convolve into G_vec, then add into f_vec 
    {
      my_convolve(G_vec,skn_vec,G_vec,N1,ws1_r,ws1_f,prec);
      printf("k=%lu 0'th convolution = ",k);acb_printd(G_vec[0],20);printf("\n");fflush(stdout);
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
  static arb_t arb_two_pi_B,sinc_tmp,msin,mcos;
  if(!init)
    {
      init=true;
      arb_init(sinc_tmp);
      arb_init(msin);
      arb_init(mcos);
      arb_init(arb_two_pi_B);
      arb_const_pi(msin,prec);
      arb_div_d(arb_two_pi_B,msin,INTER_A,prec);
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
void arb_inter_t(arb_ptr f_res, arb_ptr fd_res, acb_t *f_vec, arb_ptr t_ptr, int64_t prec)
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

void acb_inter_t(arb_ptr f_res, acb_t *f_vec, arb_ptr t_ptr, int64_t prec)
{
  arb_inter_t(f_res, NULL, f_vec, t_ptr, prec);
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
      arb_inter_t(new_guess,newton_df,f_vec,old_guess,prec); // new=f(t) df=f'(t)
      arb_div(newton_df,new_guess,newton_df,prec);  // f/f'
      arb_div_d(newton_df,newton_df,one_over_A,prec); // f/f' normalised
      arb_sub(old_guess,old_guess,newton_df,prec); // old-f/f'
    }
  arb_set(new_guess,old_guess);
}
*/


int main(int argc, char **argv)
{
  int64_t prec,num_its,it;
  int64_t i,j,k,int_step,last_max=0;
  double t0,t1,step,t;
  fpos_t fpos;
  arb_t tmp,arb_2_pi;
  acb_t omega;
  time_t last_time=time(NULL),this_time;

  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=5)
    {
      printf("Usage:- %s <prec in bits> <t0> <# interations> <step size>.\n",argv[0]);
      return 0;
    }

  prec=atol(argv[1]);

  t0=atof(argv[2]);
  num_its=atol(argv[3]);
  if(num_its<=0)
    {
      printf("Number of iterations must be > 0 \n");
      return 0;
    }
  step=atof(argv[4]);
  if(step<=0)
    {
      printf("step must be > 0.\n");
      return 0;
    }
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
    {
      printf("step*A must be integral %f.\n",step/one_over_A);
      return 0;
    }
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

  printf("Aiming to do %ld iteration(s) starting at t=%13.11e.\n",num_its,t0-step/2.0);

  arb_init(arb_2_pi);
  arb_const_pi(arb_2_pi,prec);
  arb_mul_2exp_si(arb_2_pi,arb_2_pi,1);
  for(it=0;it<num_its;it++,t0+=step)
    {
      if(it==0)
	{
	  init(t0,prec);
	  arb_init(tmp);
	  acb_init(omega);
	  arb_div_ui(tmp,arb_2_pi,N,prec);
	  arb_sin_cos(acb_imagref(omega),acb_realref(omega),tmp,prec);
	  arb_clear(tmp);
	}
      else
	re_init(t0,prec);

      this_time=time(NULL);
      printf("Time to (re-)initialise = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;

      for(k=0;k<K;k++)
	{
	  G_k(k,prec);
	  printf("G_k done.\n");fflush(stdout);
	  do_conv(k,prec);
	  printf("Convolve done.\n");fflush(stdout);
	}
      this_time=time(NULL);
      printf("Time to convolve = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;

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


      this_time=time(NULL);
      printf("Time for final iFFT = %G seconds.\n",difftime(this_time,last_time));
      last_time=this_time;



      // remove the Gaussian from the region of interest
      for(j=0,i=N/2+1;i<=N/2+TURING_LEN+int_step/2+INTER_SPACING*Ns;j++,i++)
	{
	  arb_div_d(acb_realref(f_vec[i]),acb_realref(f_vec[i]),B,prec);
	  arb_add(acb_realref(f_vec[i]),acb_realref(f_vec[i]),acb_realref(ftwiderr),prec);
	  arb_mul(acb_realref(f_vec[i]),acb_realref(f_vec[i]),exps1[j],prec);
	}
      for(j=0,i=N/2-1;i>=N/2-TURING_LEN-int_step/2-INTER_SPACING*Ns;j++,i--)
	{
	  arb_div_d(acb_realref(f_vec[i]),acb_realref(f_vec[i]),B,prec);
	  arb_add(acb_realref(f_vec[i]),acb_realref(f_vec[i]),acb_realref(ftwiderr),prec);
	  arb_mul(acb_realref(f_vec[i]),acb_realref(f_vec[i]),exps1[j],prec);
	}
      
      arb_div_d(acb_realref(f_vec[N/2]),acb_realref(f_vec[N/2]),B,prec);
      arb_add(acb_realref(f_vec[N/2]),acb_realref(f_vec[N/2]),acb_realref(ftwiderr),prec);

      last_max=turing(f_vec,N/2-int_step/2,t0-int_step/2*one_over_A,N/2+int_step/2,t0+int_step/2*one_over_A,last_max,prec);
      this_time=time(NULL);
      printf("Time to isolate zeros = %G seconds.\n",difftime(this_time,last_time));
      fflush(stdout);

      arb_clear(tmp);
      arb_clear(arb_2_pi);
      for(i=0,t=one_over_A;i<TURING_LEN+int_step/2+Ns*INTER_SPACING;t+=one_over_A,i++)
	arb_clear(exps1[i]);

      last_time=this_time;

    }
}
