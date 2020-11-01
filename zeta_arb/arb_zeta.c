/*
Code to check RH

See Booker: Artin's Conjecture, Turing's Method, and the Riemann Hypothesis - Exp. Math.
    Platt: Isolating some non-trivial zeros of zeta - Math. Comp.
    Platt & Trudgian: The Riemann hypothesis is true up to 3\cdot 10^{12} - To appear 
*/

// arb_zeta.c
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
#include "win_zeta.h"
#include "arb_fft.h"
#include "inter.h"
#include "turing.h"


int64_t HI_PREC;
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
arb_t misin;


// inits and sets z to [-d,d]+[-d,d]i
void acb_set_err(acb_t z, double d)
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
int conj_order(int i)
{
  return N1-i;
}

// on entry s_vec[i] contains log((i+1)Pi)/2Pi
// on exit                    "" - u_m
void calc_buck(int64_t prec)
{
  arb_t tmp,arb_B;
  arb_init(tmp);
  arb_init(arb_B);
  arb_set_d(arb_B,B);
  int i;
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
  arb_clear(tmp);
  arb_clear(arb_B);
}

arb_t pi_hi,half_log_2_pi;

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
  arb_t init_tmp1,sqrt_pi,pi_hi;
  arb_init(init_tmp1);
  arb_init(sqrt_pi);
  arb_init(pi_hi);
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
  //arf_mul_2exp_si(two_101,two_101,OP_ACC);
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
  arb_init(half_log_2_pi);
  arb_const_pi(sqrt_pi,HI_PREC);
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
  acb_set_err(fhatsumerr,fhatsumerr_d);
  acb_set_err(fhattwiderr,fhattwiderr_d);
  acb_set_err(ftwiderr,ftwiderr_d);
  arb_set_ui(A,1);
  arb_div_d(A,A,one_over_A,prec); // don't use A as a spare variable anymore
  arb_set_ui(A1,1);
  arb_div_d(A1,A1,one_over_A1,prec);
  calc_buck(prec);
  arb_clear(init_tmp1);
  arb_clear(sqrt_pi);
  arb_clear(pi_hi);
}

// run with a new value of t0
// just need to set g_vec and n_vec up
// rest was already done by init
void re_init(double t0, int64_t prec)
{
  static bool init=false;
  static arb_t init_tmp1;
  if(!init)
    {
      init=true;
      arb_init(init_tmp1);
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

  for(i=0;i<N1;i++)
    acb_add(G_vec[i],g_vec[i],gtwiderr,prec); // G=g~

  fft(G_vec,N1,ws1_f,prec);

  for(i=0;i<=N1/2;i++)
    {
      acb_div_arb(G_vec[i],G_vec[i],A1,prec);
      acb_add(G_vec[i],G_vec[i],Gtwiderr,prec);
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
  int i;
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

// convolve v1 and v2 into res
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

int main(int argc, char **argv)
{
  long int prec,num_its,it;
  long int i,j,k,int_step,last_max=0;
  double t0,step,t;
  arb_t tmp;
  bool first=true;
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
	}
      else
	re_init(t0,prec);

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

      // upsample by "zero" padding
      for(i=N1/2;i<=N/2;i++)
	acb_set(f_vec[i],Fmaxerr);

      // add various error terms
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

      // use Turing's method to check all zeros accounted for
      last_max=turing(f_vec,N/2-int_step/2,t0-int_step/2*one_over_A,N/2+int_step/2,t0+int_step/2*one_over_A,last_max,prec);

      this_time=time(NULL);
      printf("Time to isolate zeros = %G seconds.\n",difftime(this_time,last_time));
      fflush(stdout);
      last_time=this_time;

    }
}
