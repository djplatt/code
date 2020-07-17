//
// win_zeta1.1.c
//
// Windowed FFT based Zeta-function calculator
// See Booker - Artin, Turing ....
//
// Version 1.0 Initial implementation
//         1.1 Vary A (and N) between g/f
//
// Created: 28 September 2010
// Last Modified: 19 November 2010
//
// DJ Platt
// University of Bristol
//

#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "../includes/mpfi_fft.h"
#include "../includes/fft_defs.h"


#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define debug printf("Reached line number %d\n",__LINE__)
#define bool int

// parameters for t0<=3*10^9
#define T0_MIN (10000.0) // mimimal input t0
#define UPSAM (32) // the rate of upsampling from g to f
#define N ((int) 1<<20) // final FFT length
#define N1 ((int) N/UPSAM) // intermediate FFT length
#define one_over_A ((double) 21.0/4096.0) // output spacing
#define one_over_A1 ((double) one_over_A*UPSAM) // g/G spacing
#define h ((double) 86.148) // Gaussian width for f(t)
#define B ((double) N*one_over_A)
#define M ((int) 103000) // number of terms in F(x) sum

#define K ((int) 44) // number of Taylor terms
#define gtwiderr_d ((double) 3.2e-82)
#define Gtwiderr_d ((double) 1e-104) // actually much smaller
#define fhatsumerr_d ((double) 1.7e-83)
#define tayerr_d ((double) 6.7e-82)
#define fhattwiderr_d ((double) 1e-104) // actually much smaller
#define ftwiderr_d ((double) 1e-104) // ditto
#define Fmaxerr_d ((double) 1e-104) // Max mod of F(N1/(2B)) actually much smaller

#define TURING_WIDTH (42.0)
#define TURING_LEN ((int) TURING_WIDTH/one_over_A)

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
mpfi_c_t g_vec[N1],G_vec[N1],f_vec[N],skn_vec[N1],ws_r[N/2],ws1_r[N1/2],ws1_f[N1/2],n_vec[M];
mpfi_c_t gtwiderr,Gtwiderr,c_tmp,fhatsumerr,fhattwiderr,ftwiderr,tayerr,Fmaxerr;
int buck[M]; // which bucket does log(nsqrt(Pi)/2pi) go in.
mpfi_t *exps1;

void print_usage()
{
  printf("Usage:- win_zeta prec t0 t1 step outfile.\n");
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
mpfr_t tm_mpfr;
mpfi_t tm_t0,tm_t1;
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
      
// on exit g[n]=g((n-N/2)/A;0)
void init(double t0)
{
  int i;
  double t;
  mpfi_init(A);
  mpfi_init(A1);
  mpfi_c_init(c_tmp);
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
      mpfi_init(exps[i]);
      //printf("Taking lngamma of 0.25+%ei\n",(t+t0)/2.0);
      //mpfi_set_d(g_vec[i]->re,-t-t0);mpfi_set_ui(g_vec[i]->im,i);
      mpfi_c_lngamma_hi(g_vec[i],0.25,(t+t0)/2.0);    // lngamma(1/4+(t+t0)/2*I)
      mpfi_mul_d(G_vec[i]->re,mpfi_pi,(t+t0)/4.0); // Pi(t+t0)/4
      mpfi_add(g_vec[i]->re,g_vec[i]->re,G_vec[i]->re);
      mpfi_set_d(G_vec[i]->re,t);
      mpfi_div_d(G_vec[i]->re,G_vec[i]->re,h);
      mpfi_sqr(G_vec[i]->re,G_vec[i]->re);
      mpfi_div_ui(exps[i],G_vec[i]->re,2); //  t^2/(2h^2)
      mpfi_sub(g_vec[i]->re,g_vec[i]->re,exps[i]);
      mpfi_c_exp(g_vec[i],g_vec[i]);
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
      mpfi_init(logisqrtpi[i]);
      mpfi_c_init(n_vec[i]);
      mpfi_mul_ui(s_vec[i],mpfi_sqrt_pi,i+1); // n*sqrt(Pi)
      mpfi_log(s_vec[i],s_vec[i]);        // log (n*sqrt(Pi))
      mpfi_set(logisqrtpi[i],s_vec[i]);
      mpfi_mul_d(n_vec[i]->re,s_vec[i],-t0); // -t0*log(n*sqrt(Pi))
      mpfi_sin(n_vec[i]->im,n_vec[i]->re);
      mpfi_cos(n_vec[i]->re,n_vec[i]->re);  // exp(-i*t0*log(n*sqrt(Pi)))
      mpfi_set_ui(A,1);
      mpfi_div_ui(A,A,i+1);
      mpfi_sqrt(A,A);
      mpfi_init(sqrts[i]);
      mpfi_set(sqrts[i],A);
      mpfi_c_mul_i(n_vec[i],n_vec[i],A);    //  / sqrt(n)
      mpfi_div(s_vec[i],s_vec[i],mpfi_2_pi); // log(n*sqrt(Pi))/(2*Pi)
    }

  // set up the omega vectors for FFTs
  initfft(N,ws_r); // ws_r[i]=e(i/N)
  for(i=0;i<N1/2;i++)
    {
      mpfi_c_init(ws1_f[i]);
      mpfi_c_init(ws1_r[i]);
      mpfi_c_set(ws1_r[i],ws_r[i*UPSAM]); // ws1_r[i]=e(i/N1)
      mpfi_c_conj(ws1_f[i],ws1_r[i]);     // ws1_f[i]=e(-i/N1)
    }

  set_err(gtwiderr,gtwiderr_d);
  set_err(Gtwiderr,Gtwiderr_d);
  set_err(tayerr,tayerr_d);
  set_err(Fmaxerr,Fmaxerr_d);
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
  for(i=0,t=-N1/2*one_over_A1;i<N1;i++,t+=one_over_A1)
    {
      //mpfi_set_d(g_vec[i]->re,-t-t0);mpfi_set_ui(g_vec[i]->im,i);
      mpfi_c_lngamma_hi(g_vec[i],0.25,(t+t0)/2.0);    // lngamma(1/4+(t+t0)/2*I)
      mpfi_mul_d(c_tmp->re,mpfi_pi,(t+t0)/4.0);
      mpfi_add(g_vec[i]->re,g_vec[i]->re,c_tmp->re);
      mpfi_sub(g_vec[i]->re,g_vec[i]->re,exps[i]);
      mpfi_c_exp(g_vec[i],g_vec[i]);
    }
  for(i=0;i<M;i++)
    {
      //mpfi_mul_ui(n_vec[i]->re,mpfi_sqrt_pi,(i+1)); // ((i+1)*sqrt(Pi))
      //mpfi_log(n_vec[i]->re,n_vec[i]->re);
      mpfi_mul_d(n_vec[i]->re,logisqrtpi[i],-t0);
      mpfi_sin(n_vec[i]->im,n_vec[i]->re);
      mpfi_cos(n_vec[i]->re,n_vec[i]->re);  // exp(-i*t0*log((i+1)*sqrt(Pi)))
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

typedef int sign_t;
#define POS (0)
#define NEG (1)
#define UNK (-1)

inline sign_t sign(mpfi_ptr x)
{
  if(mpfi_is_neg(x))
    return(NEG);
  if(mpfi_is_pos(x))
    return(POS);
  return(UNK);
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

int num_zeros(mpfi_c_t *f_vec, int start, int end)
{
  int i=start,res=0;
  sign_t last_sign=sign(f_vec[i]->re),this_sign;
  while(last_sign==UNK) // skip initial section of UNKs, will break if all UNKs!
    last_sign=sign(f_vec[++i]->re);
  for(i++;i<=end;i++)
    {
      this_sign=sign(f_vec[i]->re);
      if((this_sign==last_sign)||(this_sign==UNK))
	continue;
      res++;  // A sign change
      last_sign=this_sign;
    }
  return(res);
}

// int Im loggamma(1/4+it)
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

double Nleft_int(int t0_ptr, int t1_ptr, mpfi_c_t *f_vec, double delta)
{
  int res=0;
  int ptr=t0_ptr,last_ptr;
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

double Nright_int(int t0_ptr, int t1_ptr, mpfi_c_t *f_vec, double delta)
{
  int res=0;
  int ptr=t0_ptr;
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
  assert(mpfi_cmp_d(t,168.0*M_PI+0.1)>0);
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
int turing_max(mpfi_c_t *f_vec, int a_ptr, double a, int b_ptr, double b)
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
  mpfi_get_right(tm_mpfr,tm_res);
  return(mpfr_get_si(tm_mpfr,GMP_RNDD));
}

int turing_min(mpfi_c_t *f_vec, int a_ptr, double a, int b_ptr, double b)
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
  mpfi_get_left(tm_mpfr,tm_res);
  return(mpfr_get_si(tm_mpfr,GMP_RNDU));
}

// Use Turing method to estimate and then check number of zeros
int turing(mpfi_c_t *f_vec, int a_ptr, double a, int b_ptr, double b, int last_max)
{
  // etimate maximum for N(b)
  int min_lo,max_hi=turing_max(f_vec,b_ptr,b,b_ptr+TURING_LEN,b+TURING_WIDTH);
  if(last_max==0) // No previous run or previous run failed, so estimate min N(a)
    min_lo=turing_min(f_vec,a_ptr-TURING_LEN,a-TURING_WIDTH,a_ptr,a);
  else // If previous run succeeded, zeros to height N(a) is known
    min_lo=last_max;
  int num_exp=max_hi-min_lo;
  int num_found=num_zeros(f_vec,a_ptr,b_ptr);
  if(num_exp==num_found)
    {
      printf("All zeros found (without stats) in region %f to %f.\n",a,b);
      return(max_hi); // use this for next iteration
    }
  int stats=num_stat_pts(f_vec,a_ptr,b_ptr); 
  if(num_exp==(num_found+stats))
    {
      printf("All zeros found (with stats) in region %f to %f.\n",a,b);
      return(max_hi);
    }
  printf("Missed %d zeros in region %f to %f even after %d stats.\n",num_exp-num_found-stats,a,b,stats);
  return(0);
}

int main(int argc, char **argv)
{
  int prec;
  int i,j,k,int_step,last_max=0;
  double t0,t1,step,t;
  FILE *outfile;
  mpfi_t tmp;
  bool first=TRUE;

  if(argc!=6)
    print_usage();

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    print_usage();
  printf("Running at %d bits of precision.\n",prec);
    
  mpfi_c_setup(prec);
  mpfi_init(tmp);
  /*
  mpfi_c_init(c_tmp);
  mpfi_c_lngamma_hi(c_tmp,0.25,500000);
  printf("ln_gamma(0.25+500000i)=");
  mpfi_c_printn(c_tmp,100);
  printf("abs error re=%d\n",mpfi_abs_error(c_tmp->re));
  printf("rel error re=%d\n",mpfi_rel_error(c_tmp->re));

  exit(0);
  */
  outfile=fopen(argv[5],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[5]);
      exit(1);
    }

  t0=atof(argv[2]);
  assert(t0>=T0_MIN);
  t1=atof(argv[3]);
  step=atof(argv[4]);
  int_step=step/one_over_A;
  assert((int_step&1)==0); // must be even
  assert(step/one_over_A-(double) int_step ==0.0); // must be an exact number of steps
  printf("Aiming to get %d values per run.\n",int_step+1);
  exps1=(mpfi_t *) malloc(int_step/2*sizeof(mpfi_t));
  if(!exps1)
    {
      printf("Failed to allocate memory for exps1. Exiting.\n");
      exit(0);
    }
  for(i=0,t=one_over_A;i<int_step/2;t+=one_over_A,i++)
    {
      mpfi_init(exps1[i]);
      mpfi_set_d(exps1[i],t);
      mpfi_div_d(exps1[i],exps1[i],h);
      mpfi_sqr(exps1[i],exps1[i]);
      mpfi_div_ui(exps1[i],exps1[i],2);
      mpfi_exp(exps1[i],exps1[i]);
    }

  //main loop
  while(t0<=t1)
    {
      printf("Running centred at t0=%f.\n",t0);
      if(first)
	{
	  first=FALSE;
	  printf("Calling init.\n");
	  system("date");
	  init(t0);
	  printf("Init finished.\n");
	  system("date");
	}
      else
	{
	  printf("Calling re-init.\n");
	  system("date");
	  re_init(t0);
	  printf("re-init finished.\n");
	  system("date");
	}

      for(k=0;k<K;k++)
	{
	  G_k(k);
	  do_conv(k);
	}
      printf("convolutions finished\n");
      system("date");

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

      printf("Final iFFT\n");
      fft(f_vec,N,ws_r);
      printf("Final iFFT finished.\n");
      system("date");

      
      for(j=0,i=N/2+1;i<=N/2+int_step/2;j++,i++)
	{
	  mpfi_div_d(f_vec[i]->re,f_vec[i]->re,B);
	  mpfi_add(f_vec[i]->re,f_vec[i]->re,ftwiderr->re);
	  mpfi_mul(f_vec[i]->re,f_vec[i]->re,exps1[j]);
	}
      for(j=0,i=N/2-1;i>=N/2-int_step/2;j++,i--)
	{
	  mpfi_div_d(f_vec[i]->re,f_vec[i]->re,B);
	  mpfi_add(f_vec[i]->re,f_vec[i]->re,ftwiderr->re);
	  mpfi_mul(f_vec[i]->re,f_vec[i]->re,exps1[j]);
	}
      mpfi_div_d(f_vec[N/2]->re,f_vec[N/2]->re,B);
      mpfi_add(f_vec[N/2]->re,f_vec[N/2]->re,ftwiderr->re);
      /*
      for(i=N/2-int_step/2;i<=N/2+int_step/2;i+=int_step/64)
      */

      last_max=turing(f_vec,N/2-int_step/2,t0-int_step/2*one_over_A,N/2+int_step/2,t0+int_step/2*one_over_A,last_max);
      t0+=step;
    }
  fclose(outfile);
}
