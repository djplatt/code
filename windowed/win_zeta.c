//
// win_zeta.c
//
// Windowed FFT based Zeta-function calculator
// See Booker - Artin, Turing ....
//
// Vesrion 1.0 Initial implementation
//
// Created: 28 September 2010
// Last Modified: 19 October 2010
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

#define N ((int) 1<<19)
#define one_over_A ((double) 21.0/512.0)
#define h ((double) 180.911)
#define B ((double) N*one_over_A)
#define M ((int) 88000)

#define K ((int) 36)
#define gtwiderr_d ((double) 1.5e-62)
#define Gtwiderr_d ((double) 1e-104) // actually much smaller
#define fhatsumerr_d ((double) 1e-104)
#define tayerr_d ((double) 10.0e-64)
#define fhattwiderr_d ((double) 1e-104) // actually much smaller
#define ftwiderr_d ((double) 1e-104) // ditto

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

mpfi_t A,two_pi_t[N],s_vec[M],sk_vec[M],exps[N],sqrts[M],logisqrtpi[M];
mpfi_c_t g_vec[N],G_vec[N],f_vec[N],skn_vec[N],ws_r[N/2],ws_f[N/2],n_vec[M];
mpfi_c_t gtwiderr,Gtwiderr,c_tmp,fhatsumerr,fhattwiderr,ftwiderr,tayerr;
int buck[M]; // which bucket does log(nsqrt(Pi)/2pi) go in.

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
  return(N-i);
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
  assert(buck[M-1]<N/2); // ditto
  for(i=0;i<M;i++)
    buck[i]=conj_order(buck[i]);
  //printf("buck[0]=%d\nbuck[M-1]=%d\n",buck[0],buck[M-1]);
  //printf("n=1 term is in bucket %d\n",buck[0]);
  mpfi_clear(tmp);
}
      
      
// on exit g[n]=g((n-N/2)/A;0)
void init(double t0)
{
  int i;
  double t;
  mpfi_init(A);
  mpfi_c_init(c_tmp);
  for(i=0,t=-N/2*one_over_A;i<N;i++,t+=one_over_A)
    {
      //if ((i&4095)==0)
      //printf("t=%e\n",t);
      mpfi_init(two_pi_t[i]);
      mpfi_mul_d(two_pi_t[i],mpfi_2_pi,-t);
      mpfi_c_init(g_vec[i]);
      mpfi_c_init(G_vec[i]);
      mpfi_c_init(f_vec[i]);
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
  
  
  initfft(N,ws_r);
  for(i=0;i<N/2;i++)
    {
      mpfi_c_init(ws_f[i]);
      mpfi_c_conj(ws_f[i],ws_r[i]);
    }
  set_err(gtwiderr,gtwiderr_d);
  set_err(Gtwiderr,Gtwiderr_d);
  set_err(tayerr,tayerr_d);
  mpfi_set_ui(A,M);
  mpfi_sqrt(A,A);
  mpfi_mul_ui(A,A,2);
  mpfi_sub_ui(A,A,1); // 2sqrt(M)-1
  mpfi_c_mul_i(tayerr,tayerr,A);
  printf("Total Taylor Error set to ");mpfi_printn(tayerr->re,30);
  set_err(fhatsumerr,fhatsumerr_d);
  set_err(fhattwiderr,fhattwiderr_d);
  set_err(ftwiderr,ftwiderr_d);
  mpfi_set_ui(A,1);
  mpfi_div_d(A,A,one_over_A); // don't use A as a spare variable anymore
  calc_buck();

}

// run with a new value of t0
// just need to set g_vec and n_vec up
void re_init(double t0)
{
  int i;
  double t;
  //mpfi_init(A);
  for(i=0,t=-N/2*one_over_A;i<N;i++,t+=one_over_A)
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

  for(i=0;i<N;i++)
    mpfi_c_add(G_vec[i],g_vec[i],gtwiderr); // G=g~

  //simple_idft(0);
  //simple_idft(965);
  //simple_idft(N-1);

  //nifft(G_vec,N,ws_r);
  fft(G_vec,N,ws_f);

  for(i=0;i<=N/2;i++)
    {
      mpfi_c_div_i(G_vec[i],G_vec[i],A);
      mpfi_c_add(G_vec[i],G_vec[i],Gtwiderr);
      if(i&1)
	mpfi_c_neg(G_vec[i],G_vec[i]);
    }
  for(i=N/2+1;i<N;i++)
    mpfi_c_set_ui(G_vec[i],0,0);

  if(k<K-1)
    for(i=0;i<N;i++)
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
  for(i=0;i<N;i++)
    mpfi_c_zero(skn_vec[i]);
      
  switch(k)
    {
    case 0: // set sk_vec = s_vec, skn_vec=n_vec
      for(i=0;i<M;i++)
	{
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
    
  if(k==0)
    my_convolve(f_vec,skn_vec,G_vec,N,ws_r,ws_f);
  else
    {
      my_convolve(G_vec,skn_vec,G_vec,N,ws_r,ws_f);
      for(i=0;i<=N/2;i++)
	mpfi_c_add(f_vec[i],f_vec[i],G_vec[i]);
    }
}

int main(int argc, char **argv)
{
  int prec;
  int i,j,k,int_step;
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
  t1=atof(argv[3]);
  step=atof(argv[4]);
  int_step=step/one_over_A;
  assert((int_step&1)==0); // must be even
  assert(step/one_over_A-(double) int_step ==0.0); // must be an exact number of steps
  /*
  assert(t0>=1000000.0);
  HI_BERNOULLI_K=8;
  if(t0>=1e7)
    HI_BERNOULLI_K=7;
  if(t0>=1e8)
    HI_BERNOULLI_K=6;
  if(t0>=1e9)
    HI_BERNOULLI_K=5;
  */
  printf("Aiming to get %d values per run.\n",int_step);

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
	  printf("Processing k=%d\n",k);
	  /*
	  for(i=0;i<N;i+=N/64)
	    {
	      printf("g_vec[%d]=",i);
	      mpfi_c_printn(g_vec[i],30);
	    }
	  */
	  G_k(k);
	  /*
	  for(i=0;i<N;i++)
	    if(mpfi_rel_error(G_vec[i]->re)<-30)
	      break;
	  for(j=i;j<i+32;j++)
	    {
	      printf("G_vec[%d]=",j);
	      mpfi_c_printn(G_vec[j],30);
	    }
	  exit(0);
	  */
	  do_conv(k);
	  for(i=0;i<N/32;i+=N/(32*32))
	    {
	      printf("f_vec[%d]=",i);
	      mpfi_c_printn(f_vec[i],30);
	    }
	  exit(0);
	}
      printf("convolutions finished\n");
      system("date");

      for(i=0;i<=N/2;i++)
	{
	  mpfi_c_add(f_vec[i],f_vec[i],fhatsumerr);
	  mpfi_c_add(f_vec[i],f_vec[i],tayerr);
	  mpfi_c_add(f_vec[i],f_vec[i],fhattwiderr);
	}

      for(i=N/2+1;i<N;i++)
	mpfi_c_conj(f_vec[i],f_vec[N-i]);

      for(i=1;i<N;i+=2)
	mpfi_c_neg(f_vec[i],f_vec[i]);

      printf("Final iFFT\n");
      fft(f_vec,N,ws_r);
      printf("Final iFFT finished.\n");
      system("date");

      /*
      for(i=N/2-int_step/2,t=-step/2.0;i<=N/2+int_step/2;i++,t+=one_over_A)
	{
	  mpfi_div_d(f_vec[i]->re,f_vec[i]->re,B);
	  mpfi_add(f_vec[i]->re,f_vec[i]->re,ftwiderr->re);
       
	  // now get rid of the Gaussian
	  mpfi_set_d(tmp,t);
	  mpfi_div_d(tmp,tmp,h);
	  mpfi_sqr(tmp,tmp);
	  mpfi_div_ui(tmp,tmp,2);
	  mpfi_exp(tmp,tmp);
	  mpfi_mul(tmp,f_vec[i]->re,tmp);

	  // this is it!
	  printf("f(%10f)=",t0+t);mpfi_printn(tmp,30);
	}
      */

      
      for(i=0,t=-one_over_A*N/2;i<=N/2;i+=N/256,t+=one_over_A*N/256)
	{
	  mpfi_div_d(f_vec[i]->re,f_vec[i]->re,B);
	  mpfi_add(f_vec[i]->re,f_vec[i]->re,ftwiderr->re);
       
	  // now get rid of the Gaussian
	  mpfi_set_d(tmp,t);
	  mpfi_div_d(tmp,tmp,h);
	  mpfi_sqr(tmp,tmp);
	  mpfi_div_ui(tmp,tmp,2);
	  mpfi_exp(tmp,tmp);
	  mpfi_mul(tmp,f_vec[i]->re,tmp);
	  printf("%6d: Rel error = %4d Abs error = %4d f(%18.10f)=",
		 N/2-i,mpfi_rel_error(tmp), mpfi_abs_error(tmp),t+t0);
	  mpfi_printn(tmp,40);
	}
      

      t0+=step;
    }
  fclose(outfile);
}
