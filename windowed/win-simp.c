//
// win-simp.c
//
// Windowed FFT based Zeta-function calculator
//
//
// Vesrion 1.0 Initial implementation
//
// Last Modified: 21 September 2010
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
#define MAX_N (24) // largest length of FFT array = 2^(MAX_N)

#define one_over_A ((double) 6143.0/131072.0)
#define h ((double) 15.0)
#define h_sqr ((double) h*h)
#define H ((double) 12347.0/65536.0)
#define N (51)
#define N0 (1<<15)
#define B ((double) (N0*one_over_A))
#define K (15)
#define M (10)

void print_usage()
{
  printf("Usage:- win-simp prec t0.\n");
  exit(1);
}

mpfi_t sqrt_2_pi;

void G(mpfi_ptr res, mpfi_ptr u)
{
  mpfi_mul(res,mpfi_pi,u);
  mpfi_mul_d(res,res,h);
  mpfi_sqr(res,res);
  mpfi_mul_d(res,res,-2.0);
  mpfi_exp(res,res);
  mpfi_mul(res,res,sqrt_2_pi);
  mpfi_mul_d(res,res,h);
}

mpfi_t f_hat_tmp,f_hat_tmp1;
mpfi_t logs[M+1],sqrts[M+1];
mpfi_c_t nits[M+1],f_hat_tmp2;

void init_f_hat(double t0)
{
  int i;
  mpfi_t tmp;
  mpfi_init(tmp);
  mpfi_init(f_hat_tmp);
  mpfi_init(f_hat_tmp1);
  mpfi_c_init(f_hat_tmp2);
  for(i=1;i<=M;i++)
    {
      mpfi_init(logs[i]);
      mpfi_init(sqrts[i]);
      mpfi_c_init(nits[i]);
      mpfi_set_ui(logs[i],i);
      mpfi_sqrt(sqrts[i],logs[i]);
      mpfi_log(logs[i],logs[i]);
      mpfi_mul_d(tmp,logs[i],-t0);
      mpfi_cos(nits[i]->re,tmp);
      mpfi_sin(nits[i]->im,tmp);
      mpfi_div(logs[i],logs[i],mpfi_2_pi);
      //printf("log %d ",i);mpfi_print(logs[i]);
      //printf("sqrt %d ",i);mpfi_print(sqrts[i]);
      //printf("nits %d ",i);mpfi_c_print(nits[i]);
    }
  mpfi_clear(tmp);
}

// m<0
void f_hat_bar(mpfi_c_ptr res, int m)
{
  int i;
  mpfi_set_d(f_hat_tmp1,m);
  mpfi_div_d(f_hat_tmp1,f_hat_tmp1,B);
  mpfi_c_set_ui(res,0,0);
  for(i=1;i<=M;i++)
    {
      mpfi_sub(f_hat_tmp,f_hat_tmp1,logs[i]); // m/B-log(i)/2Pi
      G(f_hat_tmp,f_hat_tmp);                  // G("")
      mpfi_div(f_hat_tmp,f_hat_tmp,sqrts[i]);  // G("")/sqrt(i)
      mpfi_c_conj(f_hat_tmp2,nits[i]);
      mpfi_c_mul_i(f_hat_tmp2,f_hat_tmp2,f_hat_tmp); // n^{it0}n^{-1/2}G(u+
      mpfi_c_inc(res,f_hat_tmp2);
    }

  if((m&255)==0)
    {
      mpfi_c_div(f_hat_tmp2,f_hat_tmp2,res);
      mpfi_c_norm(f_hat_tmp1,f_hat_tmp2);
      printf("Last term ratio sqr=");
      mpfi_printn(f_hat_tmp1,10);
    }
  mpfi_c_conj(res,res);
}

// x=m/B
void f_hat(mpfi_c_ptr res, int m)
{
  int i;
  mpfi_set_d(f_hat_tmp1,m);
  mpfi_div_d(f_hat_tmp1,f_hat_tmp1,B);
  mpfi_c_set_ui(res,0,0);
  for(i=1;i<=M;i++)
    {
      mpfi_add(f_hat_tmp,f_hat_tmp1,logs[i]); // m/B+log(i)/2Pi
      G(f_hat_tmp,f_hat_tmp);                  // G("")
      mpfi_div(f_hat_tmp,f_hat_tmp,sqrts[i]);  // G("")/sqrt(i)
      mpfi_c_mul_i(f_hat_tmp2,nits[i],f_hat_tmp); // n^{-it0}n^{-1/2}G(u+
      mpfi_c_inc(res,f_hat_tmp2);
    }

  if((m&255)==0)
    {
      mpfi_c_div(f_hat_tmp2,f_hat_tmp2,res);
      mpfi_c_norm(f_hat_tmp1,f_hat_tmp2);
      printf("Last term ratio sqr=");
      mpfi_printn(f_hat_tmp1,10);
    }

}

int main(int argc, char **argv)
{

  int i,prec;
  double t0;
  mpfi_c_t *Fft_vec,*ws;
  mpfi_t A,t;

  if(argc!=3)
    print_usage();

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    print_usage();
  printf("Running at %d bits of precision.\n",prec);
    


  t0=atof(argv[2]);
  if(t0<=0.0)
    print_usage();
  printf("t0 set to %f\n",t0);
  
  if(!(Fft_vec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*N0)))
    {
      printf("Fatal error allocating memory for Fft_vec. Exiting.\n");
      exit(FAILURE);
    }

  if(!(ws=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*N0/2)))
    {
      printf("Fatal error allocating memory for ws. Exiting.\n");
      exit(FAILURE);
    }

  mpfi_c_setup(prec);

  mpfi_init(A);
  mpfi_set_ui(A,1);
  mpfi_div_d(A,A,one_over_A);
  printf("1/A=%20.18e\n",one_over_A);
  printf("A set to ");mpfi_printn(A,10);
  printf("N0=%d\n",N0);
  printf("h=%10.8e\n",h);
  printf("B=%10.8e\n",B);

  mpfi_init(sqrt_2_pi);
  mpfi_sqrt(sqrt_2_pi,mpfi_2_pi);


  initfft(N0,ws);
  init_f_hat(t0);
  

  for(i=0;i<=N0/2;i++)
    {
      mpfi_c_init(Fft_vec[i]);
      f_hat(Fft_vec[i],i);
      if(i&1)
	mpfi_c_neg(Fft_vec[i],Fft_vec[i]);
    }
  for(i=N0/2+1;i<N0;i++)
    {
      mpfi_c_init(Fft_vec[i]);
      f_hat_bar(Fft_vec[i],i-N0);
      if(i&1)
	mpfi_c_neg(Fft_vec[i],Fft_vec[i]);
    }

  
  for(i=0;i<N0;i++)
    if((i&255)==0)
      {
	printf("%6d ",i);
	mpfi_c_print(Fft_vec[i]);
      }
  

  nifft(Fft_vec,N0,ws);
  mpfi_init(t);
  for(i=0;i<N0;i++)
    {
      mpfi_c_div_d(Fft_vec[i],Fft_vec[i],B);
      mpfi_set_d(t,B/2.0-i*one_over_A);
      mpfi_sqr(t,t);
      mpfi_div_d(t,t,h_sqr*2.0);
      mpfi_exp(t,t);
      mpfi_c_mul_i(Fft_vec[i],Fft_vec[i],t);
      if((i&255)==0)
	{
	  printf("exp(t^2/2h^2)=");mpfi_printn(t,10);
	  printf("t=%15.12e ",t0+(i-N0/2)*one_over_A);
	  mpfi_c_printn(Fft_vec[i],30);
	}
    }
  mpfi_clear(t);
  mpfi_clear(A);
}


