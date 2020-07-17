// window1/test.c
//
// try out niave windowed FFT algorithm on q=5
// chi(1)=1, chi(2)=-1, chi(3)=-1, chi(4)=1
// an even character
#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "../includes/mpfi_fft.h"
#include "../includes/fft_defs.h"
#include "../low_q/f_defs.h"

#define TRUE (0==0)
#define FALSE (1==0)
#define FAILURE (1)
#define debug printf("Reached line number %d\n",__LINE__)
#define bool int
#define one_over_A ((double) 5.0/64.0)

void print_usage()
{
  printf("Usage:- test <prec> <n> <sigma> <t0>\n");
  printf("Size of FFT will be 2^n.\n");
  printf("sigma >1\n");
  exit(1);
}

inline double calc_h(double B)
{
  return(B/20.0);
}

// M is the last term summed
inline void sum_error(mpfi_c_ptr err, unsigned int M, double sigma)
{
  mpfi_set_ui(err->re,M);
  mpfi_log(err->re,err->re);
  mpfi_neg(err->re,err->re);
  mpfi_mul_d(err->re,err->re,sigma-1.0);
  mpfi_exp(err->re,err->re);
  mpfi_div_d(err->re,err->re,sigma);
  mpfi_neg(err->im,err->re);
  mpfi_put(err->re,err->im);
  mpfi_set(err->im,err->re);
}

mpfi_c_t do_sum_t;
mpfi_t n_sigma;
void init_do_sum()
{
  mpfi_c_init(do_sum_t);
  mpfi_init(n_sigma);
}

inline void do_sum (mpfi_c_ptr res, unsigned int M, double sigma, unsigned int q, mpfi_t log_q, mpfi_c_t *chis, double t, double t0, bool *co_primes, mpfi_c_ptr sum_err)
{
  unsigned int n;
  mpfi_c_set(res,sum_err);
  for(n=1;n<=M;n++)
    {
      if(!co_primes[n%q])
	continue;
      mpfi_set_ui(n_sigma,n);
      mpfi_log(n_sigma,n_sigma);
      mpfi_neg(n_sigma,n_sigma);
      mpfi_mul_ui(do_sum_t->re,n_sigma,2);
      mpfi_sub(do_sum_t->re,do_sum_t->re,mpfi_log_pi);
      mpfi_add(do_sum_t->re,do_sum_t->re,log_q);
      mpfi_set_d(do_sum_t->im,t0);
      mpfi_add_d(do_sum_t->im,do_sum_t->im,t);
      mpfi_div_ui(do_sum_t->im,do_sum_t->im,2);
      mpfi_mul(do_sum_t->re,do_sum_t->re,do_sum_t->im);
      mpfi_sin(do_sum_t->im,do_sum_t->re);
      mpfi_cos(do_sum_t->re,do_sum_t->re);
      mpfi_mul_d(n_sigma,n_sigma,sigma); 
      mpfi_exp(n_sigma,n_sigma); // n^(-sigma)
      mpfi_c_mul_i(do_sum_t,do_sum_t,n_sigma);
      mpfi_c_mul(do_sum_t,do_sum_t,chis[n%q]);
      mpfi_c_add(res,res,do_sum_t);
    }
}


int main(int argc, char **argv)
{

  int prec,i,j,N,M;
  unsigned int q=5;
  bool *co_primes;
  mpfi_c_t *chis,*fft_vec,*ws,sum_err,sum;
  double B,t0,h,sigma;
  mpfi_t log_q;

  if(argc!=6)
    print_usage();

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    print_usage();
  printf("Running at %d bits of precision.\n",prec);
    
  mpfi_c_setup(prec);

  i=atoi(argv[2]);
  if((i<=0)||(i>32))
    print_usage();
  N=1;
  for(j=0;j<i;j++,N<<=1);
  printf("N set to %d.\n",N);
  B=N*one_over_A;
  printf("B set to %f.\n",B);

  h=calc_h(B);
  printf("h set to %f\n",h);

  sigma=atof(argv[3]);
  if(sigma<=1.0)
    print_usage();
  printf("sigma set to %f\n",sigma);

  t0=atof(argv[4]);
  if(t0-B/2<=0.0)
    print_usage();
  printf("t0 set to %f\n",t0);

  M=atoi(argv[5]);
  if(M<1)
    print_usage();
  printf("M set to %d\n",M);

  mpfi_c_init(sum_err);
  sum_error(sum_err,M,sigma);
  printf("Sum_error(M=%d,sigma=%f)=",M,sigma);mpfi_c_print(sum_err);
  if(!(ws=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*N/2)))
    {
      printf("Failed to allocate memory for ws.Exiting\n");
      exit(1);
    }

  if(!(fft_vec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*N)))
    {
      printf("Failed to allocate memory for fft_vec.Exiting\n");
      exit(1);
    }

  initfft(N,ws);
  for(i=0;i<N;i++)
    mpfi_c_init(fft_vec[i]);

  if(!(co_primes=(bool *) malloc (sizeof(bool)*q)))
    {
      printf("Failed to allocate memory for co_primes. Exiting\n");
      exit(1);
    }
  co_primes[0]=FALSE;
  for(i=1;i<q;i++)
    co_primes[i]=TRUE;

  if(!(chis=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*q)))
    {
      printf("Failed to allocate memory for chis.Exiting\n");
      exit(1);
    }
  for(i=0;i<q;i++)
    mpfi_c_init(chis[i]);
  mpfi_c_set_d(chis[1],1.0,0.0);
  mpfi_c_set_d(chis[2],-1.0,0.0);
  mpfi_c_set_d(chis[3],-1.0,0.0);
  mpfi_c_set_d(chis[4],1.0,0.0);

  mpfi_c_init(sum);
  init_do_sum();
  mpfi_init(log_q);
  mpfi_set_ui(log_q,q);
  mpfi_log(log_q,log_q);
  do_sum(sum,M,sigma,q,log_q,chis,0.0,t0,co_primes,sum_err);
  printf("sum(0.0)=");mpfi_c_print(sum);

}


