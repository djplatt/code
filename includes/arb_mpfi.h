#include "inttypes.h"

void mpfi_arb_set(mpfi_ptr x, arb_ptr y)
{
  static bool init=false;
  static mpfr_t l,r;
  if(!init)
    {
      init=true;
      mpfr_init(l);
      mpfr_init(r);
    }
  arb_get_interval_mpfr(l,r,y);
  mpfi_interv_fr(x,l,r);
}

mpfi_c_t *w;
mpfi_c_t *X;

void init_arb_mpfi_fft(uint64_t N, int64_t prec)
{
  printf("Initialising MPFI.\n");
  mpfi_c_setup(prec);
  X=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*N);
  for(uint64_t n=0;n<N;n++)
    mpfi_c_init(X[n]);
  w=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*N/2);
  initfft(N,w);
}

void arb_mpfi_fft(acb_t *res, uint64_t N, acb_t *ww, int64_t prec)
{
  if(N!=_MPFI_FFT_LEN)
    {
      printf("Length mismatch between ARB and MPFI FFT's %lu and %lu resp. Exiting.\n",N,_MPFI_FFT_LEN);
      exit(0);
    }
  static bool init=false;
  static mpfr_t fft_l,fft_r;
  if(!init)
    {
      init=true;
      mpfr_init(fft_l);
      mpfr_init(fft_r);
    }

  for(uint64_t n=0;n<N;n++)
    {
      arb_get_interval_mpfr(fft_l,fft_r,acb_realref(res[n]));
      mpfi_interv_fr(X[n]->re,fft_l,fft_r);
      arb_get_interval_mpfr(fft_l,fft_r,acb_imagref(res[n]));
      mpfi_interv_fr(X[n]->im,fft_l,fft_r);
    }
  fft(X,N,w);
  for(uint64_t n=0;n<N;n++)
    {
      mpfi_get_left(fft_l,X[n]->re);
      mpfi_get_right(fft_r,X[n]->re);
      arb_set_interval_mpfr(acb_realref(res[n]),fft_l,fft_r,prec);
      mpfi_get_left(fft_l,X[n]->im);
      mpfi_get_right(fft_r,X[n]->im);
      arb_set_interval_mpfr(acb_imagref(res[n]),fft_l,fft_r,prec);
    }
}

void arb_mpfi_ifft(acb_t *res, uint64_t N, acb_t *ww, int64_t prec)
{
  if(N!=_MPFI_FFT_LEN)
    {
      printf("Length mismatch between ARB and MPFI FFT's %lu and %lu resp. Exiting.\n",N,_MPFI_FFT_LEN);
      exit(0);
    }
  static bool init=false;
  static mpfr_t fft_l,fft_r;
  if(!init)
    {
      init=true;
      mpfr_init(fft_l);
      mpfr_init(fft_r);
    }

  for(uint64_t n=0;n<N;n++)
    {
      arb_get_interval_mpfr(fft_l,fft_r,acb_realref(res[n]));
      mpfi_interv_fr(X[n]->re,fft_l,fft_r);
      arb_get_interval_mpfr(fft_l,fft_r,acb_imagref(res[n]));
      mpfi_interv_fr(X[n]->im,fft_l,fft_r);
    }
  ifft(X,N,w);
  for(uint64_t n=0;n<N;n++)
    {
      mpfi_get_left(fft_l,X[n]->re);
      mpfi_get_right(fft_r,X[n]->re);
      arb_set_interval_mpfr(acb_realref(res[n]),fft_l,fft_r,prec);
      mpfi_get_left(fft_l,X[n]->im);
      mpfi_get_right(fft_r,X[n]->im);
      arb_set_interval_mpfr(acb_imagref(res[n]),fft_l,fft_r,prec);
    }
}
