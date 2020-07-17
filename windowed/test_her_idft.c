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

#define N ((int) 8)

mpfi_c_t Fft[N],ws_r[N/2],ws_f[N/2],ws_h[N/4];

void init()
{
  int i;
  for(i=0;i<N/4;i++)
    {
      mpfi_c_init(Fft[i]);
      mpfi_c_init(ws_r[i]);
      mpfi_c_init(ws_h[i]);
    }
  for(;i<N/2;i++)
    {
      mpfi_c_init(Fft[i]);
      mpfi_c_init(ws_r[i]);
    }
  initfft(N,ws_r);
  for(i=0;i<N/2;i++)
    {
      mpfi_c_init(Fft[N-i-1]);
      mpfi_c_init(ws_f[i]);
      mpfi_c_conj(ws_f[i],ws_r[i]);
    }
  for(i=0;i<N/4;i++)
    mpfi_c_set(ws_h[i],ws_r[i+i]);
}

int main(int argc, char **argv)
{
  int prec;
  int count,i,j;
  mpfi_c_t omega;

  if(argc!=3)
    exit(0);

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    exit(0);
  printf("Running at %d bits of precision.\n",prec);
    
  mpfi_c_setup(prec);

  init();

  mpfi_c_init(omega);
  mpfi_div_ui(omega->im,mpfi_2_pi,N);
  mpfi_cos(omega->re,omega->im);
  mpfi_sin(omega->im,omega->im);
  
  for(i=0;i<N;i++)
    {
      mpfi_set_ui(Fft[i]->re,i);
      mpfi_set_ui(Fft[i]->im,0);
    }
  mpfi_set_ui(Fft[0]->re,2);
  mpfi_set_ui(Fft[1]->re,3);
  mpfi_set_ui(Fft[2]->re,5);
  mpfi_set_ui(Fft[3]->re,7);
  mpfi_set_ui(Fft[4]->re,11);
  mpfi_set_ui(Fft[5]->re,13);
  mpfi_set_ui(Fft[6]->re,17);
  mpfi_set_ui(Fft[7]->re,19);
  for(i=0;i<N;i++)
    mpfi_c_div_ui(Fft[i],Fft[i],N);
  fft(Fft,N,ws_f);
  for(i=0;i<N;i++)
    {
      printf("Fft[%d]=",i);mpfi_c_printn(Fft[i],10);
    }
  hermidft(Fft,N/2,ws_h,omega);
  for(i=0;i<N;i++)
    {
      printf("Fft[%d]=",i);mpfi_c_printn(Fft[i],10);
    }
  return(0);
}


