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

#define N ((int) 1<<18)

mpfi_c_t Fft[N],ws_r[N/2],ws_f[N/2];

void init()
{
  int i;
  for(i=0;i<N/2;i++)
    {
      mpfi_c_init(Fft[i]);
      mpfi_c_set_ui(Fft[i],i,i);
      mpfi_c_init(ws_r[i]);
    }
  for(;i<N;i++)
    {
      mpfi_c_init(Fft[i]);
      mpfi_c_set_ui(Fft[i],i,i);
    }
  initfft(N,ws_r);
  for(i=0;i<N/2;i++)
    {
      mpfi_c_init(ws_f[i]);
      mpfi_c_conj(ws_f[i],ws_r[i]);
    }
}

int main(int argc, char **argv)
{
  int prec;
  int count,i,j;

  if(argc!=3)
    exit(0);

  prec=atoi(argv[1]);
  if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    exit(0);
  printf("Running at %d bits of precision.\n",prec);
    
  mpfi_c_setup(prec);

  count=atoi(argv[2]);
  init();

  for(i=0;i<count;i++)
    {
      fft(Fft,N,ws_f);
      fft(Fft,N,ws_r);
      for(j=0;j<N;j++)
	mpfi_c_div_ui(Fft[i],Fft[i],N);
    }
  return(0);
}


