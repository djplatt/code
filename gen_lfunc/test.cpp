#include "inttypes.h"
#include "acb_fft.h"

#define N (4)
#define PREC (200)
acb_t x[N],y[N],w[2],res[N];

int main()
{
  for(uint64_t i=0;i<N/2;i++)
    {
      acb_init(w[i]);
      acb_init(x[i+i]);
      acb_init(x[i+i+1]);
      acb_init(y[i+i]);
      acb_init(y[i+i+1]);
      acb_init(res[i+i]);
      acb_init(res[i+i+1]);
    }
  acb_initfft(w,N,PREC);
  for(uint64_t i=0;i<N;i++)
    {
      acb_set_ui(x[i],i);
      acb_set_ui(y[i],10-i);
    }
  acb_convolve(res,x,y,N,w,PREC);

  for(uint64_t i=0;i<N;i++)
    {
      printf("res[%lu]=",i);
      acb_printd(res[i],10);
      printf("\n");
    }
  return(0);
}

