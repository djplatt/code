#include "inttypes.h"
#include "../includes/int_double12.0.h"
#include "../includes/int-fft.h"


int main()
{
  _fpu_rndd();

  unsigned int conv_size[1];

  int_complex x[6];
  for(uint64_t i=0;i<6;i++)
    {
      x[i].real=i;
      x[i].imag=0;
    }
  init_ws(_ws);
  init_bluestein_fft(6,conv_size,bs,b_star_conjs);

  bluestein_fft(1,x,6,a,bs,b_star_conjs,conv_size[0],b_spare);

  for(uint64_t i=0;i<6;i++)
    print_int_complex_str("Res: ",x[i]);

  return(0);
}
