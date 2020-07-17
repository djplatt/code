#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "../includes/mpfi.h"
#include "../includes/mpfi_io.h"
#include "../includes/mpfi_c.h"

int main()
{
  mpfi_c_setup(200);
  mpfi_t x,y;
  mpfi_init(x);mpfi_init(y);
  FILE *outfile=fopen("foo.dat","wb");
  mpfi_set_ui(x,2);
  mpfi_sqrt(x,x);
  mpfi_write_bin(outfile,x);
  fclose(outfile);

  FILE *infile=fopen("foo.dat","rb");
  mpfi_read_bin(infile,y);
  fclose(infile);
  mpfi_print_str("We read in ",y);
  return(0);
}
  
