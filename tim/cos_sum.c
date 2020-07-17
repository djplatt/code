//
// File: cos_sum.c
//
// Code to do Tim Trudgian's cosine sum
//
// DJ Platt 
//
#include "stdlib.h"
#include "stdio.h"
#include "malloc.h"
#include "inttypes.h"
#include "mpfi.h"
#include "mpfi_io.h"

//#define h ((uint64_t) 30)

mpfr_t temp;

write_mpfi(mpfi_ptr x, FILE *outfile)
{
  double xx[2];
  mpfi_get_left(temp,x);
  xx[0]=mpfr_get_d(temp,GMP_RNDD);
  mpfi_get_right(temp,x);
  xx[1]=-mpfr_get_d(temp,GMP_RNDU);
  fwrite(xx,sizeof(double),2,outfile);
  //printf("Output = %f %f\n",xx[0],xx[1]);
}

int main(int argc, char** argv)
{
  if(argc!=5)
    {
      printf("Usage %s <nth file> <outfile> <log2 num files> <h>\n",argv[0]);
      exit(0);
    }

  FILE *outfile=fopen(argv[2],"wb");
  int64_t log2nfiles=atol(argv[3]);
  if((log2nfiles<=0)||(log2nfiles>10))
    {
      printf("Log2 num files must be in [1,10].\n");
      exit(0);
    }
  uint64_t nfiles=1<<log2nfiles;
  uint64_t nth=atol(argv[1]);
  if((nth<0)||(nth>=nfiles))
    {
      printf("n must be in [0,%lu].\n",nfiles);
      exit(0);
    }
  if(!outfile)
    {
      printf("Error opening file %s for binary output.\n",argv[2]);
      exit(0);
    }
  int64_t h=atol(argv[4]);
  if((h<=0)||(h>32))
    {
      printf("h must be in [1,32].\n");
      exit(0);
    }


  mpfr_set_default_prec(200);
  mpfr_init(temp);
  mpfi_t mypi,rpi,rpi_div,cosr,res;
  mpfi_init(mypi);
  mpfi_init(rpi);
  mpfi_init(res);
  mpfi_init(rpi_div);
  mpfi_init(cosr);

  uint64_t r,n;

  mpfi_const_pi(mypi);
  for(r=(1<<(h-log2nfiles))*nth;r<(1<<(h-log2nfiles))*(nth+1);r++)
    {
      mpfi_mul_ui(rpi,mypi,r);
      mpfi_cos(res,rpi);
      for(n=1;n<h;n++)
	{
	  mpfi_div_2exp(rpi_div,rpi,n);
	  mpfi_cos(cosr,rpi_div);
	  mpfi_add(res,res,cosr);
	}
      write_mpfi(res,outfile);
    }
  return(0);
}
