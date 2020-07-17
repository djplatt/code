#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "../includes/mpfi.h"
#include "../includes/mpfi_io.h"

#define one_over_twoB ((double) 5.0/64.0)

int main(int argc, char** argv)
{
  if(argc!=4)
    {
      printf("Usage:- q_it <q> <t0> <outfile>.\n");
      exit(0);
    }
  unsigned long int q=atol(argv[1]);
  double t0=atof(argv[2]);
  FILE *outfile=fopen(argv[3],"wb");
  if(!outfile)
    {
      printf("Failed to open file %s for binary output. Exiting.\n",argv[3]);
      exit(0);
    }

  mpfr_set_default_prec(100);
  mpfi_t t,logq,s,c,tlogq;
  mpfr_t m;
  long unsigned int step;
  double x[4]={1.0,1.0,0.0,0.0}; // [1,1]+[0,0]i
  mpfi_init(t);
  mpfi_init(s);
  mpfi_init(c);
  mpfi_init(tlogq);
  mpfi_init(logq);
  mpfr_init(m);
  mpfi_set_ui(t,q);
  mpfi_log(logq,t);
  mpfi_set_d(t,one_over_twoB);
  //printf("writing %g %g %g %g\n",x[0],x[1],x[2],x[3]);
  fwrite(x,sizeof(double),4,outfile);
  for(step=1;step<=t0/one_over_twoB;step++)
    {
      mpfi_mul(tlogq,logq,t);
      mpfi_sin(s,tlogq);
      mpfi_cos(c,tlogq);
      mpfi_get_left(m,c);
      x[0]=mpfr_get_d(m,GMP_RNDD);
      mpfi_get_right(m,c);
      x[1]=mpfr_get_d(m,GMP_RNDU);
      mpfi_get_left(m,s);
      x[2]=mpfr_get_d(m,GMP_RNDD);
      mpfi_get_right(m,s);
      x[3]=mpfr_get_d(m,GMP_RNDU);
      printf("writing %20.18e %20.18e %20.18e %20.18e\n",x[0],x[1],x[2],x[3]);
      fwrite(x,sizeof(double),4,outfile);
      mpfi_add_d(t,t,one_over_twoB);
    }
  printf("\n\n");

  unsigned short int old_cw,new_cw;
  unsigned int old__SSE_cw,new__SSE_cw;
#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  :: "m" (cw))
#define __SSE_getcw(cw) asm ("stmxcsr %0" : "=m" (cw))
#define __SSE_setcw(cw) asm ("ldmxcsr %0" :: "m" (cw))

  // must leave fpu rounding to crlibm
  //_fpu_getcw(old_cw);
  //new_cw=old_cw&0x3FF;
  //new_cw=new_cw|0x400;
  //_fpu_setcw(new_cw);

  __SSE_getcw(old__SSE_cw);
  new__SSE_cw=old__SSE_cw&0x1FBF; // zeros FTZ(15),RC(14,13) and DAZ(6) bits
  new__SSE_cw=new__SSE_cw|0x2000; // sets RC to Round Down (14=0, 13=1)
  __SSE_setcw(new__SSE_cw);

  mpfi_set_d(t,one_over_twoB);
  //printf("writing %g %g %g %g\n",x[0],x[1],x[2],x[3]);
  fwrite(x,sizeof(double),4,outfile);
  for(step=1;step<=t0/one_over_twoB;step++)
    {
      mpfi_mul(tlogq,logq,t);
      mpfi_sin(s,tlogq);
      mpfi_cos(c,tlogq);
      mpfi_get_left(m,c);
      x[0]=mpfr_get_d(m,GMP_RNDD);
      mpfi_get_right(m,c);
      x[1]=mpfr_get_d(m,GMP_RNDU);
      mpfi_get_left(m,s);
      x[2]=mpfr_get_d(m,GMP_RNDD);
      mpfi_get_right(m,s);
      x[3]=mpfr_get_d(m,GMP_RNDU);
      printf("writing %20.18e %20.18e %20.18e %20.18e\n",x[0],x[1],x[2],x[3]);
      fwrite(x,sizeof(double),4,outfile);
      mpfi_add_d(t,t,one_over_twoB);
    }



  fclose(outfile);
  return(0);
}
