// gcc bench.c -lgmp -lmpfr -lmpfi
#include"cycle.h"
#include<stdio.h>
#include<gmp.h>
#include<mpfr.h>
#include<mpfi.h>
#include<mpfi_io.h>


int main(/*int argc;char *argv[]*/){
  int i;ticks x,y;
  int p=300; /* 3328 bits of precision == 1066 digits precision */
  mpf_t g,b,a;
  mpfr_prec_t q=p; /* digits of precision */
  mpfr_t h,d,c;
  mp_prec_t r=p; /* digits of precision */
  mpfi_t j,e,f;
  double x1=1.0,x2=1.00000001;

  x=getticks();for(i=0;i<10000;i++){ x1*=x2;}y=getticks();
  printf("double time= % 17.15g\n",elapsed(y,x));printf("%e\n",x1);
  
  
  mpf_init2(g,p);mpf_init2(b,p);mpf_init2(a,p);
  mpf_set_d(b,2.0);mpf_set_d(a,2.0);
  x=getticks();for(i=0;i<10000;i++){ mpf_mul(g,b,a);}y=getticks();
  printf("mpf time= % 20.15g\n",elapsed(y,x));
  mpf_clear(g);mpf_clear(b);mpf_clear(a);

  mpfr_init2(h,q);mpfr_init2(d,q);mpfr_init2(c,q);
  mpfr_set_d(d,2.0,GMP_RNDZ);mpfr_set_d(c,2.0,GMP_RNDD);
  x=getticks();for(i=0;i<10000;i++){ mpfr_mul(h,d,c,GMP_RNDD);} y=getticks();
  printf("mpfr time=% 20.15g rounding towards zero\n",elapsed(y,x));
  mpfr_clear(h);mpfr_clear(d);mpfr_clear(c);

  mpfr_init2(h,q);mpfr_init2(d,q);mpfr_init2(c,q);
  mpfr_set_d(d,2.0,GMP_RNDU);mpfr_set_d(c,2.0,GMP_RNDD);
  x=getticks();for(i=0;i<10000;i++){ mpfr_mul(h,d,c,i%2?GMP_RNDU:GMP_RNDD);}
  y=getticks();printf("mpfr time=% 20.15g alternate rounding\n",elapsed(y,x));
  mpfr_clear(h);mpfr_clear(d);mpfr_clear(c);

  mpfi_init2(j,r);mpfi_init2(f,r);mpfi_init2(e,r);
  mpfi_set_d(f,2.0);mpfi_set_d(e,2.0);mpfi_div_ui(f,f,3);mpfi_div_ui(e,e,3);
  x=getticks();for(i=0;i<10000/2;i++){ mpfi_mul(j,f,e);}y=getticks();
  printf("mpfi time=% 20.15g\n",elapsed(y,x));
  mpfi_clear(j);mpfi_clear(f);mpfi_clear(e);

  printf("But beware of the following\n");
  mpf_init2(g,p);mpf_init2(b,p);mpf_init2(a,p);
  mpf_set_d(b,2.0);mpf_set_d(a,2.0);
  x=getticks();for(i=0;i<10000;i++){ mpf_mul(a,b,a);}y=getticks();
  printf("mpf time= % 20.15g using inline arguments\n",elapsed(y,x));
  mpf_clear(g);mpf_clear(b);mpf_clear(a);

  return(0);
}
