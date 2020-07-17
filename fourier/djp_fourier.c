#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpfi.h"

#define MAXN    (1<<18)   /* max FFT size was 21 */
#define prec    100       /* bits of working precision */

typedef struct {
	mpfi_t re,im;
} complex;

mpfi_t temp;
complex ctemp;

#define complex_init(x){\
	mpfi_init(x.re); mpfi_init(x.im);}
#define complex_set(x,y){\
	mpfi_set(x.re,y.re); mpfi_set(x.im,y.im);}
#define complex_conj(x,y){\
	mpfi_set(x.re,y.re); mpfi_neg(x.im,y.im);}
#define complex_set_ui(x,n){\
	mpfi_set_ui(x.re,n); mpfi_set_ui(x.im,0);}
#define complex_set_d(x,y){\
        mpfi_set_d(x.re,y);mpfi_set_ui(x.im,0);}
#define complex_add(z,x,y){\
	mpfi_add(z.re,x.re,y.re); mpfi_add(z.im,x.im,y.im);}
#define complex_sub(z,x,y){\
	mpfi_sub(z.re,x.re,y.re); mpfi_sub(z.im,x.im,y.im);}
#define complex_div_ui(x,y,n){\
	mpfi_div_ui(x.re,y.re,n); mpfi_div_ui(x.im,y.im,n);}
#define complex_swap(x,y){\
        mpfi_swap(x.re,y.re); mpfi_swap(x.im,y.im);}
 

/* z = x*y; z must be distinct from x and y! */
#define complex_mul(z,x,y){\
mpfi_mul(z.re,x.re,y.re);\
mpfi_mul(temp,x.im,y.im);\
mpfi_sub(z.re,z.re,temp);\
mpfi_mul(z.im,x.re,y.im);\
mpfi_mul(temp,x.im,y.re);\
mpfi_add(z.im,z.im,temp);}

void cprint(complex x) {
  mpfr_t temp_fr;
  int e1,e;

  mpfr_init(temp_fr);
  mpfi_diam_abs(temp_fr,x.re);
  e1 = mpfr_get_exp(temp_fr);
  mpfr_init(temp_fr);
  mpfi_diam_abs(temp_fr,x.im);
  e = mpfr_get_exp(temp_fr);
  if (e1 > e) e = e1;
  printf("%.14g,%.14g %d\n",
    mpfi_get_d(x.re),mpfi_get_d(x.im),e);
  mpfr_clear(temp_fr);
}

void cwrite(complex c)
{
  mpfi_out_str(stdout,10,0,c.re);
  mpfi_out_str(stdout,10,0,c.im);
};

void mwrite(mpfi_t x,FILE *fp) {
	int e = prec>>2;
	mpfr_out_str(fp,16,e,&x->left,GMP_RNDD);  fprintf(fp,"\n");
	mpfr_out_str(fp,16,e,&x->right,GMP_RNDU); fprintf(fp,"\n");

}

complex w[MAXN>>1];
void initfft() 
{
	int i,k;

	complex_init(w[0]);
	complex_set_ui(w[0],1);
	complex_init(w[1]);
	mpfi_const_pi(temp);
	mpfi_div_ui(temp,temp,MAXN>>1);
	mpfi_cos(w[1].re,temp);
	mpfi_sin(w[1].im,temp);
	for (k=2;k<MAXN/2;k<<=1) {
		complex_init(w[k]);
		complex_mul(w[k],w[k>>1],w[k>>1]);
		for (i=1;i<k;i++) {
			complex_init(w[i+k]);
			complex_mul(w[i+k],w[i],w[k]);
		}
	}
	/*	
	  printf("Outputing various w values.\n");
	  for(i=1;i<MAXN>>1;i<<=1)
	  {
	  printf("%d ",i);
	  cwrite(w[i]);
	  printf("\n");
	  };
	*/

}

/* perform an in place FFT */ 
void fft(complex *x,int n) {
	int i,j,k,l;
	complex *p,*xend=x+n;

	/* swap each element with one with bit-reversed index */
	for (i=0,l=n>>1;i<l;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k) j |= 1;
		}
		if (i < j)
			complex_swap(x[i],x[j])
		else if (i > j)
			complex_swap(x[n-1-i],x[n-1-j]);
		++i, j |= l;
		complex_swap(x[i],x[j]);
	}

	for (k=1,l=MAXN/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<MAXN/2;j+=l,p++) {
				complex_mul(ctemp,p[k],w[j]);
				complex_sub(p[k],p[0],ctemp);
				complex_add(p[0],p[0],ctemp);
			}
}

/* perform an in place inverse FFT */
void ifft(complex *x,int n) {
	int i,l=n>>1;

	fft(x,n);
	complex_div_ui(x[0],x[0],n);
	complex_div_ui(x[l],x[l],n);
	for (i=1;i<l;i++) {
		complex_div_ui(x[i],x[i],n);
		complex_div_ui(x[n-i],x[n-i],n);
		complex_swap(x[i],x[n-i]);
	}
}

/* circular convolution; f and g must be distinct */
void convolve(complex *result,complex *f,complex *g,int n) {
	int i;
	fft(f,n); fft(g,n);
	for (i=0;i<n;i++) {
		complex_mul(ctemp,f[i],g[i]);
		complex_set(result[i],ctemp);
	}
	ifft(result,n);
}


void complex_read(complex *x,FILE *fp) {
	mpfr_inp_str(&x->re->left,fp,16,GMP_RNDD);
	mpfr_inp_str(&x->re->right,fp,16,GMP_RNDU);
	mpfr_inp_str(&x->im->left,fp,16,GMP_RNDD);
	mpfr_inp_str(&x->im->right,fp,16,GMP_RNDU);
}

void mpfi_read(mpfi_t x,FILE *fp) {
	mpfr_inp_str(&x->left,fp,16,GMP_RNDD);
	mpfr_inp_str(&x->right,fp,16,GMP_RNDU);
}

int main(int argc,char **argv) 
{

  complex *fft_vec;
  int i;

  mpfr_set_default_prec(prec);

  mpfi_init(temp);
  complex_init(ctemp);

  if(!(fft_vec=malloc(sizeof(complex)*MAXN)))
    {
      printf("Error allocating memory for fft_vec.\n");
      exit(0);
    };

  for(i=0;i<MAXN;i++)
    {
      complex_init(fft_vec[i]);
      complex_set_ui(fft_vec[i],1);
      complex_div_ui(fft_vec[i],fft_vec[i],(i+1));
    };

  for(i=1;i<MAXN;i=i<<1)
    {printf("%d ",i);cwrite(fft_vec[i]);printf("\n");};


  initfft();

  fft(fft_vec,MAXN);

  for(i=1;i<MAXN;i=i<<1)
    {printf("%d ",i);cwrite(fft_vec[i]);printf("\n");};

  ifft(fft_vec,MAXN);

  for(i=1;i<MAXN;i=i<<1)
    {printf("%d ",i);cwrite(fft_vec[i]);printf("\n");};


  return 0;
}
