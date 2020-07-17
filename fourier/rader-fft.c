#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <malloc.h>
#include "mpfi.h"

#define MAXN    (16)   /* max FFT size was 21 */
#define prec    53       /* bits of working precision */

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

complex w[16>>1];
void initfft() 
{
	int i,k;

	complex_init(w[0]);
	complex_set_ui(w[0],1);
	complex_init(w[1]);
	mpfi_const_pi(temp);
	mpfi_div_ui(temp,temp,16>>1);
	mpfi_cos(w[1].re,temp);
	mpfi_sin(w[1].im,temp);
	for (k=2;k<16/2;k<<=1) {
		complex_init(w[k]);
		complex_mul(w[k],w[k>>1],w[k>>1]);
		for (i=1;i<k;i++) {
			complex_init(w[i+k]);
			complex_mul(w[i+k],w[i],w[k]);
		}
	}
	/*	
	  printf("Outputing various w values.\n");
	  for(i=1;i<16>>1;i<<=1)
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

	for (k=1,l=16/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<16/2;j+=l,p++) {
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

void rader_fft (unsigned int p, unsigned int g, complex *ivec, complex *ovec)
{
  unsigned int *g_q,q,conv_size;
  complex *a,*b,*c,old_x_0;
  mpfi_t two_pi_n,ctemp;

  complex_init(old_x_0);
  complex_set(old_x_0,ivec[0]);

  printf("Input\n");
  for(q=0;q<p;q++)
    {cwrite(ivec[q]);printf("\n");};
  printf("Output\n");
  for(q=0;q<p;q++)
    {cwrite(ovec[q]);printf("\n");};


  if(!(g_q=malloc(sizeof(unsigned int)*(p-1))))
    {
      printf("Error allocating memory for g_q. Exiting.\n");
      exit(0);
    };
  g_q[0]=1;
  g_q[1]=g;
  for(q=2;q<(p-1);q++)
    g_q[q]=(g_q[q-1]*g)%p;

  q=2*(p-1)-1;
  conv_size=1;
  while(conv_size<q)
    conv_size<<=1;
  printf("Using convolution size %d.\n",conv_size);

  if(!(a=malloc(sizeof(complex)*conv_size)))
  {
    printf("Error allocating memory for a. Exiting.\n");
    exit(0);
  };
  if(!(b=malloc(sizeof(complex)*conv_size)))
  {
    printf("Error allocating memory for b. Exiting.\n");
    exit(0);
  };
  for(q=0;q<conv_size;q++)
    {
      complex_init(a[q]);
      complex_init(b[q]);
    };

  mpfi_init(two_pi_n);
  mpfi_init(ctemp);
  mpfi_const_pi(two_pi_n);
  mpfi_mul_ui(two_pi_n,two_pi_n,2);
  mpfi_div_ui(two_pi_n,two_pi_n,p);

  a[0]=ivec[1];
  mpfi_cos(b[0].re,two_pi_n);
  mpfi_sin(b[0].im,two_pi_n);
  mpfi_neg(b[0].im,b[0].im);

  for(q=1;q<p-1;q++)
    {
      a[q]=ivec[g_q[q]];
      mpfi_mul_ui(ctemp,two_pi_n,g_q[(p-q-1)]);
      mpfi_cos(b[q].re,ctemp);
      mpfi_sin(b[q].im,ctemp);
      mpfi_neg(b[q].im,b[q].im);
    };

  for(q=p-1;q<conv_size;q++)
    {
      mpfi_set_ui(a[q].re,0);
      mpfi_set_ui(a[q].im,0);
      mpfi_set_ui(b[q].re,0);
      mpfi_set_ui(b[q].im,0);
    };

  for(q=0;q<conv_size;q++)
    {
      printf("a[%d]: ",q);cwrite(a[q]);
      printf("\nb[%d]: ",q);cwrite(b[q]);
      printf("\n");
    };

  convolve(a,a,b,conv_size);

  for(q=0;q<conv_size;q++)
    {
      printf("a[%d]: ",q);
      cwrite(a[q]);
      printf("\n");
    };

  complex_set(ovec[0],ivec[0]);
  printf("ovec[0] = ");cwrite(ovec[0]);
  for(q=1;q<p;q++)
    {
      complex_add(ovec[0],ovec[0],ivec[q]);
      printf("\novec[0] = ");cwrite(ovec[0]);
    };
  printf("\n");

  for(q=1;q<p;q++)
    complex_add(ovec[q],old_x_0,a[q-1]);
};





int main(int argc,char **argv) 
{

  complex *ifft_vec,*offt_vec;
  int i;

  mpfr_set_default_prec(prec);

  mpfi_init(temp);
  complex_init(ctemp);

  if(!(ifft_vec=malloc(sizeof(complex)*7)))
    {
      printf("Error allocating memory for fft_vec.\n");
      exit(0);
    };
  if(!(offt_vec=malloc(sizeof(complex)*7)))
    {
      printf("Error allocating memory for fft_vec.\n");
      exit(0);
    };

  for(i=0;i<7;i++)
    {
      complex_init(ifft_vec[i]);
      complex_init(offt_vec[i]);
      complex_set_ui(ifft_vec[i],i+1);
      //      complex_div_ui(fft_vec[i],fft_vec[i],(i+1));
    };

  //  for(i=0;i<7;i++)
  //{printf("%d ",i);cwrite(fft_vec[i]);printf("\n");};


    initfft();

  //  fft(fft_vec,MAXN);

  rader_fft(7,3,ifft_vec,offt_vec);

  for(i=0;i<7;i++)
    {printf("%d ",i);cwrite(offt_vec[i]);printf("\n");};

  //  ifft(fft_vec,MAXN);

  //for(i=0;i<7;i++)
  //{printf("%d ",i);cwrite(fft_vec[i]);printf("\n");};


  return 0;
}
