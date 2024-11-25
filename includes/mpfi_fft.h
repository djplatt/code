#include <stdio.h>
#include <string.h>
#include <math.h>
//#include "mpfi.h"

mpfi_c_t ctemp,htemp;

void initfft(unsigned int MAXN, mpfi_c_t *w) 
{
	int i;
	mpfi_t temp,temp1;

	mpfi_init(temp);
	mpfi_init(temp1);
	mpfi_c_init(ctemp);
	mpfi_c_init(htemp);
	mpfi_c_init(w[0]);
	mpfi_c_set_ui(w[0],1,0);
	mpfi_const_pi(temp);
	mpfi_div_ui(temp,temp,MAXN>>1); // 2Pi/N
	for(i=1;i<MAXN/2;i++)
	  {
	    mpfi_mul_ui(temp1,temp,i);
	    mpfi_c_init(w[i]);
	    mpfi_cos(w[i]->re,temp1);
	    mpfi_sin(w[i]->im,temp1);
	  }
	/*
	mpfi_div_ui(temp,temp,MAXN>>1);
	mpfi_cos(w[1]->re,temp);
	mpfi_sin(w[1]->im,temp);
	for (k=2;k<MAXN/2;k<<=1) {
		mpfi_c_init(w[k]);
		mpfi_c_mul(w[k],w[k>>1],w[k>>1]);
		for (i=1;i<k;i++) {
			mpfi_c_init(w[i+k]);
			mpfi_c_mul(w[i+k],w[i],w[k]);
		}
	}
	*/
	mpfi_clear(temp);
	mpfi_clear(temp1);
}

/* perform an in place FFT */ 
void fft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w) {
	int i,j,k,l;
	mpfi_c_t *p,*xend=x+n;

	/* swap each element with one with bit-reversed index */
	for (i=0,l=n>>1;i<l;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k) j |= 1;
		}
		if (i < j)
		  mpfi_c_swap(x[i],x[j]);
		else if (i > j)
			mpfi_c_swap(x[n-1-i],x[n-1-j]);
		++i, j |= l;
		mpfi_c_swap(x[i],x[j]);
	}

	for (k=1,l=n/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<n/2;j+=l,p++) {
				mpfi_c_mul(ctemp,p[k],w[j]);
				mpfi_c_sub(p[k],p[0],ctemp);
				mpfi_c_add(p[0],p[0],ctemp);
			}
}

/* perform an in place inverse FFT */
void ifft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w) {
	int i,l=n>>1;

	fft(x,n,w);
	mpfi_c_div_ui(x[0],x[0],n);
	mpfi_c_div_ui(x[l],x[l],n);
	for (i=1;i<l;i++) {
		mpfi_c_div_ui(x[i],x[i],n);
		mpfi_c_div_ui(x[n-i],x[n-i],n);
		mpfi_c_swap(x[i],x[n-i]);
	}
}

/* perform a non-normalised inverse FFT */
void nifft(mpfi_c_t *x,unsigned int n, mpfi_c_t *w) {
	int i,l=n>>1;

	fft(x,n,w);
	for (i=1;i<l;i++)
	  mpfi_c_swap(x[i],x[n-i]);
}


/* circular convolution; f and g must be distinct */
void convolve(mpfi_c_t *result,mpfi_c_t *f,mpfi_c_t *g,unsigned int n, mpfi_c_t *w) {
	int i;
	fft(f,n,w); fft(g,n,w);
	for (i=0;i<n;i++) {
		mpfi_c_mul(ctemp,f[i],g[i]);
		mpfi_c_set(result[i],ctemp);
	}
	ifft(result,n,w);
}

/* Do iDFT on a Hermitian vector of length 2N, i.e. one that produces real values */
/* N is the reduced length and the ws relate to that N */
/* ws[n]=e(n/N) for n=0..N/2-1 */
/* x is of length 2N with 0..N input values and (N+1..2N-1) anything */
/* omega=exp(2 pi I/2/N) */
/* on exit x[n]->re are set, x[n]->im are garbage. */
void hermidft(mpfi_c_t *x, int N, mpfi_c_t *ws, mpfi_c_ptr omega)
{
  int n,m;
  mpfi_c_t *res=&x[N]; // use back "half" of vector for iDFT

  mpfi_sub(res[0]->im,x[0]->re,x[N]->re);
  mpfi_add(res[0]->re,x[0]->re,x[N]->re); // x[N] now destroyed
  for(n=1;n<N;n++)
    {
      // compute e(n/(2N))*i=e((n+N/2)/(2N))
      mpfi_c_conj(ctemp,x[N-n]);
      mpfi_c_set(res[n],ctemp);
      mpfi_c_add(ctemp,ctemp,x[n]);
      //printf("ctemp=");mpfi_c_printn(ctemp,10);
      mpfi_c_sub(res[n],x[n],res[n]);
      //printf("res[%d]=",n);mpfi_c_printn(res[n],10);
      mpfi_c_muli(res[n]);
      m=(n/2)%N;
      if(m>=N/2)
	{
	  //printf("multiplying by -");mpfi_c_printn(ws[m-N/2],10);
	  mpfi_c_mul(res[n],res[n],ws[m-N/2]);
	  mpfi_c_neg(res[n],res[n]);
	  //printf("res[%d] now=",n);mpfi_c_printn(res[n],10);
	}
      else
	{
	  //printf("multiplying by ");mpfi_c_printn(ws[m],10);
	  mpfi_c_mul(res[n],res[n],ws[m]);
	  //printf("res[%d] now=",n);mpfi_c_printn(res[n],10);
	}
      if(n&1)
	{
	  //printf("multiplying by omega=");mpfi_c_printn(omega,10);
	  mpfi_c_mul(res[n],res[n],omega);
	  //printf("res[%d] now=",n);mpfi_c_printn(res[n],10);
	}
      //printf("ctemp=");mpfi_c_printn(ctemp,10);
      //printf("res[%d]=",n);mpfi_c_printn(res[n],10);
      mpfi_c_add(res[n],res[n],ctemp);
    }
  /*
  for(n=0;n<N;n++)
    {
      printf("res[%d]=",n);mpfi_c_printn(res[n],10);
    }
  */
  fft(res,N,ws);
  /*
    for(n=0;n<N;n++)
    {
      printf("res[%d]=",n);mpfi_c_printn(res[n],10);
    }
  */
  for(n=0,m=0;n<N+N;n+=2,m++)
    {
      mpfi_set(x[n]->re,res[m]->re);
      mpfi_set(x[n+1]->re,res[m]->im);
    }
}
