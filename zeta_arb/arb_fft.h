#ifndef ARB_FFT
#define ARB_FFT
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "acb.h"

#define MAX_SIMPLE_DFT (31)

//acb_t ctemp,htemp;

void initfft(unsigned int MAXN, acb_t *w, int64_t prec) 
{
  int i;
  arb_t temp,temp1;
  acb_t ctemp,htemp;  

  arb_init(temp);
  arb_init(temp1);
  acb_init(ctemp);
  acb_init(htemp);
  acb_init(w[0]);
  arb_set_ui(acb_realref(w[0]),1);
  arb_set_ui(temp,2);
  arb_div_ui(temp,temp,MAXN,prec); // 2/N
  for(i=1;i<MAXN/2;i++)
    {
      arb_mul_ui(temp1,temp,i,prec); // 2*i/N
      acb_init(w[i]);
      arb_sin_cos_pi(acb_imagref(w[i]),acb_realref(w[i]),temp1,prec);
    }
  arb_clear(temp);
  arb_clear(temp1);
  acb_clear(ctemp);
  acb_clear(htemp);
}

/* perform an in place FFT */ 
void fft(acb_t *x,unsigned int n, acb_t *w, int64_t prec) 
{
  static bool init=false;
  static acb_t ctemp;
  if(!init)
    {
      init=true;
      acb_init(ctemp);
    }
  int i,j,k,l;
  acb_t *p,*xend=x+n;

	/* swap each element with one with bit-reversed index */
	for (i=0,l=n>>1;i<l;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k) j |= 1;
		}
		if (i < j)
		  acb_swap(x[i],x[j]);
		else if (i > j)
			acb_swap(x[n-1-i],x[n-1-j]);
		++i, j |= l;
		acb_swap(x[i],x[j]);
	}

	for (k=1,l=n/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<n/2;j+=l,p++) {
			  acb_mul(ctemp,p[k],w[j],prec);
			  acb_sub(p[k],p[0],ctemp,prec);
			  acb_add(p[0],p[0],ctemp,prec);
			}
}

/* perform an in place inverse FFT */
void ifft(acb_t *x,unsigned int n, acb_t *w, int64_t prec) {
	int i,l=n>>1;

	fft(x,n,w,prec);
	acb_div_ui(x[0],x[0],n,prec);
	acb_div_ui(x[l],x[l],n,prec);
	for (i=1;i<l;i++) {
	  acb_div_ui(x[i],x[i],n,prec);
	  acb_div_ui(x[n-i],x[n-i],n,prec);
	  acb_swap(x[i],x[n-i]);
	}
}

/* perform a non-normalised inverse FFT */
void nifft(acb_t *x,unsigned int n, acb_t *w, int64_t prec) {
	int i,l=n>>1;

	fft(x,n,w,prec);
	for (i=1;i<l;i++)
	  acb_swap(x[i],x[n-i]);
}


/* circular convolution; f and g must be distinct */
void convolve(acb_t *result,acb_t *f,acb_t *g,unsigned int n, acb_t *w, int64_t prec) 
{
  int i;
  fft(f,n,w,prec); fft(g,n,w,prec);
  for (i=0;i<n;i++) {
    acb_mul(result[i],f[i],g[i],prec);
    //acb_mul(ctemp,f[i],g[i],prec);
    //acb_set(result[i],ctemp);
  }
  ifft(result,n,w,prec);
}

/* Do iDFT on a Hermitian vector of length 2N, i.e. one that produces real values */
/* N is the reduced length and the ws relate to that N */
/* ws[n]=e(n/N) for n=0..N/2-1 */
/* x is of length 2N with 0..N input values and (N+1..2N-1) anything */
/* omega=exp(2 pi I/2/N) */
/* on exit x[n]->re are set, x[n]->im are garbage. */
void hermidft(acb_t *x, int NN, acb_t *ws, acb_ptr omega, int64_t prec)
{
  static bool init=false;
  static acb_t ctemp;
  if(!init)
    {
      init=true;
      acb_init(ctemp);
    }
  int n,m;
  acb_t *res=&x[NN]; // use back "half" of vector for iDFT

  arb_sub(acb_imagref(res[0]),acb_realref(x[0]),acb_realref(x[NN]),prec);
  arb_add(acb_realref(res[0]),acb_realref(x[0]),acb_realref(x[NN]),prec); // x[NN] now destroyed
  for(n=1;n<NN;n++)
    {
      // compute e(n/(2NN))*i=e((n+NN/2)/(2NN))
      acb_conj(ctemp,x[NN-n]);
      acb_set(res[n],ctemp);
      acb_add(ctemp,ctemp,x[n],prec);
      //printf("ctemp=");mpfi_c_printn(ctemp,10);
      acb_sub(res[n],x[n],res[n],prec);
      //printf("res[%d]=",n);mpfi_c_printn(res[n],10);
      acb_mul_onei(res[n],res[n]);
      m=(n/2)%NN;
      if(m>=NN/2)
	{
	  //printf("multiplying by -");mpfi_c_printn(ws[m-NN/2],10);
	  acb_mul(res[n],res[n],ws[m-NN/2],prec);
	  acb_neg(res[n],res[n]);
	  //printf("res[%d] now=",n);mpfi_c_printn(res[n],10);
	}
      else
	{
	  //printf("multiplying by ");mpfi_c_printn(ws[m],10);
	  acb_mul(res[n],res[n],ws[m],prec);
	  //printf("res[%d] now=",n);mpfi_c_printn(res[n],10);
	}
      if(n&1)
	{
	  //printf("multiplying by omega=");mpfi_c_printn(omega,10);
	  acb_mul(res[n],res[n],omega,prec);
	  //printf("res[%d] now=",n);mpfi_c_printn(res[n],10);
	}
      //printf("ctemp=");mpfi_c_printn(ctemp,10);
      //printf("res[%d]=",n);mpfi_c_printn(res[n],10);
      acb_add(res[n],res[n],ctemp,prec);
    }
  /*
  for(n=0;n<NN;n++)
    {
      printf("res[%d]=",n);mpfi_c_printn(res[n],10);
    }
  */
  fft(res,NN,ws,prec);
  /*
    for(n=0;n<NN;n++)
    {
      printf("res[%d]=",n);mpfi_c_printn(res[n],10);
    }
  */
  for(n=0,m=0;n<NN+NN;n+=2,m++)
    {
      arb_set(acb_realref(x[n]),acb_realref(res[m]));
      arb_set(acb_realref(x[n+1]),acb_imagref(res[m]));
    }
}


#endif
