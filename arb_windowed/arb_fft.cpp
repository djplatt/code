#ifndef ARB_FFT
#define ARB_FFT
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "acb.h"

#define MAX_SIMPLE_DFT (31)

acb_t ctemp,htemp;

void initfft(unsigned int MAXN, acb_t *w, int64_t prec) 
{
	int i;
	acb_t temp,temp1;

	arb_init(temp);
	arb_init(temp1);
	acb_init(ctemp);
	acb_init(htemp);
	acb_init(w[0]);
	arb_set_ui(acb_realref(w[0]),1,0);
	arb_set_ui(temp,2);
	arb_div_ui(temp,temp,MAXN,prec); // 2/N
	for(i=1;i<MAXN/2;i++)
	  {
	    arb_mul_ui(temp1,temp,i,prec); // 2*i/N
	    acb_init(w[i]);
	    acb_sin_cos_pi(acb_imagref(w[i]),acb_realref(w[i]),temp1,prec);
	  }
	arb_clear(temp);
	arb_clear(temp1);
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

mpfi_c_t *simple_ws[MAX_SIMPLE_DFT];
mpfi_c_t dft_spare[MAX_SIMPLE_DFT+1];


void simple_dft_init()
{
  unsigned long int n;
  for(n=0;n<=MAX_SIMPLE_DFT;n++)
    mpfi_c_init(dft_spare[n]);
}

// do a simple dft by O(n^2) method
// input in x[0],x[stride].., output to x 
void simple_dft(mpfi_c_t *x, unsigned long int N, unsigned long int stride)
{
  long unsigned int n,m,ptr;
  //printf("In simple_dft.\n");
  if(!simple_ws[N-1]) // never seen before
    {
      simple_ws[N-1]=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*N);
      if(!simple_ws[N-1])
	fatal_error("Failed to allocate memory for simple_ws.");
      mpfi_c_init(simple_ws[N-1][0]);
      mpfi_c_one(simple_ws[N-1][0]);
      //printf("Done [0].\n");
      mpfi_const_pi(dft_spare[1]->re);
      //mpfi_neg(dft_spare[1]->re,dft_spare[1]->re);
      mpfi_mul_ui(dft_spare[1]->im,dft_spare[1]->re,2); // -2 Pi
      for(n=1;n<N;n++)
	{
	  //printf("Doing %lu\n",n);
	  mpfi_set_ui(dft_spare[0]->re,n);
	  mpfi_div_ui(dft_spare[0]->im,dft_spare[0]->re,N);
	  mpfi_mul(dft_spare[0]->re,dft_spare[0]->im,dft_spare[1]->im); // -2pi n/N
	  mpfi_c_init(simple_ws[N-1][n]);
	  mpfi_cos(simple_ws[N-1][n]->re,dft_spare[0]->re);
	  mpfi_sin(simple_ws[N-1][n]->im,dft_spare[0]->re);
	}
    }
  //printf("ws made.\n");
  /*
  for(n=0;n<N;n++)
    {
      printf("ws[%lu]=",n);
      mpfi_c_print(simple_ws[N-1][n]);
    }
  */
  mpfi_c_set(dft_spare[0],x[0]);
  for(m=1;m<N;m++)
    mpfi_c_add(dft_spare[0],dft_spare[0],x[m*stride]);
  for(n=1;n<N;n++)
    {
      mpfi_c_set(dft_spare[n],x[0]); // x[0]*e(0)
      for(m=1;m<N;m++)
	{
	  ptr=(m*n)%N;
	  mpfi_c_mul(dft_spare[MAX_SIMPLE_DFT],x[m*stride],simple_ws[N-1][ptr]); // x[m]+e(-mn/N)
	  //printf("term for n=%lu m=%lu is",n,m);mpfi_c_print(dft_spare[MAX_SIMPLE_DFT]);
	  mpfi_c_add(dft_spare[n],dft_spare[n],dft_spare[MAX_SIMPLE_DFT]);
	}
    }
  for(n=0;n<N;n++)
    mpfi_c_set(x[n*stride],dft_spare[n]);
}

void dft(mpfi_c_t *x, unsigned long int len, unsigned long int stride)
{
  //printf("In dft with len=%lu and stride=%lu\n",len,stride);
  if(len==2)
    {
      mpfi_c_add(dft_spare[0],x[0],x[stride]);
      mpfi_c_sub(dft_spare[1],x[0],x[stride]);
      mpfi_c_set(x[0],dft_spare[0]);
      mpfi_c_set(x[stride],dft_spare[1]);
      return;
    }
  if(len<=MAX_SIMPLE_DFT)
    {
      simple_dft(x,len,stride);
      return;
    }
  fatal_error("Can't do arbitrary lengths yet.");
}

void ndft(mpfi_c_t *x, unsigned long int N, unsigned long int num_ffts, unsigned long int *lengths)
{
  unsigned long int i,k,l,stride=N,depth;
  //printf("ndft called with #ffts=%lu\n",num_ffts);
  //for(i=0;i<num_ffts;i++) printf("length[%lu]=%lu.\n",i,lengths[i]);

  for(i=0;i<num_ffts;i++)
    {
      depth=N/stride;
      stride/=lengths[i];
      for(k=0;k<depth*lengths[i]*stride;k+=lengths[i]*stride)
	for(l=0;l<stride;l++)
	  {
	    //printf("Doing length %lu dft starting at %lu with stride %lu.\n",lengths[i],k+l,stride);
	    dft(x+k+l,lengths[i],stride);
	  }
    }
}

inline int co_prime(long unsigned int, long unsigned int);

void prep_omegas_nd(long unsigned int *offsets, long unsigned int q, mpfi_c_t *omegas,long unsigned int n_dims,long unsigned int *dims,
					long unsigned int phi_q)
{
  long unsigned int n,n1;
  mpfi_const_pi(dft_spare[0]->re);
  //printf("pi=");mpfi_print(dft_spare[0]->re);
  mpfi_mul_ui(dft_spare[0]->re,dft_spare[0]->re,2);
  mpfi_div_ui(dft_spare[0]->im,dft_spare[0]->re,q); // 2 pi/q
  
  n1=0;
  for(n=1;n<q;n++)
    if(co_prime(n,q))
      {
	mpfi_mul_ui(dft_spare[0]->re,dft_spare[0]->im,n); // 2 pi n/q
	//printf("2pi %lu/%lu=",n,q);mpfi_print(dft_spare[0]->re);
	mpfi_sin(omegas[offsets[n1]]->im,dft_spare[0]->re);
	mpfi_cos(omegas[offsets[n1]]->re,dft_spare[0]->re);
	n1++;
      }
  /*
  for(n=0;n<phi_q;n++)
    {printf("pre fft ");mpfi_c_print(omegas[n]);}
  */
  //printf("calling ndft.\n");
  ndft(omegas,phi_q,n_dims,dims);
  
  //printf("ndft returned.\n");
  /*
  for(n=0;n<phi_q;n++)
    {printf("pre fft ");mpfi_c_print(omegas[n]);}
  exit(0);
  */
  mpfi_set_ui(dft_spare[0]->re,1);
  mpfi_set_ui(dft_spare[0]->im,q);
  mpfi_sqrt(dft_spare[1]->re,dft_spare[0]->im);
  mpfi_div(dft_spare[1]->im,dft_spare[0]->re,dft_spare[1]->re); // 1/sqrt(q)

  for(n=0;n<phi_q;n++)
    {
      mpfi_mul(omegas[n]->re,omegas[n]->re,dft_spare[1]->im);
      mpfi_mul(omegas[n]->im,omegas[n]->im,dft_spare[1]->im);
    }
}

void finish_omega(mpfi_c_ptr omega, int neg_one) // do the i^(-delta) and sqrt'ing
{
  
  if(neg_one) // mult by -i and conj
    {
      mpfi_swap(omega->re,omega->im);
    }
  else // conj
    mpfi_neg(omega->im,omega->im);
  
  mpfi_c_sqrt(omega,omega);
  if(mpfi_is_neg(omega->re))
    {
      printf("omega in left half plane. negating.\n");
      mpfi_c_neg(omega,omega);
    }
}



#endif
