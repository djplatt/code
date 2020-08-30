#include <acb.h>

acb_t *fft_init(long n,long prec) {
	long i,halfn,*base;
	acb_t *x,*w;

	base = (long *)malloc(2*sizeof(long)+3*n/2*sizeof(acb_t));
	if (!base) return (acb_t *)0;
	base[0] = n;
	base[1] = prec;
	x = (acb_t *)&base[2];
	w = x+n;
	halfn = n>>1;

	for (i=0;i<n;i++)
		acb_init(x[i]);
	for (i=0;i<halfn;i++) {
		acb_init(w[i]);
		acb_set_ui(w[i],i);
		acb_div_ui(w[i],w[i],halfn,prec);
		acb_exp_pi_i(w[i],w[i],prec);
	}

	return x;
}

void fft_clear(acb_t *x) {
	long *base = (long *)x-2,i;
	for (i=base[0]*3/2-1;i>=0;i--)
		acb_clear(x[i]);
	free(base);
}

void fft(acb_t *x) {
	long *base=(long *)x-2,i,j,k,l;
	long n=base[0],prec=base[1],halfn=n>>1;
	acb_t *p,*w=x+n;
	static acb_t ctemp;
	static int init;

	if (!init) {
		acb_init(ctemp);
		init = 1;
	}

	/* swap each element with one with bit-reversed index */
	for (i=0;i<halfn;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k) j |= 1;
		}
		if (i < j)
			acb_swap(x[i],x[j]);
		else if (i > j)
			acb_swap(x[n-1-i],x[n-1-j]);
		++i, j |= halfn;
		acb_swap(x[i],x[j]);
	}

	for (k=1,l=halfn;k<n;k<<=1,l>>=1)
		for (p=x;p<w;p+=k)
			for (j=0;j<halfn;j+=l,p++) {
				acb_mul(ctemp,p[k],w[j],prec);
				acb_sub(p[k],p[0],ctemp,prec);
				acb_add(p[0],p[0],ctemp,prec);
			}
}

void ifft(acb_t *x) {
	long *base=(long *)x-2,i;
	long n=base[0],prec=base[1],halfn=n>>1;

	fft(x);
	acb_div_ui(x[0],x[0],n,prec);
	acb_div_ui(x[halfn],x[halfn],n,prec);
	for (i=1;i<halfn;i++) {
		acb_div_ui(x[i],x[i],n,prec);
		acb_div_ui(x[n-i],x[n-i],n,prec);
		acb_swap(x[i],x[n-i]);
	}
}
