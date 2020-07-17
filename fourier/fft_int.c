#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpfi.h"

#define MAXN    (1<<21)   /* max FFT size */
#define prec    100       /* bits of working precision */
#define logZ32  0.96026   /* > log(zeta(3/2)) */
#define error_terms 1

int N,r0,r1,r,eta,samples,ntaylor,q,skip;
double nmax;
#define range (nmax+1)
mpfi_t residue;
double C,alpha;

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

void mwrite(mpfi_t x,FILE *fp) {
	int e = prec>>2;
	mpfr_out_str(fp,16,e,&x->left,GMP_RNDD);  fprintf(fp,"\n");
	mpfr_out_str(fp,16,e,&x->right,GMP_RNDU); fprintf(fp,"\n");

}

complex w[MAXN/2];
void initfft(void) {
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
}

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

complex *ftemp,*gtemp,*L;
mpfi_t ferror_term,gerror_term;

void compute_samples(void) {
	int i,k;
	complex *f=L,*g=&L[ntaylor*samples];

	for (k=0;k<ntaylor;k++,f+=samples,g+=samples) {
		/* ftemp = f reversed */
		complex_set_ui(ftemp[samples],0);
		complex_set(ftemp[0],f[0]);
		for (i=1;i<samples;i++) {
			complex_set_ui(ftemp[i],0);
			complex_set(ftemp[2*samples-i],f[i]);
		}

		/* gtemp = g */
		for (i=0;i<samples;i++) {
			complex_set(gtemp[i],g[i]);
			complex_set_ui(gtemp[i+samples],0);
		}

		convolve(ftemp,ftemp,gtemp,2*samples);

		/* f = ftemp */
		for (i=0;i<samples;i++)
			complex_set(f[i],ftemp[i]);
	}
	for (i=0,f=L;i<samples;i++,f++)
		complex_set(ftemp[i],f[0]);
	for (k=1;k<ntaylor;k++)
		for (i=0;i<samples;i++,f++)
			complex_add(ftemp[i],ftemp[i],f[0]);

	/* add convolution error term */
#ifdef error_terms
	mpfi_mul(temp,ferror_term,gerror_term);
	for (i=0;i<samples;i++) {
		mpfi_add(ftemp[0].re,ftemp[0].re,temp);
		mpfi_add(ftemp[0].im,ftemp[0].im,temp);
	}
#endif
}

mpfi_t epsilon,Fhatfac,error_term;
double epsd,delta,mu,logK,cprime;

void init_truncation_error(void) {
	delta = 0.5*M_PI*(1.0-eta*0.001);
	mu = (2.0*(r1+1)-r)/(2.0*r);
	logK = 0.5*((r+3)*M_LN2-log(r*delta)+(r-1)*delta)+log(C);
	cprime = 0.5*r*(mu+alpha+0.5);
	if (cprime < 0.0) cprime = 0.0;

	/* Fhatfac = 1/(1-exp(-Pi*A)) */
	mpfi_init(Fhatfac);       mpfi_mul_ui(temp,epsilon,q*skip);
	mpfi_neg(temp,temp);      mpfi_exp(temp,temp);
	mpfi_ui_sub(temp,1,temp); mpfi_inv(Fhatfac,temp);

	mpfi_init(error_term);
}

/* compute truncation error from sample j from mth term on:
error terms are computed by their logarithms, in double precision,
with a small amount added to overcompensate for round-off error */
void truncation_error(int j,double m) {
	double e,X;

	e = (log(m/sqrt((double)N))+2*j*epsd)*2/r;

	/* limit the size of X to prevent overflow */
	if (e > 100.0 || (X=M_PI*r*delta*exp(-delta)*exp(e)) > 1.0e9)
		X = 1.0e9;
	e = logK+(alpha-0.5)*log(m)-X+0.25*r1*r/X+log(1.0+0.5*r/(X-cprime));

	mpfi_set_d(temp,e);
	mpfi_exp(error_term,temp);
	mpfi_neg(temp,error_term);
	mpfi_put(error_term,temp);
}

double labsgamma(double x,double y) /* should be ln(|GAMMA(z)|) */
{
   return(0.0);
};

void Fupper_bound(int i) {
	double A,B,t,logQ,logP,logG,beta;

	B = M_PI/(epsd*skip);
	A = q / B;
	t = i/A;

	/* compute log|Q(1/2+it)| */
	logQ = 0.5*(r0*log(0.25+t*t)+r1*log(2.25+t*t))-r*log(2*M_PI);

	/* compute polar factor */
	logP = 0.0;
	if (mpfi_cmp_ui(residue,0))
		logP = 1.5*log(1.0+2.0/(0.25+t*t));

	/* compute log|gamma(1/2+it)| */
	logG = r0*(labsgamma(0.25,t*0.5)-0.25*log(M_PI))
	     + r1*(labsgamma(0.75,t*0.5)-0.75*log(M_PI));

	mpfi_set_d(temp,r*logZ32+0.5*(logQ+logP)+logG+0.25*M_PI*r*eta/1000*t);
	mpfi_exp(error_term,temp);

	beta = 0.25*M_PI*r-0.5*(r0*atan(0.5/t)+r1*atan(1.5/t))
	     - 4.0/(M_PI*M_PI)*(r0/(t*t-0.25)+r1/(t*t-2.25));
	mpfi_set_d(temp,-B*(beta-t/fabs(t)*0.25*M_PI*r*eta/1000));
	mpfi_set_d(temp,-B*(beta-0.25*M_PI*r*eta/1000));
	mpfi_exp(temp,temp);
	mpfi_ui_sub(temp,1,temp);
	mpfi_div(error_term,error_term,temp);

	mpfi_neg(temp,error_term);
	mpfi_put(error_term,temp);
}

void compute_Lfunction(void) {
	int i,j;
	double m;
	mpfi_t e1,e2,e3,e4;

	init_truncation_error();

	/* set sample points */
	for (i=j=0;j<samples;i++,j+=skip) {
		complex_set(L[i],ftemp[j]);

		/* add truncation error for this sample */
#ifdef error_terms
		m = ceil(range*exp(-2*j*epsd));
		if (m > nmax) m = nmax+1;
		truncation_error(j,m);
		mpfi_add(L[i].re,L[i].re,error_term);
		mpfi_add(L[i].im,L[i].im,error_term);
#endif
	}

	/* pad with zeros */
	for (;i<=q/2;i++)
		complex_set_ui(L[i],0);

	/* add residual terms */
	/* ctemp = residue*sech(Pi*A/2)*exp(I*Pi*r*eta/8) */
	/* e1 = exp(x/2), e2 = 1/e1, e3 = exp((x-Pi*A)/2), e4 = 1/e3 */
	mpfi_const_pi(temp);
	mpfi_mul_ui(temp,temp,r*eta);
	mpfi_div_ui(temp,temp,8*1000);
	mpfi_sin(ctemp.im,temp);
	mpfi_cos(ctemp.re,temp);
	mpfi_init(e1); mpfi_init(e2);
	mpfi_mul_ui(temp,epsilon,skip);
	mpfi_exp(e1,temp); mpfi_inv(e2,e1);
	mpfi_init(e3); mpfi_init(e4);
	mpfi_mul_si(temp,epsilon,-q*skip/2);
	mpfi_exp(e3,temp); mpfi_inv(e4,e3);
	mpfi_add(temp,e3,e4);
	mpfi_div(temp,residue,temp);
	mpfi_mul(ctemp.re,ctemp.re,temp);
	mpfi_mul(ctemp.im,ctemp.im,temp);
	for (i=0;i<=q/2;i++) {
		mpfi_add(temp,e3,e4);
		mpfi_mul(temp,temp,ctemp.re);
		mpfi_sub(L[i].re,L[i].re,temp);
		mpfi_sub(temp,e3,e4);
		mpfi_mul(temp,temp,ctemp.im);
		mpfi_sub(L[i].im,L[i].im,temp);
		mpfi_mul(e3,e3,e1); mpfi_mul(e4,e4,e2);
	}

	/* add pre-FFT discretization error */
#ifdef error_terms
	for (i=j=0;i<=q/2;i++,j+=skip) {
		truncation_error((j < samples) ? j+q*skip : j,1.0);
		mpfi_set(ctemp.re,error_term);
		truncation_error(q*skip-j,1.0);
		mpfi_add(error_term,error_term,ctemp.re);
		mpfi_mul(error_term,error_term,Fhatfac);
		mpfi_add(L[i].re,L[i].re,error_term);
		mpfi_add(L[i].im,L[i].im,error_term);
	}
#endif

	/* Fhat(-x) = conjugate(Fhat(x)) (functional equation) */
	for (;i<q;i++)
		complex_conj(L[i],L[q-i]);

	fft(L,q);

	/* multiply by 2*Pi/B */
	mpfi_mul_ui(temp,epsilon,2*skip);
	for (i=0;i<q;i++) {
		mpfi_mul(L[i].re,L[i].re,temp);
		mpfi_mul(L[i].im,L[i].im,temp);
	}

	/* add post-FFT discretization error */
#ifdef error_terms
	for (i=0;i<q/2;i++) {
		Fupper_bound(i+q);
		mpfi_set(ctemp.re,error_term);
		Fupper_bound(i-q);
		mpfi_add(error_term,error_term,ctemp.re);
		mpfi_add(L[i].re,L[i].re,error_term);
		mpfi_add(L[i].im,L[i].im,error_term);
	}
#endif
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

int main(int argc,char *argv[]) {
	int i,k,n;
	FILE *ffile,*gfile,*fp;
	double Ainv,t;
	char buf[128];

	mpfr_set_default_prec(prec);

	/* read L-function description file */
	if (argc != 4 || !(fp = fopen(argv[1],"r")) ) {
		printf("usage: %s description_file ffile gfile\n",argv[0]);
		return 0;
	}
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%d",&N);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%d",&r0);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%d",&r1);
	r = r0+r1;
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%d",&eta);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%lf",&nmax);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%d",&samples);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%d",&ntaylor);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%d",&q);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%d",&skip);
	ftemp = malloc(2*samples*sizeof(ftemp[0]));
	gtemp = malloc(2*samples*sizeof(gtemp[0]));
	n = 2*samples*ntaylor; if (n < q) n = q;
	L = malloc(n*sizeof(L[0]));
	for (i=0;i<n;i++) complex_init(L[i]);
	fgets(buf,sizeof(buf),fp); mpfi_init(residue);
	mpfr_set_str(&residue->left,buf,10,GMP_RNDD);
	mpfr_set_str(&residue->right,buf,10,GMP_RNDU);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%lf",&C);
	fgets(buf,sizeof(buf),fp); sscanf(buf,"%lf",&alpha);
	fclose(fp);

	/* epsilon = log(range)/(2*samples) */
	complex_init(ctemp); mpfi_init(temp);
	mpfi_init(epsilon);  mpfi_set_d(temp,range);
	mpfi_log(temp,temp); mpfi_div_ui(epsilon,temp,2*samples);
	epsd = mpfi_get_d(epsilon);

fprintf(stderr,"N=%d, r0=%d, r1=%d, eta=%d, nmax=%f\n",N,r0,r1,eta,nmax);
fprintf(stderr,"residue=%f, C=%f, alpha=%f\n",mpfi_get_d(residue),C,alpha);
fprintf(stderr,"ntaylor=%d, samples=%d, q=%d, skip=%d\n",
ntaylor,samples,q,skip);

	ffile = fopen(argv[2],"r");
	gfile = fopen(argv[3],"r");

	for (i=0;i<samples;i++)
		for (k=0;k<ntaylor;k++) {
			/* use same memory for L samples and f,g */
			mpfi_read(L[k*samples+i].re,ffile);
			mpfi_set_ui(L[k*samples+i].im,0);
			complex_read(&L[(k+ntaylor)*samples+i],gfile);
		}
	for (i=2*samples*ntaylor;i<q;i++) complex_init(L[i]);
	mpfi_init(ferror_term);
	mpfi_read(ferror_term,ffile);
	mpfi_init(gerror_term);
	mpfi_read(gerror_term,gfile);
	fclose(ffile); fclose(gfile);

	initfft();
	for (i=0;i<2*samples;i++) complex_init(ftemp[i]);
	for (i=0;i<2*samples;i++) complex_init(gtemp[i]);
fprintf(stderr,"computing samples\n");
	compute_samples();
fprintf(stderr,"computing L-function\n");
	compute_Lfunction();
	Ainv = M_PI/(epsd*q*skip);
	for (i=0;i<q/2;i++) {
		t = i*Ainv;
		//printf("%.14f ",t);
		t = exp(-r0*(labsgamma(0.25,t*0.5)-0.25*log(M_PI))
			-r1*(labsgamma(0.75,t*0.5)-0.75*log(M_PI))
			-0.25*M_PI*r*eta/1000*t);
		mpfi_mul_d(temp,L[i].re,t);
		//printf("%.14g\n",mpfi_get_d(temp));
		mwrite(temp,stdout);
	}

	return 0;
}
