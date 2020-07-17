//
// win_zeta_dbl.cpp
//
// Windowed FFT based Zeta-function calculator
// See Booker - Artin, Turing ....
//
// Vesrion 1.0 Initial implementation
//
// Created: 1 November 2010
// Last Modified: 1 November 2010
//
// DJ Platt
// University of Bristol
//

#include "stdlib.h"
#include "stdio.h"
#include "assert.h"
#include "../includes/int_double12.0.h"

#define debug printf("Reached line number %d\n",__LINE__)

#define N ((int) 1<<18)
#define one_over_A ((double) 41.0/1024.0)
#define h ((double) 176.6036422)
#define B ((double) N*one_over_A)
#define M ((int) 90000)
#define log2M ((int) 17)  // 2^log2M>=M

#define K ((int) 26)
#define gtwiderr_d ((double) 8.3e-104)
#define Gtwiderr_d ((double) 1e-104) // actually much smaller
#define fhatsumerr_d ((double) 2.8e-57)
#define tayerr_d ((double) 1.2e-45)
#define fhattwiderr_d ((double) 1e-104) // actually much smaller
#define ftwiderr_d ((double) 1e-104) // ditto

//
// s_vec[n-1] contains log(n sqrt(Pi))/2Pi)-u_m n=1..M
// sk_vec[n-1] contains (log(n sqrt(Pi))/2Pi)-u_m)^k
// buck[n] tells us which sk to put s_vec[n] into
// n_vec[n-1] contains 1/sqrt(n)*(nsqrt(Pi))^(-it0) n=1..M
// g_vec contains g(t)*(2*Pi*I*t)^k/k!
// two_pi_t contains (-2*Pi*t)
// G_vec contains (-2*Pi*i*t)^k*g(t)/k!, then G^(k)(x)/k!, then sum m=1..M G^(k)(x+u_m)S^(k)_m/k!
// f_vec contains f^(x), then f(t)
// sqrts[i] contains (i+1)^-0.5
//

int_double A,two_pi_t[N],s_vec[M],sk_vec[M],exps[N],sqrts[M],logisqrtpi[M],d_sqrt_pi;
int_complex g_vec[N],G_vec[N],f_vec[N],skn_vec[N],ws_r[N/2],ws_f[N/2],n_vec[M];
int_complex gtwiderr,Gtwiderr,c_tmp,fhatsumerr,fhattwiderr,ftwiderr,tayerr;
int buck[M]; // which bucket does log(nsqrt(Pi)/2pi) go in.

void print_usage()
{
  printf("Usage:- win_zeta_dbl t0 t1 step outfile.\n");
  exit(1);
}

void print_vec (char *str, int_complex *vec)
{
  for(int i=0;i<N;i+=N/64)
    {
      printf("%s[%d]=",str,i);
      print_int_complex(vec[i]);
      printf("\n");
    }
}

// perform an in place FFT ,n a power of 2
// F[k]=sum_{j=0..N-1} F[j]e(-kj/N)
void fft(int_complex *x,unsigned int n,int_complex *w)
{
	unsigned int i,j,k,l;
	int_complex *p,*xend=x+n,ctemp;

	/* swap each element with one with bit-reversed index */
	for (i=0,l=n>>1;i<l;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k)
				j |= 1;
		}
		if (i < j)
		{
			ctemp=x[i];
			x[i]=x[j];
			x[j]=ctemp;
		}
		else if (i > j)
		{
			ctemp=x[n-1-i];
			x[n-1-i]=x[n-1-j];
			x[n-1-j]=ctemp;
		}
		++i, j |= l;
		ctemp=x[i];
		x[i]=x[j];
		x[j]=ctemp;
	}

	for (k=1,l=n/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<n/2;j+=l,p++) {
				ctemp=p[k]*w[j];
				p[k]=p[0]-ctemp;
				p[0]+=ctemp;
			}
}


int_complex set_err(double d)
{
  return(int_complex(int_double(-d,d),int_double(-d,d)));
}

// reverse skn_vec
inline int conj_order(int i)
{
  return(N-i);
}

// on entry s_vec[i] contains log((i+1)Pi)/2Pi
// on exit                    "" - u_m
void calc_buck()
{
  int i;
  double x;
  int_double tmp;
  //printf("s_vec[0]=");mpfi_printn(s_vec[0],10);
  for(i=0;i<M;i++)
    {
      buck[i]=s_vec[i].left*B+0.5;
      s_vec[i]-=int_double(buck[i])/B;
    }
  //printf("Before shuffle, n=1 term is in bucket %d\n",buck[0]);
  assert(buck[0]>0); // keep conj_order simple
  assert(buck[M-1]<N/2); // ditto
  for(i=0;i<M;i++)
    buck[i]=conj_order(buck[i]);
  //printf("buck[0]=%d\nbuck[M-1]=%d\n",buck[0],buck[M-1]);
  //printf("n=1 term is in bucket %d\n",buck[0]);
}
      
      
// on exit g[n]=g((n-N/2)/A;0)
void init(double t0)
{
  int i;
  double t;
  int_complex z=int_complex(int_double(0.25),d_zero);
  d_sqrt_pi=sqrt(d_pi);
  for(i=0,t=-N/2*one_over_A;i<N;i++,t+=one_over_A)
    {
      z.imag=int_double((t+t0)/2.0);
      //if ((i&4095)==0)
      //printf("t=%e\n",t);
      two_pi_t[i]=d_two_pi*(-t);
      exps[i]=int_double(t)/h;
      exps[i]=exps[i]*exps[i]/2.0;
      g_vec[i]=lngamma(z)+d_pi*(t+t0)/4-exps[i];
      //print_int_complex_str("before exp=",g_vec[i]);
      g_vec[i]=exp(g_vec[i]);
      //print_int_complex_str("after exp=",g_vec[i]);
      //exit(0);
    }
  //print_vec("g_vec",g_vec);
  //exit(0);
  for(i=0;i<M;i++)
    {
      sqrts[i]=sqrt(int_double(i+1));
      s_vec[i]=log(d_sqrt_pi*(i+1));
      logisqrtpi[i]=s_vec[i];
      sin_cos(s_vec[i]*(-t0),&n_vec[i].imag,&n_vec[i].real);
      n_vec[i]*=sqrts[i];
      n_vec[i]/=d_two_pi;
    }
  
  int_double omega=d_two_pi/N;
  for(i=0;i<N/2;i++)
    {
      sin_cos(omega*i,&ws_r[i].imag,&ws_r[i].real);
      ws_f[i]=conj(ws_r[i]);
    }
  gtwiderr=set_err(gtwiderr_d);
  Gtwiderr=set_err(Gtwiderr_d);
  tayerr=set_err(tayerr_d);
  tayerr*=sqrt(int_double(M))*2.0-1.0;
  fhatsumerr=set_err(fhatsumerr_d);
  fhattwiderr=set_err(fhattwiderr_d);
  ftwiderr=set_err(ftwiderr_d);
  A=d_one/one_over_A;
  calc_buck();

}

// run with a new value of t0
// just need to set g_vec and n_vec up
void re_init(double t0)
{
}

// on entry g_vec[i]=g(t)*(-2*Pi*t*I)^k/k!=g(t;k)/k!
// on exit G_vec[i]=G^(k)(i/B)/k!
void G_k(int k)
{
  int i;
  double t,max_t;

  for(i=0;i<N;i++)
    G_vec[i]=g_vec[i]+gtwiderr;

  fft(G_vec,N,ws_f);

  for(i=0;i<=N/2;i++)
    {
      G_vec[i]=G_vec[i]/A+Gtwiderr;
      if(i&1)
	G_vec[i]=-G_vec[i];
    }
  for(i=N/2+1;i<N;i++)
    G_vec[i]=c_zero;

  if(k<K-1)
    for(i=0;i<N;i++)
      {
	int_double a=two_pi_t[i]/(k+1),b;
	b=g_vec[i].real*a;
	g_vec[i].real=-g_vec[i].imag*a;
	g_vec[i].imag=b;
      }
}


// on entry sk_vec = (log(nsqrt(Pi))/2Pi-u_m)^k (or undefined if k==0)
// n_vec[i] = (i+1)^(-1/2)*(sqrt(Pi)*(i+1))^-it0
// on exit skn_vec = sum nsqrt(pi)^(-ito)/sqrt(n)*(log(nsqrt(Pi))/2Pi-u_m)^k
//         sk_vec = (log..-u_m)^(k+1)
void make_skn(int k)
{
  int i,pos;
  for(i=0;i<N;i++)
    skn_vec[i]=c_zero;
      
  switch(k)
    {
    case 0: // set sk_vec = s_vec, skn_vec=n_vec
      for(i=0;i<M;i++)
	{
	  skn_vec[buck[i]]+=n_vec[i];
	  sk_vec[i]=s_vec[i];
	}
      break;
    case K-1: // no need to increment sk_vec
      for(i=0;i<M;i++)
	{
	  skn_vec[buck[i]]+=n_vec[i]*sk_vec[i];
	}
      break;
    default:
      for(i=0;i<M;i++)
	{
	  skn_vec[buck[i]]+=n_vec[i]*sk_vec[i];
	  sk_vec[i]*=s_vec[i];
	}
    }
}

void my_convolve (int_complex *res, int_complex *v1, int_complex *v2, int n, int_complex *ws_r, int_complex *ws_f)
{
  int i;
  fft(v1,n,ws_r);
  fft(v2,n,ws_r);
  for(i=0;i<n;i++)
    res[i]=v1[i]*v2[i];
  fft(res,n,ws_f);
  for(i=0;i<n;i++)
    res[i]/=N;
}

// f_vec+=G(k)(x+u_m)/k!*S^(k)_m
void do_conv (int k)
{
  int i;
  make_skn(k);
    
  if(k==0)
    my_convolve(f_vec,skn_vec,G_vec,N,ws_r,ws_f);
  else
    {
      my_convolve(G_vec,skn_vec,G_vec,N,ws_r,ws_f);
      for(i=0;i<=N/2;i++)
	f_vec[i]+=G_vec[i];
    }
}

int main(int argc, char **argv)
{
  int prec;
  int i,j,k,int_step;
  double t0,t1,step,t;
  FILE *outfile;
  int_double tmp;
  bool first=true;

  if(argc!=5)
    print_usage();

  _fpu_rndd();  

  outfile=fopen(argv[4],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[5]);
      exit(1);
    }

  t0=atof(argv[1]);
  t1=atof(argv[2]);
  step=atof(argv[3]);
  int_step=step/one_over_A;
  assert((int_step&1)==0); // must be even
  assert(step/one_over_A-(double) int_step ==0.0); // must be an exact number of steps
  printf("Aiming to get %d values per run.\n",int_step);

  //main loop
  while(t0<=t1)
    {
      printf("Running centred at t0=%f.\n",t0);
      if(first)
	{
	  first=false;
	  printf("Calling init.\n");
	  system("date");
	  init(t0);
	  printf("Init finished.\n");
	  system("date");
	}
      else
	{
	  printf("Calling re-init.\n");
	  system("date");
	  re_init(t0);
	  printf("re-init finished.\n");
	  system("date");
	}

      for(k=0;k<K;k++)
	{
	  printf("Processing k=%d\n",k);
	  G_k(k);
	  //print_vec("G_vec",G_vec);
	  do_conv(k);
	  //print_vec("F_vec",f_vec);
	  //exit(0);
	}
      printf("convolutions finished\n");
      system("date");

      for(i=0;i<=N/2;i++)
	f_vec[i]+=fhatsumerr+tayerr+fhattwiderr;

      for(i=N/2+1;i<N;i++)
	f_vec[i]=conj(f_vec[N-i]);

      for(i=1;i<N;i+=2)
	f_vec[i]=-f_vec[i];

      printf("Final iFFT\n");
      fft(f_vec,N,ws_r);
      printf("Final iFFT finished.\n");
      system("date");

      /*
      for(i=N/2-int_step/2,t=-step/2.0;i<=N/2+int_step/2;i++,t+=one_over_A)
	{
	  mpfi_div_d(f_vec[i]->re,f_vec[i]->re,B);
	  mpfi_add(f_vec[i]->re,f_vec[i]->re,ftwiderr->re);
       
	  // now get rid of the Gaussian
	  mpfi_set_d(tmp,t);
	  mpfi_div_d(tmp,tmp,h);
	  mpfi_sqr(tmp,tmp);
	  mpfi_div_ui(tmp,tmp,2);
	  mpfi_exp(tmp,tmp);
	  mpfi_mul(tmp,f_vec[i]->re,tmp);

	  // this is it!
	  printf("f(%10f)=",t0+t);mpfi_printn(tmp,30);
	}
      */

            
      for(i=0,t=-one_over_A*N/2;i<=N/2;i+=N/256,t+=one_over_A*N/256)
	{
	  f_vec[i].real/=B;
	  f_vec[i].real+=ftwiderr.real;
	  f_vec[i]*=exp(exps[i]);
	  printf("%6d: Rel error = %4d Abs error = %4d f(%18.10f)=",
		 N/2-i,rel_error(f_vec[i].real), abs_error(f_vec[i].real),t+t0);
	  print_int_double_str("",f_vec[i].real);
	}
      

      t0+=step;
    }
  fclose(outfile);
}
