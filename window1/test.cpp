/*
  File: test.cpp


  By: DJ Platt
  Bristol University

  Copyright 2010.

  This work is funded by the UK ESPRC. */

#define VERSION "1.0"
//#define PRINT
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
using namespace std;



#include "../includes/int_double11.0.h"

#define one_over_A ((double) 5.0/64.0)
#define MAX_SIGMA (50)
#define MIN_SIGMA (10)

void print_usage()
/* called when wrong arguments passed via command line */
{
  printf("usage:- test <q> <N> <t0> <M>\n");
  exit(1);
}


void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};

int_complex *init_ws(int N)
{
  int_complex *ws;
  if(!(ws=(int_complex *) _aligned_malloc(sizeof(int_complex)*N/2,16)))
    fatal_error("Failed to allocate memory for ws. Exiting.\n");
  ws[0]=c_one;
  int_double two_pi_n=-d_two_pi/N;
  for(int i=1;i<N/2;i++)
      sin_cos(two_pi_n*i,&ws[i].imag,&ws[i].real);
  return(ws);
}

// perform an in place FFT ,n a power of 2
// F[k]=sum_{j=0..N-1} F[j]e(-kj/N)
void fft(int_complex *x,unsigned int n,int_complex *w)
{
	unsigned int i,j,k,l;
	int_complex *p,*xend=x+n,ctemp;
	//printf("In FFT with n=%d\n");
	//for(i=0;i<n/2;i++)
	//print_int_complex_str("w=",w[i]);
	//for(i=0;i<n;i++)
	//print_int_complex_str("fft=",x[i]);
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

/* perform an in place inverse FFT */
void ifft(int_complex *x,unsigned int n,int_complex *w) {
	unsigned int i,l=n>>1;
	int_complex ctemp;

	fft(x,n,w);
	//x[0]=x[0];
	//x[l]=x[l];

	for (i=1;i<l;i++) {
	  //x[i]=x[i]/n;
	  //x[n-i]=x[n-i]/n;
		ctemp=x[i];
		x[i]=x[n-i];
		x[n-i]=ctemp;
	}
}


/* circular convolution; f and g must be distinct */
void convolve(int_complex *result,int_complex *f,int_complex *g,unsigned int n, int_complex *ws)
{
	unsigned int i;
	fft(f,n,ws);
	fft(g,n,ws);
	for (i=0;i<n;i++)
		result[i]=f[i]*g[i];
	ifft(result,n,ws);
}

double calc_eta(double t0, int sigma)
{
  int_complex lng=lngamma(int_complex(int_double(sigma),int_double(t0))/2.0);
  return(-lng.real.left*4.0/M_PI/t0);
}

int_complex do_sum(double t, int sigma, int q, const int_double &log_q_by_pi, int_complex *chis, bool *co_primes, int M, int_double *logs)
{
  int_complex res=c_zero,term;
  for(int n=1;n<=M;n++)
    {
      if(!co_primes[n%q])
	continue;
      sin_cos((log_q_by_pi-logs[n]*2)*t/2.0,&term.imag,&term.real);
      res+=term*exp(logs[n]*(-sigma))*chis[n%q];
    }
  return(res);
}

int_complex *_ws,*fft_vec;

int main(int argc, char **argv)
{

  _fpu_rndd();

  int q,N,i,j;
  clock_t no_clicks;
  double B;

  no_clicks=clock(); // start timing
  if(argc!=5)
    print_usage();
  
  q=atoi(argv[1]);
  if((q<3)||((q&3)==2))
    print_usage();
  q=5;
  printf("q=%d\n",q);

  i=atoi(argv[2]);
  if((i<1)||(i>30))
    print_usage();

  for(N=1,j=0;j<i;j++,N<<=1);
  printf("N=%d\n",N);

  B=N*one_over_A;
  printf("B=%f\n",B);

  double t0=atof(argv[3]);
  if(t0<0.0)
    print_usage();
  printf("t0=%f\n",t0);

  /*
  int sigma=atoi(argv[4]);
  if(sigma<=1)
    print_usage();
  printf("sigma set to %d\n",sigma);
  */
  int M=atoi(argv[4]);
  if(M<2)
    print_usage();
  printf("M=%d\n",M);

  double h=B/5.0;
  printf("h=%f\n",h);

  int sigma=2*h;
  if(sigma<MIN_SIGMA)
    {
      printf("B two narrow to allow sensible sigma (set to %d). Try increasing N. Exiting.\n",sigma);
      exit(1);
    }
  else
    {
      if(sigma>MAX_SIGMA)
	sigma=MAX_SIGMA;
      printf("sigma=%d\n",sigma);
    }

  int_double sum_error=exp(log(int_double(M))*(1-sigma))/sigma;
  sum_error.left=sum_error.right;
  print_int_double_str("sum_error=",sum_error);

  //N=64;
  printf("Initialising ws...\n");
  _ws=init_ws(N);

  if(!(fft_vec=(int_complex *)_aligned_malloc(sizeof(int_complex)*N,16)))
    fatal_error("Failed to allocate memory for fft_vec. Exiting.\n");
  /*
  for(i=0;i<N;i++)
    fft_vec[i]=c_zero+((i*i)%29);
  for(i=0;i<N;i++)
    print_int_complex_str("",fft_vec[i]);
  fft(fft_vec,N,_ws);
  for(i=0;i<N;i++)
    print_int_complex_str("",fft_vec[i]);
  ifft(fft_vec,N,_ws);
  for(i=0;i<N;i++)
    print_int_complex_str("",fft_vec[i]);
  exit(0);
  */
  int_complex *chis;
  if(!(chis=(int_complex *) _aligned_malloc(sizeof(int_complex)*q,16)))
    fatal_error("Failed to allocate memory for chis. Exiting.\n");

  chis[0]=c_zero;
  chis[1]=c_one;
  chis[2]=-c_one;
  chis[3]=-c_one;
  chis[4]=c_one;

  bool *co_primes;

  if(!(co_primes=(bool *) malloc(sizeof(bool)*q)))
    fatal_error("Failed to allocate memory for co_primes. Exiting.\n");

  co_primes[0]=false;
  for(i=1;i<q;i++)
    co_primes[i]=true;

  int_complex epsilon=c_one;

  double eta=1.0;//calc_eta(t0,sigma);
  printf("eta=%f\n",eta);

  int_double *logs;
  if(!(logs=(int_double *) _aligned_malloc(sizeof(int_double)*M,16)))
    fatal_error("Failed to allocate memory for logs. Exiting.\n");

  for(i=1;i<q;i++)
    if(co_primes[i])
      logs[i]=log(int_double(i));

  double t,t_t0;
  int_complex sum;
  int_complex lng;
  int_double log_q_by_pi=log(int_double(q)/d_pi);
  int_double two_h_sqr=int_double(h)*h*2.0;


  for(t=0.0,i=0,t_t0=t0;i<N;i++,t+=one_over_A,t_t0+=one_over_A)
    {
      //t=0.0;t_t0=t0;
      //printf("t=%f t+t0=%f\n",t,t_t0);
      sum=do_sum(t_t0,sigma,q,log_q_by_pi,chis,co_primes,M,logs)+int_complex(sum_error,sum_error);
      sum*=epsilon;
      lng=lngamma(int_complex(int_double(sigma),int_double(t_t0))/2.0);
      //print_int_complex_str("lng=",lng);
      lng.real+=d_pi*eta*(t_t0)/4.0;
      //lng.real+=log(int_double(t)+t0)*(sigma-1.0)/2.0;
      //print_int_double_str("Re(lng)+pi*eta*(t+t0)/4)=",lng.real);
      lng.real+=log_q_by_pi*(sigma-0.5)/2.0;
      lng.real+=int_double(sigma-0.5)*(sigma-0.5)/two_h_sqr;
      lng.real-=sqr(int_double(t)-B/2.0)/two_h_sqr;
      //print_int_double_str("Re(log(res))=",lng.real);
      lng.imag+=d_pi*eta*(0.5-sigma)/4.0;
      lng.imag+=int_double(t)*2.0*(sigma-0.5)/two_h_sqr;
      lng.imag+=int_double(0.5-sigma)*B/two_h_sqr;
      lng=exp(lng)*sum;
      //print_int_complex_str("res=",lng);
      fft_vec[i]=lng;
      //exit(0);
    }
  // pretend the above is 100% accurate
  /*
  for(i=0;i<N;i++)
    {
      fft_vec[i].real.right=-fft_vec[i].real.left;
      fft_vec[i].imag.right=-fft_vec[i].imag.left;
    }
  */
  printf("Starting FFT...\n");
  fft(fft_vec,N,_ws);
  printf("FFT completed.\n");
  for(i=1;i<N;i<<=1)
    {
      print_int_complex_str("fft=",fft_vec[i]);
      print_int_complex_str("fft=",fft_vec[N-i]);
      printf("\n");
    }

  for(i=1;i<N/2;i++)
    {
      fft_vec[i]*=exp(d_pi*i/B*(1-2*sigma));
      fft_vec[N-i]=conj(fft_vec[i]);
    }
  fft_vec[N/2]*=exp(d_pi*N/2/B*(1-2*sigma));
  printf("After * exp\n");
  print_int_complex_str("fft[0]=",fft_vec[0]);
  print_int_complex_str("fft[N/2]=",fft_vec[N/2]);
  for(i=1;i<N;i<<=1)
    {
      print_int_complex_str("fft=",fft_vec[i]);
      print_int_complex_str("fft=",fft_vec[N-i]);
      printf("\n");
    }
  ifft(fft_vec,N,_ws);
  printf("After iFFT\n");
  for(i=0;i<N;i++)
    {
      printf("%d ",i);print_int_complex_str("fft=",fft_vec[i]);
      //printf("\n");
    }


  return(0);
}
