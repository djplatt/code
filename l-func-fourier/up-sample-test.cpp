/*

File: up-sample.cpp

Created: 13 May 2009

Version: <v> = 1.0

Last Modified: 13 May 2009

Dialect: C++

Requires: fttw3

Implementation notes: 

Build instructions: g++ -lfftw3 -O3

Test of Upsampling principle

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define _CRT_SECURE_NO_WARNINGS  // shut Microsoft up
#define VERSION "1.0"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <fftw3.h>
#include <time.h>
#include <malloc.h>

#define _aligned_malloc(a,b) memalign(b,a)

#define debug printf("got to line %d\n",__LINE__);

using namespace std;

typedef complex<double> dcomplex;

#define N_POINTS (16)

#define LAMBDA1 (-0.0)//(-0.05) // attenuation for signal = exp(lambda*x^2)

#define LAMBDA2 (-0.00001) // attenuation for filter

#define FFT_LEN (256) // size of upsample

#define LEFT_LEAD (1) // this many zeros left of 1st sample

#define GAP (16) // this many zeros between samples


double my_sample[N_POINTS]={1.364856642,1.085844808,0.82122139,0.580146724,
0.371289402,0.202618585,0.081209293,0.013064869,
0.0029605336,0.054311533,0.169068918,0.347645448,
0.588873561,0.889996675,1.246694464,1.653142052};

void convolve(dcomplex *c, dcomplex *a, dcomplex *b, unsigned int len)
{
	unsigned int i;
	fftw_plan p;

	p = fftw_plan_dft_1d(len, reinterpret_cast<fftw_complex*>(a), 
		reinterpret_cast<fftw_complex*>(a), FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	p = fftw_plan_dft_1d(len, reinterpret_cast<fftw_complex*>(b), 
		reinterpret_cast<fftw_complex*>(b), FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	for(i=0;i<len;i++)
		c[i]=a[i]*b[i]/dcomplex(len,0.0);
	p = fftw_plan_dft_1d(len, reinterpret_cast<fftw_complex*>(c), 
		reinterpret_cast<fftw_complex*>(c), FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
}
	

int main()
{
	unsigned int i,fft_ptr;
	double x;

	fftw_plan p1,p2;

	dcomplex *signal,*result,*gauss;

	signal = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * FFT_LEN,16);

	result = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * FFT_LEN,16);

	gauss = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * FFT_LEN,16);

		p1 = fftw_plan_dft_1d(FFT_LEN, reinterpret_cast<fftw_complex*>(signal), 
		reinterpret_cast<fftw_complex*>(signal), FFTW_FORWARD, FFTW_ESTIMATE);

		p2 = fftw_plan_dft_1d(FFT_LEN, reinterpret_cast<fftw_complex*>(signal), 
		reinterpret_cast<fftw_complex*>(result), FFTW_BACKWARD, FFTW_ESTIMATE);


// prepare the Gaussian

	for(i=0;i<N_POINTS/2;i++)
	{
		gauss[i]=dcomplex(exp(LAMBDA2*(N_POINTS-i-i-1)*(N_POINTS-i-i-1)),0.0);
		gauss[N_POINTS-i-1]=gauss[i];
	}

// prepare the signal
	for(i=0;i<FFT_LEN;i++)
		signal[i]=dcomplex(0.0,0.0);

	fft_ptr=LEFT_LEAD-1;
	x=-N_POINTS+1;
	for(i=0;i<N_POINTS;i++)
	{
		signal[fft_ptr]=dcomplex(my_sample[i]*exp(LAMBDA1*x*x),0.0);
		fft_ptr+=GAP;
		x+=2;
	}

	
	for(i=LEFT_LEAD-1;i<FFT_LEN-LEFT_LEAD;i++)
			cout << signal[i].real() << endl;


// convolve them
//	convolve(result,signal,gauss,FFT_LEN);

	fftw_execute(p1); // fft the signal

	for(i=0;i<120;i++)
	{
		signal[i]=dcomplex(0.0,0.0);
		signal[FFT_LEN-i-1]=dcomplex(0.0,0.0);
	}

	fftw_execute(p2); // ifft it


	for(i=0;i<FFT_LEN;i++)
		if(i&1)
			result[i]=-result[i]/dcomplex(16.0,0.0);
		else
			result[i]=result[i]/dcomplex(16.0,0.0);

	for(i=0;i<FFT_LEN;i++)
		cout << result[i].real() << endl;

	for(i=0;i<FFT_LEN;i++)
	if(result[i].real()<0.0)
			printf("i: %ld Re(z): %10.8e\n",i,result[i].real());		
		

	free(result);free(signal);free(gauss);

	return(0);
} /* main */
