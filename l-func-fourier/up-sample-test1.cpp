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

#define N_POINTS (512)

#define FFT_LEN (N_POINTS*N_POINTS) // size of upsample

#define GAP (N_POINTS) // this many zeros between samples


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
	unsigned int i,j,fft_ptr;
	double x;

	fftw_plan p1,p2;

	dcomplex *signal,*result,*gauss;

	signal = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * FFT_LEN,16);

	result = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * FFT_LEN,16);

		p1 = fftw_plan_dft_1d(FFT_LEN, reinterpret_cast<fftw_complex*>(signal), 
		reinterpret_cast<fftw_complex*>(signal), FFTW_FORWARD, FFTW_ESTIMATE);

		p2 = fftw_plan_dft_1d(FFT_LEN, reinterpret_cast<fftw_complex*>(signal), 
		reinterpret_cast<fftw_complex*>(result), FFTW_BACKWARD, FFTW_ESTIMATE);



	FILE *infile;
	double left,right;
	infile=fopen("z_821.dat","rb");

	for(j=0;j<1024*5;j++)
		fread(&left,sizeof(double),1,infile);

	for(j=0;j<1;j++){
	for(i=0;i<FFT_LEN;i++)
		signal[i]=dcomplex(0.0,0.0);

	fft_ptr=0;
	for(i=0;i<N_POINTS;i++)
	{
		fread(&left,sizeof(double),1,infile);
		fread(&right,sizeof(double),1,infile);
		signal[fft_ptr]=dcomplex((left-right)/2.0,0.0);
		fft_ptr+=GAP;
	}
	
//	for(i=0;i<FFT_LEN;i++)
//			cout << signal[i].real() << endl;


// convolve them
//	convolve(result,signal,gauss,FFT_LEN);

	fftw_execute(p1); // fft the signal
	for(i=0;i<(N_POINTS*(N_POINTS-1)/2);i++)
	{
		signal[i]=dcomplex(0.0,0.0);
		signal[FFT_LEN-i-1]=dcomplex(0.0,0.0);
	}
	fftw_execute(p2); // ifft it


	for(i=0;i<FFT_LEN;i++)
		if(i&1)
			result[i]=-result[i]/dcomplex(N_POINTS,0.0);
		else
			result[i]=result[i]/dcomplex(N_POINTS,0.0);

	for(i=0;i<FFT_LEN;i++)
		cout << result[i].real() << endl;
	}
		

	free(result);free(signal);

	return(0);
} /* main */
