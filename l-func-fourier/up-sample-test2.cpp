/*

File: up-sample-test2.cpp

Created: 3 June 2009

Version: <v> = 1.0

Last Modified: 3 June 2009

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

// file contains Z values for t=0..119+123/128 (3072 values, each two 8 byte doubles)
#define SKIP_IN_POINTS (0)
#define SPACING (16)

#define IN_POINTS (128)

#define T0 ((IN_POINTS-1)/2.0)

// up sample by a factor of 8
#define UP_SAMPLE (32)

#define OUT_POINTS (UP_SAMPLE*IN_POINTS)

#define LAMBDA (0.002)

dcomplex c_zero=dcomplex(0.0,0.0);

double gaussian(unsigned int t)
{
	return(exp(-LAMBDA*(t-T0)*(t-T0)));
}

int main()
{
	unsigned int i,j,k;

	fftw_plan p1,p2;

	dcomplex *signal,*result;

	signal = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * IN_POINTS,16);

	result = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * OUT_POINTS,16);

	// plan p1 is the forward fft on the original data.
	p1 = fftw_plan_dft_1d(IN_POINTS, reinterpret_cast<fftw_complex*>(signal),
		reinterpret_cast<fftw_complex*>(signal), FFTW_FORWARD, FFTW_ESTIMATE);

	// plan p2 is the inverse fft on the padded frequency data
	p2 = fftw_plan_dft_1d(OUT_POINTS, reinterpret_cast<fftw_complex*>(result), 
		reinterpret_cast<fftw_complex*>(result), FFTW_BACKWARD, FFTW_ESTIMATE);

	FILE *infile;
	double left,right;
	infile=fopen("z_821.dat","rb");

	// ignore the first batch of data 
	for(j=0;j<SKIP_IN_POINTS*2;j++)
		fread(&left,sizeof(double),1,infile);

	for(i=0;i<IN_POINTS;i++)
	{
		fread(&left,sizeof(double),1,infile);
		fread(&right,sizeof(double),1,infile);
		signal[i]=dcomplex((left-right)/2.0*gaussian(i),0.0);
		for(j=1;j<SPACING;j++)
		{
			fread(&left,sizeof(double),1,infile);
			fread(&left,sizeof(double),1,infile);
		}
	}

	fclose(infile);

/*
	for(j=0;j<IN_POINTS;j++)
		cout << signal[j].real() << endl;
*/
	fftw_execute(p1); // fft the signal
	
	
	for(j=0;j<IN_POINTS;j++)
		cout << abs(signal[j]) << endl;
	


	i=0;
// put in the +ve frequencies
	for(j=0;j<IN_POINTS/2;j++)
		result[i++]=signal[j];
// zero pad
	for(k=0;k<OUT_POINTS-IN_POINTS;k++)
		result[i++]=c_zero;
// put in the -ve frequencies
	for(;j<IN_POINTS;j++)
		result[i++]=signal[j];

	fftw_execute(p2); // ifft it

/*
	for(i=0;i<OUT_POINTS;i++)
		cout << result[i].real() << endl;
*/		

	free(result);free(signal);

	return(0);
} /* main */
