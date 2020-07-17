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

// this many leading and trailing zeros
#define FRONT_ZEROS (79)

// this many zeros between points
#define GAP_ZEROS (30)

// number of points to process
#define IN_POINTS (128)

// position of "central point"
#define T0 ((IN_POINTS-1)/2.0)

// number of points to output
#define OUT_POINTS (4096)

// lambda value for Gaussian
// larger makes for more damping.
#define LAMBDA (0.002)

// c_zero <-- 0.0+0.0i
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

	if(!(signal = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * OUT_POINTS,16)))
		exit(0);
	if(!(result = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * OUT_POINTS,16)))
		exit(0);

	// plan p1 is the forward fft on the padded original data.
	if(!(p1 = fftw_plan_dft_1d(OUT_POINTS, reinterpret_cast<fftw_complex*>(signal),
			reinterpret_cast<fftw_complex*>(result), FFTW_FORWARD, FFTW_ESTIMATE)))
			exit(0);

	// plan p2 is the inverse fft on the low passed frequency data
	if(!(p2 = fftw_plan_dft_1d(OUT_POINTS, reinterpret_cast<fftw_complex*>(result), 
		reinterpret_cast<fftw_complex*>(signal), FFTW_BACKWARD, FFTW_ESTIMATE)))
		exit(0);

	FILE *infile;
	double left,right;
	infile=fopen("z_821.dat","rb");

	// ignore the first batch of data 
	for(j=0;j<SKIP_IN_POINTS*2;j++)
		fread(&left,sizeof(double),1,infile);

	for(i=0;i<FRONT_ZEROS;i++)
		signal[i]=c_zero;
	for(j=0;j<IN_POINTS-1;j++)
	{
		fread(&left,sizeof(double),1,infile);
		fread(&right,sizeof(double),1,infile);
		signal[i++]=dcomplex((left-right)/2.0*gaussian(j),0.0);
		for(k=0;k<GAP_ZEROS;k++)
			signal[i++]=c_zero;
	}
		fread(&left,sizeof(double),1,infile);
		fread(&right,sizeof(double),1,infile);
		signal[i++]=dcomplex((left-right)/2.0*gaussian(j),0.0);

	for(j=0;j<FRONT_ZEROS;j++)
		signal[i++]=c_zero;

	fclose(infile);

/*
	for(j=0;j<OUT_POINTS;j++)
		cout << signal[j].real() << endl;
*/
	fftw_execute(p1); // fft the signal
/*	
	for(j=0;j<OUT_POINTS;j++)
		cout << abs(result[j]) << endl;
*/
// get rid of the high frequencies
// keep only IN_POINTS of the lowest
	for(j=IN_POINTS/2;j<OUT_POINTS/2;j++)
	{
		result[j]=c_zero;
		result[OUT_POINTS-j]=c_zero;
	}
/*	
	for(j=0;j<OUT_POINTS;j++)
		cout << abs(result[j]) << endl;
*/

	fftw_execute(p2); // ifft it


	for(i=0;i<OUT_POINTS;i++)
		cout << signal[i].real() << endl;


	free(signal);free(result);

	return(0);
} /* main */
