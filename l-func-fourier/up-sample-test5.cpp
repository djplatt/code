/*

File: up-sample-test4.cpp

Created: 28 June 2009

Version: <v> = 1.0

Last Modified: 28 June 2009

Dialect: C++

Requires: fttw3

Implementation notes: 

Build instructions: g++ <file> -lfftw3 -O3

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
#define POINTS_IN_FILE (3072)

#define SPACING (2)

#define IN_POINTS (512)

#define SKIP_POINTS (POINTS_IN_FILE-IN_POINTS-(IN_POINTS-1)*SPACING+1)

// up sample by a factor of
#define UP_SAMPLE (32)

#define OUT_POINTS (UP_SAMPLE*IN_POINTS)

#define LAMBDA (0.00002)

dcomplex c_zero=dcomplex(0.0,0.0);

#define T0 ((IN_POINTS-1)/2.0)


double gaussian(unsigned int t)
{
	return(exp(-LAMBDA*(t-T0)*(t-T0)));
}


#define NEG 0
#define POS 1
#define ZER 2

char sign(double x)
{
	if(x<0) return(NEG);
	if(x>0) return(POS);
	printf("Returning Zero from sign.\n");
	return(ZER);
}

#define SACRIFICE (0)

unsigned int count_zeros (dcomplex *vec, unsigned int N)
{
	unsigned int i,count=0;
	char last_sign=sign(vec[(int)(N*SACRIFICE)].real()),this_sign;
	for(i=N*SACRIFICE+1;i<N*(1-SACRIFICE);i++)
	{
		this_sign=sign(vec[i].real());
		if(this_sign==ZER)
			continue;
		if(this_sign!=last_sign)
		{
			count++;
			last_sign=this_sign;
		}
	}
	return(count);
}


int main(int argc, char** argv)
{
	unsigned int i,j,k;
	int choice;
	double x; 
	fftw_plan p1,p2;
	dcomplex *signal,*result;

	if(SKIP_POINTS<0)
	{
		printf("Not enough points in file. SKIP_POINTS=%d. Exiting.\n",SKIP_POINTS);
		exit(0);
	}

	if(argc==2)
		choice=atoi(argv[1]);
	else
	choice=-1;




	signal = (dcomplex*) _aligned_malloc(sizeof(dcomplex) * IN_POINTS *2,16);

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
	for(j=0;j<SKIP_POINTS;j++)
	{
		fread(&left,sizeof(double),1,infile);
		fread(&right,sizeof(double),1,infile);
	}

	for(i=0;i<IN_POINTS-1;i++)
	{
		fread(&left,sizeof(double),1,infile);
		fread(&right,sizeof(double),1,infile);
		signal[i]=dcomplex((left-right)/2.0/**gaussian(i)*/);
		for(k=0;k<SPACING;k++)
		{
			fread(&left,sizeof(double),1,infile);
			fread(&right,sizeof(double),1,infile);
		}
	}
	signal[IN_POINTS-1]=(signal[0]+signal[IN_POINTS-2])/2.0;


	fclose(infile);

	if(choice==3)
	{
		printf("Starting at n=%d.\n",SKIP_POINTS);
		printf("%d zeros in original data\n",count_zeros(signal,IN_POINTS));
	}

	if(choice==0)
	for(j=0;j<IN_POINTS;j++)
		cout << signal[j].real() << endl;

	fftw_execute(p1); // fft the signal

	for(i=0;i<IN_POINTS;i++)
		signal[i]/=(double)IN_POINTS;


	if(choice==1)
		for(j=0;j<IN_POINTS;j++)
		{
			x=abs(signal[j]);
			if(x==0.0)
				cout << 0 << endl;
			else
				cout << log(x) << endl;
		}



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

if(choice==2)
	for(i=0;i<OUT_POINTS;i++)
		cout << result[i].real() << endl;

if(choice==3)
printf("%d zeros found in output.\n",count_zeros(result,OUT_POINTS));


	free(result);free(signal);

	return(0);
} /* main */
