#ifndef _GQD_UTIL_CPP_
#define _GQD_UTIL_CPP_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <memory.h>

#include "gqd_type.h"
#include "gqd_util.h"

using namespace std;


void rand(double* val) {
	(*val) = rand()/(double)RAND_MAX;
}

void rand(dd_real* val) {
	(*val) = ddrand();
}

void rand(qd_real* val) {
	(*val) = qdrand();
}

//copy qd_real type to the GPU_qd type on the GPU
void qd_to_gqd( const qd_real* h_in, GPU_qd* h_out, const int numElement ) {
	for( int i = 0; i < numElement; i++ ) {
		h_out[i].d1.x = h_in[i].x[0];
		h_out[i].d1.y = h_in[i].x[1];
		h_out[i].d2.x = h_in[i].x[2];
		h_out[i].d2.y = h_in[i].x[3];
	}
}

void dd_to_gdd( const dd_real* h_in, GPU_dd* h_out, const int numElement ) {
	for( int i = 0; i < numElement; i++ ) {
		h_out[i].x = h_in[i].x[0];
		h_out[i].y = h_in[i].x[1];
	}
}

//copy GPU_qd type from the GPU to the qd_real type on the CPU

void gqd_to_qd( const GPU_qd* h_in, qd_real* h_out, const int numElement ) {

	for( int i = 0; i < numElement; i++ ) {
		h_out[i].x[0] = h_in[i].d1.x;
		h_out[i].x[1] = h_in[i].d1.y;
		h_out[i].x[2] = h_in[i].d2.x;
		h_out[i].x[3] = h_in[i].d2.y;
	}
}


void gdd_to_dd( const GPU_dd* h_in, dd_real* h_out, const int numElement ) {

	for( int i = 0; i < numElement; i++ ) {
		h_out[i].x[0] = h_in[i].x;
		h_out[i].x[1] = h_in[i].y;
	}
}


//to_string sometimes has problems with nvcc compiler
char* to_char(qd_real a) {
	char* str = new char[1024];
	a.write( str, 1000 );
	return str;
}

char* to_char(dd_real a) {
	char* str = new char[1024];
	a.write( str, 1000 );
	return str;
}

int checkTwoArray( const qd_real* gold, const qd_real* ref, const int numElement ) {
	qd_real maxRelError = 0.0;
	qd_real avgRelError = 0.0;
	int maxId = 0;
	
	for( int i = 0; i < numElement; i++ ) {
		qd_real relError = abs((gold[i] - ref[i])/gold[i]);
		avgRelError += (relError/numElement);
		if( relError > maxRelError ) {
			maxRelError = relError;
			maxId = i;
		}
	}	
	
	cout << "abs. of max. relative error: " << maxRelError << endl;
	cout << "abs. of avg. relative error: " << avgRelError << endl;
	if( maxRelError > 0.0 ) {
		cout << "max. relative error elements" << endl;
		cout << "i = " << maxId << endl;
		cout << "gold = " << (gold[maxId]).to_string() << endl;
		//cout << "ref  = " << to_char(ref[maxId]) << endl;
		cout << "ref  = " << (ref[maxId]).to_string() << endl;
		cout << "rel. error = " << to_char(abs((gold[maxId] - ref[maxId])/gold[maxId])) << endl;
	} else {
		cout << "a sample:" << endl;
		const int i = rand()%numElement;
		cout << "i = " << i << endl;
		cout << "gold = " << to_char(gold[i]) << endl;
		cout << "ref  = " << to_char(ref[i]) << endl;
		maxId = i;
	}

	return maxId;
}


int checkTwoArray( const dd_real* gold, const dd_real* ref, const int numElement ) {
	dd_real maxRelError = 0.0;
	dd_real avgRelError = 0.0;
	int maxId = 0;
	
	for( int i = 0; i < numElement; i++ ) {
		dd_real relError = abs((gold[i] - ref[i])/gold[i]);
		avgRelError += (relError/numElement);
		if( relError > maxRelError ) {
			maxRelError = relError;
			maxId = i;
		}
	}	
	
	cout << "abs. of max. relative error: " << maxRelError << endl;
	cout << "abs. of avg. relative error: " << avgRelError << endl;
	if( maxRelError > 0.0 ) {
		cout << "max. relative error elements" << endl;
		cout << "i = " << maxId << endl;
		cout << "gold = " << to_char(gold[maxId]) << endl;
		printf("Components(%.16e, %.16e)\n", gold[maxId].x[0], gold[maxId].x[1]);
		cout << "ref  = " << to_char(ref[maxId]) << endl;
		printf("Components(%.16e, %.16e)\n", ref[maxId].x[0], ref[maxId].x[1]);

	} else {
		cout << "a sample:" << endl;
		const int i = rand()%numElement;
		cout << "i = " << i << endl;
		cout << "gold = " << to_char(gold[i]) << endl;
		printf("Components(%.16e, %.16e)\n", gold[maxId].x[0], gold[maxId].x[1]);
		cout << "ref  = " << to_char(ref[i]) << endl;
		printf("Components(%.16e, %.16e)\n", ref[maxId].x[0], ref[maxId].x[1]);
		maxId = i;
	}

	return maxId;
}

int checkTwoArray( const double* gold, const double* ref, const int numElement ) {
	double maxRelError = 0.0;
	double avgRelError = 0.0;
	int maxId = 0;
	
	for( int i = 0; i < numElement; i++ ) {
		double relError = abs((gold[i] - ref[i])/gold[i]);
		avgRelError += (relError/numElement);
		if( relError > maxRelError ) {
			maxRelError = relError;
			maxId = i;
		}
	}	
	
	cout << "abs. of max. relative error: " << maxRelError << endl;
	cout << "abs. of avg. relative error: " << avgRelError << endl;
	if( maxRelError > 0.0 ) {
		cout << "max. relative error elements" << endl;
		cout << "i = " << maxId << endl;
		//cout << "gold = " << (gold[maxId]) << endl;
		//cout << "ref  = " << (ref[maxId]) << endl;
		printf( "gold = %.15f\n", gold[maxId] );
		printf( "ref  = %.15f\n", ref[maxId] );
	} else {
		cout << "a sample:" << endl;
		const int i = rand()%numElement;
		cout << "i = " << i << endl;
		//cout << "gold = " << (gold[i]) << endl;
		//cout << "ref  = " << (ref[i]) << endl;
		printf( "gold = %.15f\n", gold[i] );
		printf( "ref  = %.15lf\n", ref[i] );
		maxId = i;
	}

	return maxId;
}

#endif
