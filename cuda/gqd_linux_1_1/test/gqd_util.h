
#ifndef _GQD_UTIL_H_
#define _GQD_UTIL_H_

#include "gqd_type.h"
#include <qd/qd_real.h>

void rand(double* val);
void rand(dd_real* val);
void rand(qd_real* val);


void dd_to_gdd( const dd_real* h_in, GPU_dd* d_out, const int numElement );
void qd_to_gqd( const qd_real* h_in, GPU_qd* d_out, const int numElement ); 
void gdd_to_dd( const GPU_dd* h_in, dd_real* h_out, const int numElement );
void gqd_to_qd( const GPU_qd* h_in, qd_real* h_out, const int numElement ); 


int checkTwoArray( const dd_real* gold, const dd_real* ref, const int numElement );
int checkTwoArray( const qd_real* gold, const qd_real* ref, const int numElement );

#endif //_GQD_UTIL_H_
