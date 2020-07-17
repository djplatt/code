#ifndef _GPU_QD_TYPE_CU_
#define _GPU_QD_TYPE_CU_

#include "gqd_type.h"


/** make functions */
__device__ __host__
GPU_dd make_dd( const double x, const double y ) {
	return make_double2( x, y );
}

__device__ __host__
GPU_dd make_dd( const double x ) {
	return make_double2( x, 0.0 );
}

__device__ __host__
GPU_qd make_qd( const double x, 
		const double y, 
		const double z, 
		const double w ) {
	GPU_qd a;
	a.d1.x = x;
	a.d1.y = y;
	a.d2.x = z;
	a.d2.y = w;
	return a;
}

__device__ __host__
GPU_qd make_qd( const double x ) {
	return make_qd( x, 0.0, 0.0, 0.0 );
}


/*
__inline__ __host__ __device__
void make_real( double x, GPU_dd *dd_x )
{
	dd_x[0] = make_dd( x );
}


__inline__ __host__ __device__
void make_real( double x, GPU_qd* qd_x )
{
	qd_x[0] = make_qd( x );
}

*/


//double-double
#define GPU_D_EPS (1.11e-16)


#define TEST_MAX (1<<30)
//quad-double
#define _qd_e make_qd(2.718281828459045091e+00, 1.445646891729250158e-16,  -2.127717108038176765e-33, 1.515630159841218954e-49)
#define _qd_log2 make_qd(6.931471805599452862e-01, 2.319046813846299558e-17,5.707708438416212066e-34,-3.582432210601811423e-50)
#define _qd_eps 1.21543267145725e-63 // = 2^-209
#define _qd_2pi make_qd(6.283185307179586232e+00, 2.449293598294706414e-16, -5.989539619436679332e-33, 2.224908441726730563e-49)
#define _qd_pi make_qd(3.141592653589793116e+00, 1.224646799147353207e-16, -2.994769809718339666e-33, 1.112454220863365282e-49)
#define _qd_pi2 make_qd(1.570796326794896558e+00, 6.123233995736766036e-17, -1.497384904859169833e-33, 5.562271104316826408e-50)
#define _qd_pi1024 make_qd( 3.067961575771282340e-03, 1.195944139792337116e-19,  -2.924579892303066080e-36, 1.086381075061880158e-52)

//double-double
#define _dd_eps 4.93038065763132e-32  // 2^-104
#define _dd_e make_dd(2.718281828459045091e+00, 1.445646891729250158e-16)
#define _dd_log2 make_dd(6.931471805599452862e-01, 2.319046813846299558e-17)
#define _dd_2pi make_dd(6.283185307179586232e+00, 2.449293598294706414e-16)
#define _dd_pi make_dd(3.141592653589793116e+00, 1.224646799147353207e-16)
#define _dd_pi2 make_dd(1.570796326794896558e+00, 6.123233995736766036e-17)
#define _dd_pi16 make_dd(1.963495408493620697e-01, 7.654042494670957545e-18)



#define _QD_SPLITTER 134217729.0               // = 2^27 + 1
#define _QD_SPLIT_THRESH 6.69692879491417e+299 // = 2^996

//the precision
#define _dd_ndigits 31
#define _qd_ndigits 62

#endif

