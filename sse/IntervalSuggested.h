/*

  IntervalSuggested.h

  Sample implementation of the additional instuctions suggested to
  speed up interval computations with SSE-2 and similar extensions.

  For more information about these functions, refer to the paper
  "Interval arithmetic using SSE-2", downloadable from
  http://www.brics.dk/~barnie/RealLib/

*/

#ifndef FILE_INTERVAL_SUGGESTED_H
#define FILE_INTERVAL_SUGGESTED_H

#ifdef _MSC_VER
#include <emmintrin.h>
#else
#include <xmmintrin.h>
#endif

// implementation of the _mm_shuffle_pd operation
// with a non-const shuffle parameter
static inline
__m128d myshuffle(const __m128d a, const __m128d b, int m)
{
	switch(m) {
		case 0:
			return _mm_shuffle_pd(a, b, 0);
		case 1:
			return _mm_shuffle_pd(a, b, 1);
		case 2:
			return _mm_shuffle_pd(a, b, 2);
		case 3:
		default:
			return _mm_shuffle_pd(a, b, 3);
	}
}

// ivchoice: choice instruction for the preparation
// of parameters for interval multiplication
// (used in IntervalMul and IntervalDiv)
static inline 
__m128d ivchoice(__m128d a, const __m128d b)
{
    a = _mm_xor_pd(a, _mm_set_pd(-0.0, 0.0));
    a = myshuffle(a, a, _mm_movemask_pd(b));
    return a;
}

// ivmul: fusion of ivchoice and multiplication
// (used in IntervalMul2)
static inline
__m128d ivmul(__m128d a, const __m128d b)
{
    a = _mm_xor_pd(a, _mm_set_pd(-0.0, 0.0));
    a = myshuffle(a, a, _mm_movemask_pd(b));
	a = _mm_mul_pd(a, b);
    return a;
}

// ivdiv: fusion of ivchoice and division
// (used in IntervalDiv2)
static inline
__m128d ivdiv(__m128d a, const __m128d b)
{
    a = _mm_xor_pd(a, _mm_set_pd(-0.0, 0.0));
    a = myshuffle(a, a, _mm_movemask_pd(b));
	a = _mm_div_pd(a, b);
    return a;
}

// mulpn: multiply positive-negative
// helpful for multiplication of two intervals that
// are known to be positive
static inline
__m128d mulpn(const __m128d a, const __m128d b)
{
	return _mm_mul_pd(a, _mm_xor_pd(b, _mm_set_pd(-0.0, 0.0)));
}

// ivsub: interval subtraction
// just a fusion of swapping and addition
// implementing a+(-b) for intervals
static inline
__m128d ivsub(__m128d a, __m128d b)
{
	return _mm_add_pd(a, _mm_shuffle_pd(b, b, 1));
}

// IntervalMul: interval multiplication
// using the ivchoice operation
static inline
__m128d IntervalMul(__m128d x, __m128d y) {
	__m128d a, b;
	a = _mm_shuffle_pd(x, x, 1);
	b = _mm_shuffle_pd(y, y, 1);
	y = ivchoice(y, a);
	b = ivchoice(b, x);
	y = _mm_mul_pd(y, a);
	b = _mm_mul_pd(b, x);
	y = _mm_min_pd(b, y);
	return y;
}

// IntervalMul2: interval multiplication
// using the ivmul operation
static inline
__m128d IntervalMul2(__m128d x, __m128d y) {
	__m128d a, b;
	a = _mm_shuffle_pd(x, x, 1);
	b = _mm_shuffle_pd(y, y, 1);
	y = ivmul(y, a);
	b = ivmul(b, x);
	y = _mm_min_pd(b, y);
	return y;
}

// IntervalDiv: interval division
// using the ivchoice operation
__m128d IntervalDiv(__m128d y, __m128d x) {
	__m128d a, b;
	if (_mm_movemask_pd(x)==3)
		throw RealLib::PrecisionException("IntervalDiv");
	a = _mm_shuffle_pd(x, x, 1);
	b = _mm_shuffle_pd(y, y, 1);
	y = ivchoice(y, a);
	b = ivchoice(b, x);
	y = _mm_div_pd(y, a);
	b = _mm_div_pd(b, x);
	y = _mm_min_pd(y, b);
	return y;
}

// IntervalDiv2: interval division
// using the ivdiv operation
__m128d IntervalDiv2(__m128d y, __m128d x) {
	__m128d a, b;
	if (_mm_movemask_pd(x)==3)
		throw RealLib::PrecisionException("IntervalDiv2");
	a = _mm_shuffle_pd(x, x, 1);
	b = _mm_shuffle_pd(y, y, 1);
	y = ivdiv(y, a);
	b = ivdiv(b, x);
	y = _mm_min_pd(y, b);
	return y;
}

#endif // FILE
