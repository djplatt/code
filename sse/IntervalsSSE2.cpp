#include <cfloat>
#include <cassert>
#include <exception>
#include <string.h>
#include "IntervalsSSE2.h"

namespace RealLib {

// this repeats the code in "RealEstimate.cpp"
RealLibException::RealLibException(const char *wht) throw()
{ 
     if (wht) { 
        strncpy(m_what, wht, 127);
        m_what[127] = 0;
     } else m_what[0] = 0;
} 

__m128d MachineEstimate::signmask = _mm_set_pd(0.0, -1.0 * 0.0);
__m128d MachineEstimate::mdelta = _mm_set1_pd(-DBL_MIN);
__m128d MachineEstimate::half = _mm_set1_pd(0.5);
__m128d MachineEstimate::mhalf = _mm_set_pd(-0.5, 0.5);
__m128d MachineEstimate::zero = _mm_setzero_pd();
__m128d MachineEstimate::mone = _mm_set1_pd(-1.0);
__m128d MachineEstimate::sqrt_corr;
__m128d MachineEstimate::sign = _mm_set1_pd(-0.0);

// how do we do initialization and finalization,
// i.e. all rounding mode changes, so that we can 
// be fine in a multithreaded program?

// problems:
//		-- a nested initialization should not simply clear the 
//			SSE exception flags. It needs to check them first
//		-- cannot use global/class-static variables, because they're
//			shared among threads, while the processor state isn't
//		-- when a computation is initialized the calling function does 
//			not know the context (higher-level comp object)

// a very big solution would be to record thread id's and save states
// according to them

// we prefer to keep it simple, and just use the SSE2 rounding mode as an
// indication whether or not this is a nested call, and we're just saving 
// what the rounding mode was to restore it after the call.

// this has a few restrictions:
//		-- other SSE-2 code has to restore rounding-to-nearest before
//		   our initialization gets called
//		-- wrong destruction order (not LIFO) will cause trouble


int MachineEstimate::BeginComputation()
{
	static bool initialized = false;

	if (!initialized) {
		double z = 0.0;
		double minusinf = -1.0 / z;
		sqrt_corr = //_mm_mul_pd(
			_mm_set_pd(0.0, minusinf);//, zero);
		initialized = true;
	}
	int rm = _MM_GET_ROUNDING_MODE();

#ifdef REALLIB_RELY_ON_SSE_EXCEPTIONS
	// if the initialization is already done at some higher level,
	// we need to check if an exception has occured already.
	if ((rm == _MM_ROUND_DOWN) && 
		(_MM_GET_EXCEPTION_STATE() & 
		(_MM_EXCEPT_INVALID | _MM_EXCEPT_DIV_ZERO | _MM_EXCEPT_OVERFLOW)))
		throw PrecisionException("SSE exception");
	_MM_SET_EXCEPTION_STATE(0);
#endif

	if (rm != _MM_ROUND_DOWN) 
		_MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);

	return rm;
}

void MachineEstimate::FinishComputation(int rm)
{
	assert(_MM_GET_ROUNDING_MODE()==_MM_ROUND_DOWN);
	_MM_SET_ROUNDING_MODE(rm);

#ifdef REALLIB_RELY_ON_SSE_EXCEPTIONS
	// throw exception if something has not been evaluated correctly
	// (only if we're not called by an exception)
	if (!std::uncaught_exception() &&
		(_MM_GET_EXCEPTION_STATE() & 
		(_MM_EXCEPT_INVALID | _MM_EXCEPT_DIV_ZERO | _MM_EXCEPT_OVERFLOW)))
		throw PrecisionException("SSE exception");
	_MM_SET_EXCEPTION_STATE(0);
#endif
}


}	// namespace

