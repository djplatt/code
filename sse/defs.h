#include "ieeefp.h"
#ifndef FILE_DEFS_H
#define FILE_DEFS_H

namespace RealLib {

typedef int i32;
typedef unsigned int u32;

#define I32_MAX INT_MAX
#define I32_MIN INT_MIN

#if defined(__GNUC__) || defined(__MWERKS__)
// doesn't define NDEBUG in non-debug builds
// comment this if you want a debug build
#define NDEBUG

typedef long long i64;
typedef unsigned long long u64;

#define NAMESPACE_STD std

#else	// MS Visual C++
typedef __int64 i64;
typedef unsigned __int64 u64;

#define NAMESPACE_STD 
#endif

#if defined(__MWERKS__)
typedef i32 exp_type;
#define MINIMUM_EXPONENT (-(1<<28))
#else
typedef i64 exp_type;
#define MINIMUM_EXPONENT (-I32_MAX)
#endif

#define I32_MIN INT_MIN		// to be used as -inf
#define I32_MAX INT_MAX		// to be used as +inf

} // namespace

#if defined(__GNUC__) || defined(__MWERKS__)
#include "GCChelper.h"
#endif

#if defined(__MWERKS__)
#define FUNCTION(x) &x
#else
#define FUNCTION(x) &(x)
#endif

// if this is set, multiplications of large number will
// be done via faster convolution
// recommended to keep this on
#define MULTIPLY_BY_CONVOLUTION

#ifdef MULTIPLY_BY_CONVOLUTION
// direct multiplication will be used for precisions below
// this threshold and convolution for larger precisions
#define CONVOLUTION_THRESHOLD 60
#else
#define CONVOLUTION_THRESHOLD INT_MAX
#endif

// without this the system will not limit its recursion depth
// may run slightly faster but will probably cause errors
// for longer computations on the class Real
// recommended to keep this on
#define LIMIT_RECURSION

#ifdef LIMIT_RECURSION
// the size of the chunk that is processed recursively
#define EvaluationDepth 500
#else
#define EvaluationDepth INT_MAX
#endif

#define LOG_2_10 3.32192809488736234787031942948939017586483139

#if !defined(NDEBUG) && defined(_MSC_VER)
#define CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif

#ifndef NDEBUG
// with this the system prints a message every time the
// precision is increased
#define REALLIB_SHOW_PRECISION_INCREASES
#endif

#if (_MSC_VER >= 1400)
#pragma warning (disable: 4996)
#endif

// the SSE2 version of MachineEstimate can rely on SSE2
// exception flags to recognize invalid operations
// and overflows
// with this turned off it will do more checks instead
// (slower performance)
#define REALLIB_RELY_ON_SSE_EXCEPTIONS

#ifndef NDEBUG
// relying on SSE exceptions does less checks
// and gives somewhat less precise information
// about the origin of the exception
// so use the other mode for debug builds

#undef REALLIB_RELY_ON_SSE_EXCEPTIONS
#endif

#ifdef REALLIB_RELY_ON_SSE_EXCEPTIONS
#define REALLIB_MACHEST_INVALID_CHECK(x) false
#else
#define REALLIB_MACHEST_INVALID_CHECK(x) !(x).IsValueValid()
#endif

#endif // file
