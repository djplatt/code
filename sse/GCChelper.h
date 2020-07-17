#ifndef FILE_GCCHELPER_H
#define FILE_GCCHELPER_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

namespace RealLib {

/* Return values for fpclass. */
#define	_FPCLASS_SNAN	0x0001	/* Signaling "Not a Number" */
#define	_FPCLASS_QNAN	0x0002	/* Quiet "Not a Number" */
#define	_FPCLASS_NINF	0x0004	/* Negative Infinity */
#define	_FPCLASS_NN	0x0008	/* Negative Normal */
#define	_FPCLASS_ND	0x0010	/* Negative Denormal */
#define	_FPCLASS_NZ	0x0020	/* Negative Zero */
#define	_FPCLASS_PZ	0x0040	/* Positive Zero */
#define	_FPCLASS_PD	0x0080	/* Positive Denormal */
#define	_FPCLASS_PN	0x0100	/* Positive Normal */
#define	_FPCLASS_PINF	0x0200	/* Positive Infinity */

static inline bool __isnan(const double x)
{
	return isnan(x);
	//return x != x;
}

static inline bool __isinf(const double x)
{
	return isinf(x);
	//return x == x + 1.0;
}

static inline bool _finite(const double x)
{
	return finite(x);
	//return x != x + 0.1;
}

static inline int _fpclass(const double x)
{
	if (__isnan(x)) return _FPCLASS_QNAN;
	else if (__isinf(x)) return x > 0 ? _FPCLASS_PINF : _FPCLASS_NINF;
	else if (x == 0.0) return _FPCLASS_PZ;
	else return x > 0 ? _FPCLASS_PN : _FPCLASS_NN;
}

} // namespace

#endif
