//
// File int_double.h
//
// Define interval arithmetic on doubles
// extend to complexes using rectangles
//
//
// Use SSE op-codes for IEEE compliance and
// efficiency. 
// See "Interval Arithmetic using SSE-2" Branimir Lambov.
// In: Hertling P., Hoffmann C.M., Luther W., Revol N. (eds) 
//   Reliable Implementation of Real Number Algorithms: Theory and Practice. 
//   Lecture Notes in Computer Science, vol 5045. Springer, Berlin, Heidelberg
//
// 22/3/9 v 8.0 added use of nextbefore
// 15/4/9 v 9.0 added SSE assembler for GCC
// 30/4/9 v 10.0 added SSE for GCC + Intel Linux
// 24/7/9 v 11.0 coded log and exp in FP assembler
// 27/6/10 v 12.0 replaced FPU calls to sin, cos, atan, exp, log
//                with calls to crlibm functions. Can't trust FPU!
//                No longer mess with FPU rounding so crlibm has a
//                chance.
// 8/5/2017 v 13.0 Added logp1(x)=log(1+x)
//          v 14.0 Fixed abs(int_double)
//          v 14.1 Changed print_int_double_str in operator /
//                 to remove empty string
//          v 14.2 Changed constants to use uint32_t
//
//
// To use:-
// IMPORTANT. THIS CODE WORKS BY CHANGING THE DEFAULT ROUNDING MODE ON THE SSE
// REGISTERS TO ROUND DOWN. YOU MUST ENSURE THAT NOTHING ELSE IS RELYING ON
// SSE REGISTERS BEHAVING NORMALLY. IN PARTICULAR, WHEN BUILDING CRLIBM, 
// USE APPROPRIATE COMPILER FLAGS TO ENSURE NO SIMD INSTRUCTIONS ARE USED
// (UNLESS YOU HAPPEN TO KNOW HOW CRLIBM WILL BEHAVE IN ROUND DOWN MODE).
//
// Suggested gcc compile flags for INT_DOUBLE applications 
// -fomit-frame-pointer -O1 -msse3 -mfpmath=387 -frounding-math -finline-functions
//
// A note on the use of inline assembler
//
// We realise that using inline assembler makes the code less readable and,
// as it stands, such that it will only compile under gcc. We could have
// used compiler intrinsics to get more readable and portable code, but we
// wanted to ensure that the assembly output during compile was EXACTLY what
// we specified. We apologise if this offends anyone.

/*
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef INT_DOUBLE
#define INT_DOUBLE


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <inttypes.h>
#include "malloc.h"
#include "crlibm.h"

using namespace std;

// simple print routines. Would be nice to have a printf version
#define print_int_complex_str(str,x) {printf(str);printf(" ");print_int_complex(x);printf("\n");}
#define print_int_double_str(str,x) {printf(str);printf(" ");print_int_double(x);printf("\n");}

// nextafter(before) returns the next representable double after(before)
inline double nextafter(double);
inline double nextbefore(double);

// int_double must be 16 byte aligned for SSE instructions
// it is stored as [left,-right]
// e.g. the Interval [1,2] looks like [1.0,-2.0]
class int_double{
public:
  double __attribute__ ((aligned(16))) left;
  double right;
	
inline  int_double ()                      // constructors
  {
  };

// does no rounding, see below
inline int_double(double l) 
{
  left=l;
  right=-l;
};

// N.B. this default contructor does no rounding.
// Note that the compiler probably rounds to nearest
// e.g. third=1.0/3.0; int_third=int_double(third,third);
//      results in a point interval at nearest to 1/3.
//
// you probably want int_third=int_double(1.0)/3.0;
inline int_double(double l,double r)
{
  if(l>r)
    {
      printf("Error constructing int_double, right %20.18e < left %20.18e . Exiting.\n",r,l);
      exit(1);
    };
	
  left=l;
  right=-r;
 }

 friend int_double operator + (const int_double &lhs, const int_double &rhs);
 friend int_double operator + (const int_double &lhs, const double &rhs);
 friend int_double operator + (const double &lhs, const int_double &rhs);
 friend int_double operator - (const int_double &lhs, const int_double &rhs);
 friend int_double operator - (const int_double &lhs, const double &rhs);
 friend int_double operator - (const int_double &lhs);
 friend int_double operator * (const int_double &lhs, const int_double &rhs);
 friend int_double operator * (const int_double &lhs, const double &rhs);
 friend int_double operator * (const double &lhs, const int_double &rhs);
 friend int_double operator / (const int_double &lhs, const int_double &rhs);
 friend int_double operator / (const int_double &lhs, const double &rhs);
 friend int_double operator += (int_double &lhs, const int_double &rhs);
 friend int_double operator += (int_double &lhs, const double &rhs);
 friend int_double operator -= (int_double &lhs, const int_double &rhs);
 friend int_double operator -= (int_double &lhs, const double &rhs);
 friend int_double operator *= (int_double &lhs, const int_double &rhs);
 friend int_double operator *= (int_double &lhs, const double &rhs);
 friend int_double operator /= (int_double &lhs, const int_double &rhs);
 friend int_double operator /= (int_double &lhs, const double &rhs);
 friend int operator >= (const int_double &lhs, const int_double &rhs);
 friend int operator > (const int_double &lhs, const int_double &rhs);
 friend int operator < (const int_double &lhs, const int_double &rhs);
}; // end class

void print_int_double(const int_double &x)
{
  printf("[ %20.18e , %20.18e ]",x.left,-x.right);
};

double width(const int_double x)
{
  return(-x.right-x.left);
}


int_double times_pos(const int_double &, const int_double &);

int_double exp (const int_double &);

int_double logp1 (const int_double &);

int_double log (const int_double &);

void sin_cos(const int_double &, int_double *, int_double *);

void sin_cospi(const int_double &, int_double *, int_double *);

int_double sqr(const int_double &);

int_double pow(const int_double &, const int_double &);

int_double pow(const int_double &, const double &);

int_double atan2(const int_double &, const int_double &);

int_double sqrt(const int_double &);

int_double cosh(const int_double &);

int rel_error(const int_double &);

int abs_error(const int_double &);

int contains_zero (const int_double &);

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Define int_complex */
/* operators + += - -= * *= / defined for 
   int_complex op int_complex, int_double, double
   unary - defined for int_complex
   */
class int_complex{
 public:
  int_double real;
  int_double imag;

  inline int_complex ()
    {};
  
  inline int_complex(const int_double &re)
    {
      real=re;
      imag=int_double(0.0);
    }

  inline int_complex (const int_double &re,const int_double &im)
    {
      real=re;
      imag=im;
    };


  friend int_complex operator + (const int_complex &lhs, const int_complex &rhs);
  friend int_complex operator + (const int_complex &lhs, const int_double &rhs);
  friend int_complex operator + (const int_double &lhs, const int_complex &rhs);
  friend int_complex operator + (const int_complex &lhs, const double &rhs);
  friend int_complex operator + (const double &lhs, const int_complex &rhs);

  friend int_complex operator * (const int_complex &lhs, const int_complex &rhs);
  friend int_complex operator * (const int_complex &lhs, const int_double &rhs);
  friend int_complex operator * (const int_double &lhs, const int_complex &rhs);
  friend int_complex operator * (const int_complex &lhs, const double &rhs);
  friend int_complex operator * (const double &lhs, const int_complex &rhs);

  friend int_complex operator - (const int_complex &lhs, const int_complex &rhs);
  friend int_complex operator - (const int_complex &lhs, const int_double &rhs);
  friend int_complex operator - (const int_complex &lhs, const double &rhs);
  friend int_complex operator - (const int_complex &lhs);

  friend int_complex operator / (const int_complex &lhs, const int_complex &rhs);
  friend int_complex operator / (const int_complex &lhs, const int_double &rhs);
  friend int_complex operator / (const int_complex &lhs, const double &rhs);

  friend int_complex operator += (int_complex &lhs, const int_complex &rhs);
  friend int_complex operator += (int_complex &lhs, const int_double &rhs);
  friend int_complex operator += (int_complex &lhs, const double &rhs);
  friend int_complex operator -= (int_complex &lhs, const int_complex &rhs);
  friend int_complex operator -= (int_complex &lhs, const int_double &rhs);
  friend int_complex operator -= (int_complex &lhs, const double &rhs);
  friend int_complex operator *= (int_complex &lhs, const int_complex &rhs);
  friend int_complex operator *= (int_complex &lhs, const int_double &rhs);
  friend int_complex operator *= (int_complex &lhs, const double &rhs);
  friend int_complex operator /= (int_complex &lhs, const int_complex &rhs);
  friend int_complex operator /= (int_complex &lhs, const int_double &rhs);
  friend int_complex operator /= (int_complex &lhs, const double &rhs);
};

int_complex exp(const int_complex &);

int_complex e(const int_double &);

void print_int_complex(const int_complex &);

int_complex pow (const int_double &,const int_complex &);

int_complex pow1(const double &, const double &, const double &);

int_complex pow1(const int_double &, const double &, const double &);

int_complex pow (const double, const int_complex &);

int_double norm(const int_complex &);

int_complex sqrt(const int_complex &);

int_complex conj (const int_complex &);

int_double argument (const int_complex &);

int_complex log (const int_complex &);

int contains_zero (const int_complex &);

inline int_double real (const int_complex &z)
{
	return(z.real);
}

inline int_double imag (const int_complex &z)
{
	return(z.imag);
}

// Useful constants
//
// a half
const int_double d_half=int_double(0.5,0.5);
//
// Pi
uint32_t _i_pi[2]={1413754136,1074340347};   // this is what double pi less a bit looks like
uint32_t _i_pi2[2]={1413754137,1074340347};  // this is double pi plus a bit
double *_d_pi=(double *)&_i_pi;
double *_d_pi2=(double *)&_i_pi2;
int_double d_pi=int_double(_d_pi[0],_d_pi2[0]);
//
// Euler gamma
uint32_t _i_gamma[2]={0xfc6fb618,0x3fe2788c};
uint32_t _i_gamma2[2]={0xfc6fb619,0x3fe2788c};
double *_d_gamma=(double *)&_i_gamma;
double *_d_gamma2=(double *)&_i_gamma2;
int_double d_gamma=int_double(_d_gamma[0],_d_gamma2[0]);
//
// Pi/2
uint32_t _i_pi_2[2]={0x54442d18,0x3ff921fb};   // this is what double pi/2 less a bit looks like
uint32_t _i_pi2_2[2]={0x54442d19,0x3ff921fb};  // this is double pi/2 plus a bit
double *_d_pi_2=(double *)&_i_pi_2;
double *_d_pi2_2=(double *)&_i_pi2_2;
int_double d_pi_2=int_double(_d_pi_2[0],_d_pi2_2[0]);
//
// 2*Pi
uint32_t _i_2pi_2[2]={0x54442d18,0x401921fb};   // this is what double pi*2 less a bit looks like
uint32_t _i_2pi2_2[2]={0x54442d19,0x401921fb};  // this is double pi*2 plus a bit
double *_d_2pi_2=(double *)&_i_2pi_2;
double *_d_2pi2_2=(double *)&_i_2pi2_2;
int_double d_two_pi=int_double(_d_2pi_2[0],_d_2pi2_2[0]);
//
// log(pi)
int_double d_ln_pi;
//
// log(2*pi)
int_double d_ln_two_pi;
//
// The following are used as masks in SSE code
//
//
uint32_t __nze[2]={0,0x80000000}; // -0.0
double _nze = *(double*)__nze;
int_double d_zero_zero;
int_double d_zero;  // +0.0 -0.0
int_double d_neg_zero; // -0.0 -0.0
int_double d_neg_neg_zero; // -0.0 +0.0

// but these should be ok
int_double d_one=int_double(1.0);
int_complex c_zero=int_complex(int_double(0.0),int_double(0.0));
int_complex c_one=int_complex(d_one,int_double(0.0));
int_complex c_half=int_complex(int_double(0.5),int_double(0.0));

// move a double prec float towards +infinity by
// the smallest delta possible
int_double delta_int_double=int_double(DBL_MIN);
int_double delta_int_double_neg=int_double(-DBL_MIN);

inline double nextafter (double x)
{
  int_double y;
  y.right=-x;
  int_double z=x+delta_int_double;
  return(-z.right);
}

inline double nextbefore (double x)
{
  int_double y;
  y.left=x;
  int_double z=x+delta_int_double_neg;
  return(z.left);
}

int_double delta_blow=int_double(-DBL_MIN,DBL_MIN);

inline int_double blow (int_double x)
{
  return(x+delta_blow);
}

// return (log of) relative error
inline int rel_error(const int_double &x)
{
  double rerr;
  if(x.left==0.0)
    return(0);
  else
    {
      rerr=fabs((x.left+x.right)/x.left);
      if(rerr==0.0)
	return(0);
      else
	return((int)floor(log10(rerr)));
    }
}

// return (log of) absolute error
inline int abs_error(const int_double &x)
{
  double aerr;
  aerr=x.left+x.right;
  if(aerr==0.0)
    return(0);
  if(aerr<0.0)
    return((int)floor(log10(-aerr)));
  return((int)floor(log10(aerr)));
}

inline int_double operator + (const int_double &lhs, const int_double &rhs)
{
  int_double temp;
  __asm__("movapd %1,%%XMM0\n\t" // lhs.left lhs.right
	  "addpd %2,%%XMM0\n\t"  // lhs.left+rhs.left lhs.right+rhs.right
	  "movapd %%XMM0,%0\n\t"
	  :"=m" (temp)
	  : "m" (lhs), "m" (rhs)
	  :"xmm0");
  return(temp);
}

inline int_double operator + (const int_double &lhs, const double &rhs)
{
	int_double temp;
	__asm("movddup %2,%%XMM0\n\t" // rhs rhs
	  "xorpd %3,%%XMM0\n\t"   // rhs -rhs
	  "addpd %1,%%XMM0\n\t"   // lhs.left+rhs rhs.right-rhs
	  "movapd %%XMM0,%0\n\t"
	  :"=m" (temp)
      :"m" (lhs), "m" (rhs), "m" (d_zero)
	  :"xmm0");
	return(temp);
}

inline int_double operator + (const double &rhs, const int_double &lhs)
{
	int_double temp;
	__asm("movddup %2,%%XMM0\n\t" // rhs rhs
	  "xorpd %3,%%XMM0\n\t"   // rhs -rhs
	  "addpd %1,%%XMM0\n\t"   // lhs.left+rhs rhs.right-rhs
	  "movapd %%XMM0,%0\n\t"
	  :"=m" (temp)
      :"m" (lhs), "m" (rhs), "m" (d_zero)
	  :"xmm0");
	return(temp);
}

inline int_double operator += (int_double &lhs, const int_double &rhs)
{
	__asm__("movapd %0,%%XMM0\n\t" // lhs.left lhs.right
		"addpd %1,%%XMM0\n\t"  // lhs.left+rhs.left lhs.right+rhs.right
		"movapd %%XMM0,%0\n\t"
		:
		: "m" (lhs), "m" (rhs)
		:"xmm0");
	return(lhs);
}

inline int_double operator += (int_double &lhs, const double &rhs)
{
	__asm("movddup %1,%%XMM0\n\t" // rhs rhs
	  "xorpd %2,%%XMM0\n\t"   // rhs -rhs
	  "addpd %0,%%XMM0\n\t"   // lhs.left+rhs rhs.right-rhs
	  "movapd %%XMM0,%0\n\t"
	  :
      :"m" (lhs), "m" (rhs), "m" (d_zero)
	  :"xmm0");
	return(lhs);
}

inline int_double operator - (const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	__asm("movapd %2,%%xmm0\n\t"
		"shufpd $1,%%xmm0,%%xmm0\n\t" // -rhs
		"addpd %1,%%xmm0\n\t"         // lhs-rhs
		"movapd %%xmm0,%0\n\t"
		  :"=m" (temp)
		  :"m" (lhs), "m" (rhs)
		  :"xmm0");
	return(temp);
}

inline int_double operator - (const int_double &lhs, const double &rhs)
{
	int_double temp;
	__asm("movddup %2,%%XMM0\n\t"       // rhs rhs
	      "movapd %1,%%XMM1\n\t"        // l.l l.r
	      "addsubpd %%XMM0,%%XMM1\n\t"  // l.l+rhs l.r-rhs
	      "movapd %%XMM1,%0"
	      :"=m" (temp)
	      :"m" (lhs),"m" (rhs)
	      :"xmm0", "xmm1");
	return(temp);
}

inline int_double operator - (const int_double &lhs)
  {
	  int_double temp;
	  __asm("movapd %1,%%XMM0\n\t"
		  "shufpd $1,%%XMM0,%%XMM0\n\t"	// -lhs
		  "movupd %%XMM0,%0\n\t"
		:"=m" (temp)
		:"m" (lhs)
		:"xmm0");
	  return(temp);
}

inline int_double operator -= (int_double &lhs, const int_double &rhs)
{
	return(lhs=lhs-rhs);
}

inline int_double operator -= (int_double &lhs, const double &rhs)
{
	return(lhs=lhs-rhs);
}

// Lambov's version
inline int_double operator * (const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	__asm(
	"movapd %2,%%xmm0\n\t"
	"xorpd	%3,%%xmm0\n\t"			// a=xmm0
	"movapd %1,%%xmm1\n\t"		// b=xmm1
	"movapd %%xmm1,%%xmm2\n\t"
	"shufpd	$1,%%xmm2,%%xmm2\n\t"		// p=xmm2
	"movapd %%xmm0,%%xmm3\n\t"
	"cmpltpd %4,%%xmm3\n\t"			// c=xmm3
	"xorpd	%4,%%xmm2\n\t"			// p=? d=xmm2
	"movapd %%xmm3,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=xmm4
	"movapd %%xmm3,%%xmm5\n\t"
	"andpd %%xmm2,%%xmm5\n\t"		// e=xmm5
	"andnpd %%xmm1,%%xmm3\n\t"		// c=? f=xmm3
	"movapd %%xmm4,%%xmm6\n\t"
	"andnpd %%xmm1,%%xmm6\n\t"		// k=xmm6
	"andpd %%xmm4,%%xmm2\n\t"		// d=? i=xmm2
	"orpd %%xmm5,%%xmm3\n\t"		// f=? h=xmm3
	"movapd %%xmm0,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=? l=xmm4
	"orpd %%xmm6,%%xmm2\n\t"		// i=? m=xmm2
	"mulpd %%xmm3,%%xmm0\n\t"		// a=? j=xmm0
	"mulpd %%xmm4,%%xmm2\n\t"		// m=? n=xmm2
	"minpd %%xmm2,%%xmm0\n\t"		// j=? o=xmm0
	"movapd %%xmm0,%0\n\t"		// return(o)
	:"=m" (temp)
	:"m" (lhs), "m" (rhs), "m" (d_zero), "m" (d_neg_zero)
	:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6");
	return(temp);
}

// faster version of times if we know everything is +ve
inline int_double times_pos(const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	__asm("movapd %2,%%XMM0\n\t" // rhs.left rhs.right
		"xorpd %3,%%XMM0\n\t"  // rhs.left -rhs.right
		"mulpd %1,%%XMM0\n\t"  // lhs.left*rhs.left lhs.right*(-rhs.right)
		"movapd %%XMM0,%0\n\t"
		:"=m" (temp)
		:"m" (lhs), "m" (rhs), "m" (d_zero)
		:"xmm0");
	return(temp);
}

inline int_double operator *= (int_double &lhs, const int_double &rhs)
{
	return(lhs=lhs*rhs);
}

inline int_double operator * (const int_double &lhs, const double &rhs)
{
	int_double temp;
	__asm("movddup %2,%%xmm0\n\t"         // rhs rhs
		"movapd %%xmm0,%%xmm1\n\t"    // rhs rhs
		"xorpd %3,%%xmm1\n\t"         // -rhs -rhs
		"movapd %1,%%xmm2\n\t"        // lhs.left lhs.right
		"mulpd %%xmm2,%%xmm0\n\t"     // lhs.left*rhs lhs.right*rhs
		"mulpd %%xmm2,%%xmm1\n\t"     // lhs.left*(-rhs) lhs.right*(-rhs)
		"shufpd $1,%%xmm1,%%xmm1\n\t"  // lhs.right*(-rhs) lhs.left*(-rhs)
		"minpd %%xmm1,%%xmm0\n\t"     // min(lhs.left*rhs,lhs.right*(-rhs)) min(lhs.right*rhs,lhs.left*(-rhs))
		"movapd %%xmm0,%0"            // move to temp
		:"=m" (temp)
		:"m" (lhs), "m" (rhs), "m" (d_neg_zero)
		:"xmm0", "xmm1", "xmm2");
	return(temp);
}

inline int_double operator * (const double &rhs, const int_double &lhs)
{
	int_double temp;
	__asm("movddup %2,%%xmm0\n\t"         // rhs rhs
		"movapd %%xmm0,%%xmm1\n\t"    // rhs rhs
		"xorpd %3,%%xmm1\n\t"         // -rhs -rhs
		"movapd %1,%%xmm2\n\t"        // lhs.left lhs.right
		"mulpd %%xmm2,%%xmm0\n\t"     // lhs.left*rhs lhs.right*rhs
		"mulpd %%xmm2,%%xmm1\n\t"     // lhs.left*(-rhs) lhs.right*(-rhs)
		"shufpd $1,%%xmm1,%%xmm1\n\t"  // lhs.right*(-rhs) lhs.left*(-rhs)
		"minpd %%xmm1,%%xmm0\n\t"     // min(lhs.left*rhs,lhs.right*(-rhs)) min(lhs.right*rhs,lhs.left*(-rhs))
		"movapd %%xmm0,%0"            // move to temp
		:"=m" (temp)
		:"m" (lhs), "m" (rhs), "m" (d_neg_zero)
		:"xmm0", "xmm1", "xmm2");
	return(temp);
}

inline int_double operator *= (int_double &lhs, const double &rhs)
{
	return(lhs=lhs*rhs);
}

inline int_double operator / (const int_double &lhs, const int_double &rhs)
{
  int_double temp;

  // if you are sure this can't happen, then remove this test.		
  if(contains_zero(rhs))
    {
      print_int_double_str("Division by interval containing zero. Exiting.\n Interval =  ",rhs);
      exit(1);
    }
  
  __asm("movapd %2,%%xmm0\n\t"
	"xorpd %3,%%xmm0\n\t"
	"movapd %%xmm0,%%xmm1\n\t"
	"shufpd $1,%%xmm1,%%xmm1\n\t"
	"movapd %%xmm0,%%xmm2\n\t"
	"cmpnlepd %4,%%xmm2\n\t"
	"movapd %1,%%xmm3\n\t"
	"movapd %%xmm3,%%xmm4\n\t"
	"shufpd $1,%%xmm3,%%xmm3\n\t"
	"xorpd %5,%%xmm3\n\t"
	"andpd %%xmm2,%%xmm4\n\t"
	"andnpd %%xmm3,%%xmm2\n\t"
	"orpd %%xmm4,%%xmm2\n\t"
	"movapd %%xmm2,%%xmm3\n\t"
	"divpd %%xmm1,%%xmm2\n\t"
	"divpd %%xmm0,%%xmm3\n\t"
	"minpd %%xmm2,%%xmm3\n\t"
	"movapd %%xmm3,%0"
	:"=m" (temp)
	:"m" (lhs), "m" (rhs), "m" (d_zero), "m" (d_zero_zero), "m" (d_neg_zero)
	:"xmm0", "xmm1", "xmm2", "xmm3", "xmm4");
  return(temp);
}

inline int_double operator /= (int_double &lhs, const int_double &rhs)
{
	return(lhs=lhs/rhs);
}

inline int_double operator / (const int_double &lhs, const double &rhs)
{
  double t1;
  int_double temp;
  if(rhs>0.0)
    {
      __asm("movddup %2,%%xmm0\n\t"
	    "movapd %1,%%xmm1\n\t"
	    "divpd %%xmm0,%%xmm1\n\t"
	    "movapd %%xmm1,%0"
	    :"=m" (temp)
	    :"m" (lhs), "m" (rhs)
	    :"xmm0", "xmm1");
      return(temp);
    }
  // if you are sure rhs cannot = 0, then this test is nugatory
  if(rhs<0.0)
    {
      t1=-rhs;
      __asm("movddup %2,%%xmm0\n\t"		// -rhs -rhs
	    "movapd %1,%%xmm1\n\t"		// lhs.left lhs.right
	    "shufpd $1,%%xmm1,%%xmm1\n\t"	// lhs.right lhs.left
	    "divpd %%xmm0,%%xmm1\n\t"	// lhs.right/-rhs lhs.left/-rhs
	    "movapd %%xmm1,%0"
	    :"=m" (temp)
	    :"m" (lhs), "m" (t1)
	    :"xmm0", "xmm1");
      return(temp);
    }
  printf("Division by zero in int_double / double. Exiting.\n");
  exit(1);
}

inline int_double operator /= (int_double &lhs, const double &rhs)
{
	return(lhs=lhs/rhs);
}

inline int operator >= (const int_double &lhs, const int_double &rhs)
{
	return(lhs.left>=(-rhs.right));
}

inline int operator > (const int_double &lhs, const int_double &rhs)
{
	return(lhs.left>(-rhs.right));
}

inline int operator < (const int_double &lhs, const int_double &rhs)
{
	return((-lhs.right)<rhs.left);
}


//
// all trancendental routines use crlibm
// which must be built to use FPU exclusively
//
// can't trust C++ exp routine if we change rounding
// Probably not needed as we don't touch FPU
// here just in case 
inline double exp1(const double x)
{
  return(exp_rn(x));
}

inline int_double exp (const int_double &x)  // nicely increasing
{
  return(int_double(exp_rd(x.left),exp_ru(-x.right)));
}

// if right hand endpoint is one, then interval is
// widened unnecessarily. could trap it but......
inline int_double log (const int_double &x) // nicely increasing
{
  return(int_double(log_rd(x.left),log_ru(-x.right)));
}

inline int_double log1p (const int_double &x) // nicely increasing
{
  return(int_double(log1p_rd(x.left),log1p_ru(-x.right)));
}

// masks to show which quadrant our theta is in
#define Q1 (0b0000) // 1st
#define Q2 (0b1100) // 2nd
#define Q3 (0b1111)
#define Q4 (0b0011)
#define Q14 (0b0010) // straddling 1 and 4
#define Q12 (0b1000)
#define Q23 (0b1110)
#define Q34 (0b1011)
#define QZ (0b1010) // stradding zero
//in the mask, a 1 indicates negative

inline int make_mask1(double a, double b, double c, double d)
{
  int mask=0;
  if(a<0.0)
    mask|=0b1000;
  if(b<0.0)
    mask|=0b0100;
  if(c<0.0)
    mask|=0b0010;
  if(d<0.0)
    mask|=0b0001;
  return(mask);
}

inline int make_mask(const int_double &x, const int_double &y)
{
  return(make_mask1(x.left, -x.right, y.left, -y.right));
}

inline void sin_cos1(const double &cos_left, const double &cos_right,
		     const double &sin_left, const double &sin_right,
		     int_double *sin_x, int_double *cos_x)
{
  unsigned char mask=0;
		     
  if(cos_left>=0.0)
    mask=8;
  if(cos_right>=0.0)
    mask=mask|0x4;
  if(sin_left>=0.0)
    mask=mask|0x2;
  if(sin_right>=0.0)
    mask++;
  
  switch (mask)
    {
    case 0: // all -ve, 3rd quadrant
      sin_x->left=sin_right;
      sin_x->right=-nextafter(sin_left);
      cos_x->left=cos_left;
      cos_x->right=-nextafter(cos_right);
      return;
    case 2:
      sin_x->left=nextbefore(sin_right);  // sin_right < 0 and too close to 0
      sin_x->right=-nextafter(sin_left);
      cos_x->left=-1.0;
      if(cos_left>=cos_right)
	cos_x->right=-nextafter(cos_left);
      else
	cos_x->right=-nextafter(cos_right);
      return;
    case 3:
      sin_x->left=sin_right;
      sin_x->right=-nextafter(sin_left);
      cos_x->left=cos_right;
      cos_x->right=-nextafter(cos_left);
      return;
    case 4:
      sin_x->left=-1.0;
      if(sin_right>=sin_left)
	sin_x->right=-nextafter(sin_right);
      else
	sin_x->right=-nextafter(sin_left);
      cos_x->left=nextbefore(cos_left);
      cos_x->right=-nextafter(cos_right);
      return;
    case 11:
      if(sin_left<=sin_right)
	sin_x->left=sin_left;
      else
	sin_x->left=sin_right;
      sin_x->right=-1.0;
      cos_x->left=cos_right;
      cos_x->right=-nextafter(cos_left);
      return;
    case 12:
      sin_x->left=sin_left;
      sin_x->right=-nextafter(sin_right);
      cos_x->left=cos_left;
      cos_x->right=-nextafter(cos_right);
      return;
    case 13:
      cos_x->right=-1.0;
      if(cos_left<=cos_right)
	cos_x->left=cos_left;
      else
	cos_x->left=cos_right;
      sin_x->left=sin_left;
      sin_x->right=-nextafter(sin_right);
      return;
    case 15:
      sin_x->left=sin_left;
      sin_x->right=-nextafter(sin_right);
      cos_x->left=cos_right;
      cos_x->right=-nextafter(cos_left);
      return;
    default:
      printf("Weird error in sin_cos, mask = %d Exiting.\n",(int) mask);
      printf("\ncos_left: %20.18e\ncos_right: %20.18e\n",cos_left,cos_right);
      printf("sin_left: %20.18e\nsin_right: %20.18e\n",sin_left,sin_right);
      exit(1);
    }
}

// perform sin(x) into sin_x and cos(x) into cos_x
// if you want to handle sin/cos returning [-1,1] then define FULL_SIN_COS
inline void sin_cos(const int_double &x, int_double *sin_x, int_double *cos_x)
{
  double cos_left,cos_right,sin_left,sin_right;
  unsigned char mask=0;

  if(x.right+x.left<=-M_PI_2)
    {
#ifndef FULL_SIN_COS
      printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");
      exit(1);
#endif
      sin_x->left=-1;sin_x->right=-1;
      cos_x->left=-1;cos_x->right=-1;
      return;
    }

  int_double xp=x/d_pi;
  cos_left=cospi_rd(xp.left);
  cos_right=cospi_rd(-xp.right); // - probably a nop because cos is even
  sin_left=sinpi_rd(xp.left);
  sin_right=sinpi_rd(-xp.right);
  sin_cos1(cos_left,cos_right,sin_left,sin_right,sin_x,cos_x);
}

//
// sin (pi*x), cos(pi*x)
//
inline void sin_cospi(const int_double &x, int_double *sin_x, int_double *cos_x)
{
  double cos_left,cos_right,sin_left,sin_right;

  if(x.right+x.left<=-0.5)
    {
#ifndef FULL_SIN_COS
      printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");
      exit(1);
#endif
      sin_x->left=-1;sin_x->right=-1;
      cos_x->left=-1;cos_x->right=-1;
      return;

    }

  cos_left=cospi_rd(x.left);
  cos_right=cospi_rd(-x.right); // probably a nop because cos is even
  sin_left=sinpi_rd(x.left);
  sin_right=sinpi_rd(-x.right);
  sin_cos1(cos_left,cos_right,sin_left,sin_right,sin_x,cos_x);
}

// return sinc(x)=sin(x)/x
inline int_double sinc(const int_double &x)
{
  int_double s,c;
  if(contains_zero(x))
    {
      printf("Haven't written sinc to handle intervals containing zero. Exiting.");
      exit(1);
    }
  sin_cos(x,&s,&c);
  return(s/x);
}

// return sinc(x*Pi)
inline int_double sincpi(const int_double &x)
{
  int_double s,c;
   if(contains_zero(x))
    {
      printf("Haven't written sincpi to handle intervals containing zero. Exiting.");
      exit(1);
    }
   sin_cospi(x,&s,&c);
   return(s/(x*d_pi));
}


// x^y
inline int_double pow(const int_double &x, const int_double &y)
{
	return(exp(log(x)*y));
}

inline int_double pow(const int_double &x, const double &y)
{
	return(exp(log(x)*y));
}

// simple atan
inline int_double atan(const int_double &x)
{
  return(int_double(atan_rd(x.left),atan_ru(-x.right)));
}

// returns theta in [-Pi,Pi]
inline int_double atan2(const int_double &y, const int_double &x)
{
  //print_int_double_str("x=",x);
  //print_int_double_str("y=",y);

  int mask=make_mask(x,y);

  switch(mask)
    {
    case Q1:
    case Q4:
    case Q14:
      return(atan(y/x));
    case Q12: // rotate by -Pi/2
    case Q2:
      return(atan((-x)/y)+d_pi_2);
    case Q34: // rotate by Pi/2
    case Q3:
      return(atan(x/(-y))-d_pi_2);
    case Q23:
      printf("Error in atan2. Can't handle arguments on negative real axis. Exiting.\n");
      exit(1);
    case QZ:
      printf("Error in atan2. Both x and y contain zero. Exiting.\n");
      exit(1);
    default:
      printf("Bad mask in atan2. Mask=%d. Exiting.\n",mask);
      print_int_double_str("x=",x);
      print_int_double_str("y=",y);
      exit(1);
    }
}

inline int_double sqr(const int_double &x)
{
	int_double res;
	__asm("movapd %1,%%xmm0\n\t"
		"movapd %%xmm0,%%xmm1\n\t"
		"movapd %%xmm0,%%xmm2\n\t"
		"shufpd $1,%%xmm1,%%xmm1\n\t"
		"minsd %%xmm1,%%xmm0\n\t"
		"maxsd %%xmm2,%%xmm1\n\t"
		"maxsd %2,%%xmm1\n\t"
		"unpcklpd %%xmm0,%%xmm1\n\t"
		"movapd %%xmm1,%%xmm0\n\t"
		"xorpd %3,%%xmm1\n\t"
		"mulpd %%xmm1,%%xmm0\n\t"
		"movapd %%xmm0,%0"
		:"=m" (res)
		:"m" (x), "m" (d_zero_zero), "m" (d_zero)
		:"xmm0","xmm1","xmm2");
	return(res);
}

// SSE sqrt is IEEE compliant
// so we could use it, but we don't
inline int_double sqrt(const int_double &x) // strictly increasing
{
  double l=x.left,r=x.right,lhs,rhs;
  __asm("fldl %2\n\t"
	"fsqrt\n\t"
	"fstpl %0\n\t"
	"fldl %3\n\t"
	"fchs\n\t"
	"fsqrt\n\t"
	"fstpl %1\n\t"
	: "=m" (lhs), "=m" (rhs)
	:"m" (l), "m" (r)
	: "st(7)" );
  return(blow(int_double(lhs,rhs)));//int_double(nextbefore(lhs),nextafter(rhs)));
}

inline int_double abs(const int_double &x)
{
  if(x.left>=0.0) return x;
  if(x.right>=0.0) return -x;
  int_double y;
  y.left=0.0;
  y.right=min(x.left,x.right);
  return y;
}

inline int_double cosh(const int_double &x)
{
  return(0.5*(exp(x)+exp(-x)));
}

inline int contains_zero(const int_double &x)
{
  int temp;
  __asm("movapd %1,%%xmm0\n\t"
	"cmplepd %2,%%xmm0\n\t"		// cmp with 0.0,0.0 
	"movmskpd %%xmm0,%%eax\n\t" // = 0b11 iff x contains zero
	"movl %%eax,%%ecx\n\t"
	"shr %%eax\n\t"				// = 0b01 iff lhs.left<=0.0
	"and %%ecx,%%eax\n\t"		// eax = 1 iff left<=0 and right <=0, else 0
	"movl %%eax,%0"
	:"=m" (temp)
	:"m" (x), "m" (d_zero_zero)
	:"eax", "ecx", "xmm0");
  return(temp);
}

inline int contains_zero(const int_complex &z)
{
	return(contains_zero(z.real)&&contains_zero(z.imag));
}

inline double avg_int(int_double x)
{
	return((x.left-x.right)/2.0);
}


inline int_complex conj(const int_complex &rhs)
{
  int_complex temp;
  __asm("movapd %2,%%xmm0\n\t"
	"movapd %%xmm0,%0\n\t"
	"movapd %3,%%xmm1\n\t"
	"shufpd $1,%%xmm1,%%xmm1\n\t"
	"movapd %%xmm1,%1\n\t"
	:"=m" (temp.real), "=m" (temp.imag)
    	:"m" (rhs.real), "m" (rhs.imag)
	: "xmm0", "xmm1");
  return(temp);
}

inline int_complex operator + (const int_complex &lhs, const int_complex &rhs)
{
  int_complex temp;
  __asm("movapd %2,%%xmm0\n\t"
	"addpd %4,%%xmm0\n\t"
	"movapd %%xmm0,%0\n\t"
	"movapd %3,%%xmm1\n\t"
	"addpd %5,%%xmm1\n\t"
	"movapd %%xmm1,%1\n\t"
	:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real), "m" (lhs.imag), "m" (rhs.real), "m" (rhs.imag)
	:"xmm0", "xmm1");
  return(temp);
}

inline int_complex operator + (const int_complex &lhs, const int_double &rhs)
{
	int_complex temp;
	__asm("movapd %2,%%xmm0\n\t"
		"addpd %4,%%xmm0\n\t"
		"movapd %%xmm0,%0\n\t"
		"movapd %3,%%xmm1\n\t"
		"movapd %%xmm1,%1\n\t"
		:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real), "m" (lhs.imag), "m" (rhs)
	:"xmm0", "xmm1");
	return(temp);
}

inline int_complex operator + (const int_double &rhs, const int_complex &lhs)
{
	int_complex temp;
	__asm("movapd %2,%%xmm0\n\t"
		"addpd %4,%%xmm0\n\t"
		"movapd %%xmm0,%0\n\t"
		"movapd %3,%%xmm1\n\t"
		"movapd %%xmm1,%1\n\t"
		:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real), "m" (lhs.imag), "m" (rhs)
	:"xmm0", "xmm1");
	return(temp);
}

inline int_complex operator + (const int_complex &lhs, const double &rhs)
{
	int_complex temp;
	__asm("movddup %4,%%xmm0\n\t"				// [rhs,rhs]
		"xorpd %5,%%xmm0\n\t"				// [rhs,-rhs]
		"addpd %2,%%xmm0\n\t"			// real+[rhs,-rhs]
		"movapd %%xmm0,%0\n\t"
		"movapd %3,%%xmm1\n\t"
		"movapd %%xmm1,%1\n\t"
		:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real),"m" (lhs.imag), "m" (rhs), "m" (d_zero)
		:"xmm0", "xmm1");
	return(temp);
}

inline int_complex operator + (const double &rhs, const int_complex &lhs)
{
	int_complex temp;
	__asm("movddup %4,%%xmm0\n\t"				// [rhs,rhs]
		"xorpd %5,%%xmm0\n\t"				// [rhs,-rhs]
		"addpd %2,%%xmm0\n\t"			// real+[rhs,-rhs]
		"movapd %%xmm0,%0\n\t"
		"movapd %3,%%xmm1\n\t"
		"movapd %%xmm1,%1\n\t"
		:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real),"m" (lhs.imag), "m" (rhs), "m" (d_zero)
		:"xmm0", "xmm1");
	return(temp);
}

inline int_complex operator += (int_complex &lhs, const int_complex &rhs)
{
	lhs.real+=rhs.real;
	lhs.imag+=rhs.imag;
	return(lhs);
}

inline int_complex operator += (int_complex &lhs, const int_double &rhs)
{
	lhs.real+=rhs;
	return(lhs);
}

inline int_complex operator += (int_complex &lhs, const double &rhs)
{
	lhs.real+=rhs;
	return(lhs);
}

inline int_complex operator - (const int_complex &lhs)
{
  int_complex temp;
  __asm("movapd %2,%%xmm0\n\t"		// lhs.real
	"shufpd $1,%%xmm0,%%xmm0\n\t"	// -lhs.real
	"movapd %%xmm0,%0\n\t"		// to temp.real
	"movapd %3,%%xmm1\n\t"		// lhs.imag
	"shufpd $1,%%xmm1,%%xmm1\n\t"	// -lhs.imag
	"movapd %%xmm1,%1\n\t"		// to temp.imag
	:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real), "m" (lhs.imag)
	:"xmm0", "xmm1");
  return(temp);
}

inline int_complex operator - (const int_complex &lhs, const int_complex &rhs)
{
  int_complex temp;
  __asm("movapd %4,%%xmm0\n\t"	// rhs.real
	"shufpd $1,%%xmm0,%%xmm0\n\t"	// -rhs.real
	"addpd %2,%%xmm0\n\t"		// lhs.real-rhs.real
	"movapd %%xmm0,%0\n\t"
	"movapd %5,%%xmm1\n\t"	// rhs.imag
	"shufpd $1,%%xmm1,%%xmm1\n\t"	// -rhs.imag
	"addpd %3,%%xmm1\n\t"		// lhs.imag-rhs.imag
	"movapd %%xmm1,%1\n\t"
	:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real), "m" (lhs.imag), "m" (rhs.real), "m" (rhs.imag)
	:"xmm0", "xmm1");
  return(temp);
}

inline int_complex operator - (const int_complex &lhs, const int_double &rhs)
{
	int_complex temp;
	__asm("movapd %4,%%xmm0\n\t"		// rhs
		"shufpd $1,%%xmm0,%%xmm0\n\t"	// -rhs
		"addpd %2,%%xmm0\n\t"		// lhs.real-rhs
		"movapd %%xmm0,%0\n\t"		// to temp.real
		"movapd %3,%%xmm1\n\t"
		"movapd %%xmm1,%1\n\t"
		:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real), "m" (lhs.imag), "m" (rhs)
	:"xmm0", "xmm1");
	return(temp);
}

inline int_complex operator - (const int_complex &lhs, const double &rhs)
{
	int_complex temp;
	__asm("movddup %4,%%xmm0\n\t"
		"movapd %2,%%xmm1\n\t"
		"addsubpd %%xmm0,%%xmm1\n\t"
		"movapd %%xmm1,%0\n\t"
		"movapd %3,%%xmm1\n\t"
		"movapd %%xmm1,%1\n\t"
		:"=m" (temp.real), "=m" (temp.imag)
		:"m" (lhs.real), "m" (lhs.imag), "m" (rhs)
		:"xmm0", "xmm1");
	return(temp);
}

inline int_complex operator -= (int_complex &lhs, const int_complex &rhs)
{
	lhs.real-=rhs.real;
	lhs.imag-=rhs.imag;
	return(lhs);
}

inline int_complex operator -= (int_complex &lhs, const int_double &rhs)
{
	lhs.real-=rhs;
	return(lhs);
}

inline int_complex operator -= (int_complex &lhs, const double &rhs)
{
	lhs.real-=rhs;
	return(lhs);
}

inline int_complex operator * (const int_complex &lhs, const int_complex &rhs)
{
  int_complex temp;
  __asm("movapd %3,%%xmm0\n\t"
	"xorpd	%6,%%xmm0\n\t"			// a=xmm0
	"movapd %5,%%xmm1\n\t"		// b=xmm1
	"movapd %%xmm1,%%xmm2\n\t"
	"shufpd	$1,%%xmm2,%%xmm2\n\t"		// p=xmm2
	"movapd %%xmm0,%%xmm3\n\t"
	"cmpltpd %7,%%xmm3\n\t"			// c=xmm3
	"xorpd	%7,%%xmm2\n\t"			// p=? d=xmm2
	"movapd %%xmm3,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=xmm4
	"movapd %%xmm3,%%xmm5\n\t"
	"andpd %%xmm2,%%xmm5\n\t"		// e=xmm5
	"andnpd %%xmm1,%%xmm3\n\t"		// c=? f=xmm3
	"movapd %%xmm4,%%xmm6\n\t"
	"andnpd %%xmm1,%%xmm6\n\t"		// k=xmm6
	"andpd %%xmm4,%%xmm2\n\t"		// d=? i=xmm2
	"orpd %%xmm5,%%xmm3\n\t"		// f=? h=xmm3
	"movapd %%xmm0,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=? l=xmm4
	"orpd %%xmm6,%%xmm2\n\t"		// i=? m=xmm2
	"mulpd %%xmm3,%%xmm0\n\t"		// a=? j=xmm0
	"mulpd %%xmm4,%%xmm2\n\t"		// m=? n=xmm2
	"minpd %%xmm2,%%xmm0\n\t"		// j=? o=xmm0
	"movapd %%xmm0,%%xmm7\n\t"		// im*im
	"shufpd $1,%%xmm7,%%xmm7\n\t"		// -(im*im)
	
	"movapd %2,%%xmm0\n\t"
	"xorpd	%6,%%xmm0\n\t"			// a=xmm0
	"movapd %4,%%xmm1\n\t"		// b=xmm1
	"movapd %%xmm1,%%xmm2\n\t"
	"shufpd	$1,%%xmm2,%%xmm2\n\t"		// p=xmm2
	"movapd %%xmm0,%%xmm3\n\t"
	"cmpltpd %7,%%xmm3\n\t"			// c=xmm3
	"xorpd	%7,%%xmm2\n\t"			// p=? d=xmm2
	"movapd %%xmm3,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=xmm4
	"movapd %%xmm3,%%xmm5\n\t"
	"andpd %%xmm2,%%xmm5\n\t"		// e=xmm5
	"andnpd %%xmm1,%%xmm3\n\t"		// c=? f=xmm3
	"movapd %%xmm4,%%xmm6\n\t"
	"andnpd %%xmm1,%%xmm6\n\t"		// k=xmm6
	"andpd %%xmm4,%%xmm2\n\t"		// d=? i=xmm2
	"orpd %%xmm5,%%xmm3\n\t"		// f=? h=xmm3
	"movapd %%xmm0,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=? l=xmm4
	"orpd %%xmm6,%%xmm2\n\t"		// i=? m=xmm2
	"mulpd %%xmm3,%%xmm0\n\t"		// a=? j=xmm0
	"mulpd %%xmm4,%%xmm2\n\t"		// m=? n=xmm2
	"minpd %%xmm2,%%xmm0\n\t"		// j=? o=xmm0 = re*re
	"addpd %%xmm7,%%xmm0\n\t"		// re*re-im*im
	"movapd %%xmm0,%0\n\t"

	"movapd %3,%%xmm0\n\t"			// lhs.imag
	"xorpd	%6,%%xmm0\n\t"			// a=xmm0
	"movapd %4,%%xmm1\n\t"		// rhs.real=b=xmm1 
	"movapd %%xmm1,%%xmm2\n\t"
	"shufpd	$1,%%xmm2,%%xmm2\n\t"		// p=xmm2
	"movapd %%xmm0,%%xmm3\n\t"
	"cmpltpd %7,%%xmm3\n\t"			// c=xmm3
	"xorpd	%7,%%xmm2\n\t"			// p=? d=xmm2
	"movapd %%xmm3,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=xmm4
	"movapd %%xmm3,%%xmm5\n\t"
	"andpd %%xmm2,%%xmm5\n\t"		// e=xmm5
	"andnpd %%xmm1,%%xmm3\n\t"		// c=? f=xmm3
	"movapd %%xmm4,%%xmm6\n\t"
	"andnpd %%xmm1,%%xmm6\n\t"		// k=xmm6
	"andpd %%xmm4,%%xmm2\n\t"		// d=? i=xmm2
	"orpd %%xmm5,%%xmm3\n\t"		// f=? h=xmm3
	"movapd %%xmm0,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=? l=xmm4
	"orpd %%xmm6,%%xmm2\n\t"		// i=? m=xmm2
	"mulpd %%xmm3,%%xmm0\n\t"		// a=? j=xmm0
	"mulpd %%xmm4,%%xmm2\n\t"		// m=? n=xmm2
	"minpd %%xmm2,%%xmm0\n\t"		// j=? o=xmm0
	"movapd %%xmm0,%%xmm7\n\t"		// im*re

	"movapd %2,%%xmm0\n\t"			// lhs.real
	"xorpd	%6,%%xmm0\n\t"			// a=xmm0
	"movapd %5,%%xmm1\n\t"		// rhs.imag=b=xmm1 
	"movapd %%xmm1,%%xmm2\n\t"
	"shufpd	$1,%%xmm2,%%xmm2\n\t"		// p=xmm2
	"movapd %%xmm0,%%xmm3\n\t"
	"cmpltpd %7,%%xmm3\n\t"			// c=xmm3
	"xorpd	%7,%%xmm2\n\t"			// p=? d=xmm2
	"movapd %%xmm3,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=xmm4
	"movapd %%xmm3,%%xmm5\n\t"
	"andpd %%xmm2,%%xmm5\n\t"		// e=xmm5
	"andnpd %%xmm1,%%xmm3\n\t"		// c=? f=xmm3
	"movapd %%xmm4,%%xmm6\n\t"
	"andnpd %%xmm1,%%xmm6\n\t"		// k=xmm6
	"andpd %%xmm4,%%xmm2\n\t"		// d=? i=xmm2
	"orpd %%xmm5,%%xmm3\n\t"		// f=? h=xmm3
	"movapd %%xmm0,%%xmm4\n\t"
	"shufpd $1,%%xmm4,%%xmm4\n\t"		// g=? l=xmm4
	"orpd %%xmm6,%%xmm2\n\t"		// i=? m=xmm2
	"mulpd %%xmm3,%%xmm0\n\t"		// a=? j=xmm0
	"mulpd %%xmm4,%%xmm2\n\t"		// m=? n=xmm2
	"minpd %%xmm2,%%xmm0\n\t"		// j=? o=xmm0
	"addpd %%xmm0,%%xmm7\n\t"		// im*re+re*im
	"movapd %%xmm7,%1"

	:"=m" (temp.real), "=m" (temp.imag)
	:"m" (lhs.real), "m" (lhs.imag), "m" (rhs.real), "m" (rhs.imag), "m" (d_zero), "m" (d_neg_zero)
	:"xmm0","xmm1","xmm2","xmm3","xmm4","xmm5","xmm6","xmm7");
	return(temp);
}

inline int_complex operator * (const int_complex &lhs, const int_double &rhs)
{
	return(int_complex(lhs.real*rhs,lhs.imag*rhs));
}

inline int_complex operator * (const int_double &rhs, const int_complex &lhs)
{
	return(int_complex(lhs.real*rhs,lhs.imag*rhs));
}

inline int_complex operator * (const int_complex &lhs, const double &rhs)
{
	return(int_complex(lhs.real*rhs,lhs.imag*rhs));
}

inline int_complex operator * (const double &rhs, const int_complex &lhs)
{
	return(int_complex(lhs.real*rhs,lhs.imag*rhs));
}
inline int_complex operator *= (int_complex &lhs,const int_complex &rhs)
{
	return(lhs=lhs*rhs);
}

inline int_complex operator *= (int_complex &lhs, const int_double &rhs)
{
	return(lhs=lhs*rhs);
}

inline int_complex operator *= (int_complex &lhs, const double &rhs)
{
	return(lhs=lhs*rhs);
}

inline int_complex operator / (const int_complex &lhs, const int_complex &rhs)
{
	return(lhs*conj(rhs)/norm(rhs));
}

inline int_complex operator / (const int_complex &lhs, const int_double &rhs)
{
	return(int_complex(lhs.real/rhs,lhs.imag/rhs));
}

inline int_complex operator / (const int_complex &lhs, const double &rhs)
{
	return(int_complex(lhs.real/rhs,lhs.imag/rhs));
}

inline int_complex operator /= (int_complex &lhs, const int_complex &rhs)
{
	return(lhs=lhs/rhs);
}

inline int_complex operator /= (int_complex &lhs, const int_double &rhs)
{
	return(lhs=lhs/rhs);
}

inline int_complex operator /= (int_complex &lhs, const double &rhs)
{
	return(lhs=lhs/rhs);
}

inline int_complex exp(const int_complex &z)
{
	int_double xs,xc;
	sin_cos(z.imag,&xs,&xc);
	return(int_complex(xc,xs)*exp(z.real));
}

inline int_complex e(const int_double &x)
{
	int_double xs,xc;
	sin_cospi(x*2,&xs,&xc);
	return(int_complex(xc,xs));
}

inline int_complex sqrt(const int_complex &z) // returns sqrt with arg in [-pi/2,pi/2]
{
  int_double theta,mod_z;
  int_complex res;
  if(contains_zero(z.imag)&&(z.real.right>0.0)) // don't try to atan this
    return(sqrt(-z)*int_complex(0,1));
  mod_z=pow(norm(z),0.25);
  theta=atan2(z.imag,z.real)/2.0;
  sin_cos(theta,&res.imag,&res.real);
  res=res*mod_z;
  return(res);
}

inline void print_int_complex(const int_complex &z)
{
  print_int_double(z.real);
  printf("+");
  print_int_double(z.imag);
  printf("i");
}

inline int_complex pow (const int_double &x, const int_complex &s)
{
  int_double lnx;
  int_complex z;
  
  lnx=log(x);
  sin_cos(lnx*s.imag,&z.imag,&z.real);
  return(z*exp(lnx*s.real));
}

inline int_complex pow (const double x, const int_complex &s)
{
  return(pow(int_double(x),s));
}


inline int_complex pow1 (const double &x, const double &re, const double &im)
{
  return(pow1(int_double(x),re,im));
}

inline int_double norm(const int_complex &z)
{
  return(sqr(z.real)+sqr(z.imag));
}

inline int_double mod(const int_complex &z)
{
  return(sqrt(norm(z)));
}

inline int_double argument(const int_complex &z)
{
  return(atan2(z.imag,z.real));
}

inline int_complex log(const int_complex &z)
{
  return(int_complex(log(norm(z))/2,argument(z)));
}

inline int_double read_int_double(FILE *infile)
{
  int y;
  double x[2];
  y=fread(x,sizeof(double),2,infile);
  return(int_double(x[0],x[1]));
}

inline int_complex read_int_complex(FILE *infile)
{
  int y;
  double x[4];
  y=fread(x,sizeof(double),4,infile);
  return(int_complex(int_double(x[0],x[1]),int_double(x[2],x[3])));
}

#define _aligned_malloc(a,b) memalign(b,a)
#define _aligned_free(x) free(x)   // I believe GCC free can handle memory allocated by memalign
unsigned short int old_cw,new_cw;
unsigned int old__SSE_cw,new__SSE_cw;
#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  :: "m" (cw))
#define __SSE_getcw(cw) asm ("stmxcsr %0" : "=m" (cw))
#define __SSE_setcw(cw) asm ("ldmxcsr %0" :: "m" (cw))

void _fpu_rndd() 
{
  // must leave fpu rounding to crlibm
  //_fpu_getcw(old_cw);
  //new_cw=old_cw&0x3FF;
  //new_cw=new_cw|0x400;
  //_fpu_setcw(new_cw);

  crlibm_init();
  __SSE_getcw(old__SSE_cw);
  //printf("Old SSE Control Register was %X\n",old__SSE_cw);
  new__SSE_cw=old__SSE_cw&0x1FBF; // zeros FTZ(15),RC(14,13) and DAZ(6) bits
  new__SSE_cw=new__SSE_cw|0x2000; // sets RC to Round Down (14=0, 13=1)
  //printf("New SSE Control Register is %X\n",new__SSE_cw);
  __SSE_setcw(new__SSE_cw);
  //printf("Set FPU control word to %X and SSE control word to %X\n",new_cw,new__SSE_cw);
  d_zero_zero.left=0.0;
  d_zero_zero.right=0.0;
  d_zero.left=0.0;
  d_zero.right=_nze;
  d_neg_zero.left=_nze;
  d_neg_zero.right=_nze;
  d_neg_neg_zero.left=_nze;
  d_neg_neg_zero.right=0.0;
  d_ln_pi=log(d_pi);
  d_ln_two_pi=log(d_two_pi);

}

#endif
