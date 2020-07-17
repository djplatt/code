//
// File int_double
//
// Define interval arithmetic on doubles
// extend to complexes using rectangles
//
//
// Use SSE op-codes for IEEE compliance and
// efficiency. See Lambov.
//
// Could write operators as assembly routines
// in a library but wanted the efficiency of 
//  inlining.
//
// (Abridged) change log
//
// 22/3/9 v 8.0 added use of nextbefore
// 15/4/9 v 9.0 added SSE assembler for GCC
// 30/4/9 v 10.0 added SSE for GCC + Intel Linux
// 24/7/9 v 11.0 coded log and exp in FP assembler
// 27/6/10 v 12.0 replaced FPU calls to sin, cos, atan, exp, log
//                with calls to crlibm functions. Can't trust FPU!
//                No longer mess with FPU rounding so crlibm has a
//                chance.
// 14/6/13 v 13.0 Using intrinsics
//
// To do: Use GCC intrinsics to let the compiler try to optimise
//        SSE register usage.
//
// g++ compile flags:- -fomit-frame-pointer -O1 -march=nocona -msse3 -mfpmath=387 -frounding-math -finline-functions

#ifndef INT_DOUBLE13
#define INT_DOUBLE13
#define LINUX
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include "malloc.h"
#include "crlibm.h"
#include <emmintrin.h>

using namespace std;

#define debug printf("Reached line number %d\n",__LINE__)




#define _aligned_malloc(a,b) memalign(b,a)
#define _aligned_free(x) free(x)   // I believe GCC free can handle memory allocated by memalign
unsigned short int old_cw,new_cw;
unsigned int old__SSE_cw,new__SSE_cw;
#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  :: "m" (cw))
#define __SSE_getcw(cw) asm ("stmxcsr %0" : "=m" (cw))
#define __SSE_setcw(cw) asm ("ldmxcsr %0" :: "m" (cw))

__v2df d_zero_zero;                
__v2df d_zero;  // +0.0 -0.0       
__v2df d_neg_zero; // -0.0 -0.0    
__v2df d_neg_neg_zero; // -0.0 +0.0

double nextafter(const double);
double nextbefore(const double);


class int_double
{
 public:
  __v2df v;
  

  inline void set_left(const double x)
  {
    this->v=_mm_loadl_pd(this->v,&x);
    //((double *) &v)[0]=x;
  }

  inline void set_right(const double x)
  {
    this->v=_mm_loadh_pd(this->v,&x);
    //((double *) &v)[1]=x;
  }

  int_double()
    {
    }

  inline int_double (double x)
    {
      __attribute__(aligned(16)) double a[2];
      a[0]=x;a[1]=x;
      __m128d b=_mm_load_pd(a);
      this->v=_mm_xor_pd(b,d_zero);
      //((double *) &v)[0] = x;
      //((double *) &v)[1] = -x;
    }

  inline int_double(double x, double y)
    {
      __attribute__(aligned(16)) double a[2];
      a[0]=x;a[1]=y;
      __m128d c=_mm_load_pd(a);
      this->v=_mm_xor_pd(c,d_zero);
      //((double *) &v)[0] = x;
      //((double *) &v)[1] = -y;
    }

  inline int_double(__v2df v)
    {
      this->v = v;
    }

  inline friend int_double operator+(const int_double lhs, const int_double rhs)
  {
    return _mm_add_pd(lhs.v, rhs.v);
  }

  inline friend int_double operator+(const int_double lhs, const double rhs)
  {
    return _mm_add_pd(_mm_xor_pd(d_zero,_mm_load1_pd(&rhs)), lhs.v);
  }

  inline friend int_double operator+(const double lhs, const int_double rhs)
  {
    return _mm_add_pd(_mm_xor_pd(d_zero,_mm_load1_pd(&lhs)), rhs.v);
  }

  inline friend int_double operator +=(int_double &lhs, const int_double rhs)
  {
    return(lhs=lhs+rhs);
  }

  inline friend int_double operator +=(int_double &lhs, const double rhs)
  {
    return(lhs=lhs+rhs);
  }

  inline friend int_double operator-(const int_double lhs, const int_double rhs)
  {
    return _mm_add_pd(lhs.v, _mm_shuffle_pd(rhs.v, rhs.v, 1));
  }

  inline friend int_double operator-(const int_double lhs)
  {
    return _mm_shuffle_pd(lhs.v, lhs.v, 1);
  }

  inline friend int_double operator -=(int_double &lhs, const int_double rhs)
  {
    return(lhs=lhs-rhs);
  }

  inline friend int_double operator -=(int_double &lhs, const double rhs)
  {
    return(lhs=lhs-rhs);
  }


  inline friend int_double operator*(const int_double lhs, const int_double rhs)
  {
    __m128d a = _mm_xor_pd(lhs.v, d_zero);
    //__m128d b = rhs.v; // no op
    __m128d p = _mm_shuffle_pd(rhs.v, rhs.v, 1);
    __m128d c = _mm_cmplt_pd(a, d_neg_zero);
    __m128d d = _mm_xor_pd(p, d_neg_zero);
    __m128d g = _mm_shuffle_pd(c, c, 1);
    __m128d e = _mm_and_pd(c, d);
    __m128d f = _mm_andnot_pd(c, rhs.v);
    __m128d k = _mm_andnot_pd(g, rhs.v);
    __m128d i = _mm_and_pd(g, d);
    __m128d h = _mm_or_pd(e, f);
    __m128d l = _mm_shuffle_pd(a, a, 1);
    __m128d m = _mm_or_pd(i, k);
    __m128d j = _mm_mul_pd(a, h);
    __m128d n = _mm_mul_pd(l, m);
    return _mm_min_pd(j, n);
  }

  inline friend int_double operator *=(int_double &lhs, const int_double rhs)
  {
    return(lhs=lhs*rhs);
  }

  inline friend int_double operator * (const int_double lhs, const double rhs)
  {
    //int_double temp;
    __m128d a = _mm_load1_pd(&rhs);
    __m128d b = _mm_xor_pd(a,d_neg_zero);
    __m128d c = _mm_mul_pd(lhs.v,a);
    __m128d e = _mm_mul_pd(lhs.v,b);
    __m128d f = _mm_shuffle_pd(e,e,1);
    __m128d g = _mm_min_pd(f,c);
    return(g);
  }

  inline friend int_double operator *=(int_double &lhs, const double rhs)
  {
    return(lhs=lhs*rhs);
  }

  inline friend int_double operator /(const int_double lhs, const int_double rhs)
  {
    __m128d x0=_mm_xor_pd(d_zero,rhs.v);
    __m128d x1=_mm_shuffle_pd(x0,x0,1);
    __m128d x2=_mm_cmpnle_pd(x0,d_zero_zero);
    __m128d x3=lhs.v;
    __m128d x4=x3;
    x3=_mm_shuffle_pd(x3,x3,1);
    x3=_mm_xor_pd(d_neg_zero,x3);
    x4=_mm_and_pd(x2,x4);
    x2=_mm_andnot_pd(x2,x3);
    x2=_mm_or_pd(x4,x2);
    x3=x2;
    x2=_mm_div_pd(x2,x1);
    x3=_mm_div_pd(x3,x0);
    return(_mm_min_pd(x2,x3));
  }

  inline friend int_double operator /=(int_double &lhs, const int_double rhs)
  {
    return(lhs=lhs/rhs);
  }

  inline friend int_double operator /(const int_double lhs, const double rhs)
  {
    __m128d x0=_mm_load1_pd(&rhs);
    //__m128d x0=_mm_xor_pd(d_zero,x0);
    __m128d x1=_mm_shuffle_pd(x0,x0,1);
    __m128d x2=_mm_cmpnle_pd(x0,d_zero_zero);
    __m128d x3=lhs.v;
    __m128d x4=x3;
    x3=_mm_shuffle_pd(x3,x3,1);
    x3=_mm_xor_pd(d_neg_zero,x3);
    x4=_mm_and_pd(x2,x4);
    x2=_mm_andnot_pd(x2,x3);
    x2=_mm_or_pd(x4,x2);
    x3=x2;
    x2=_mm_div_pd(x2,x1);
    x3=_mm_div_pd(x3,x0);
    return(_mm_min_pd(x2,x3));
  }

  inline friend int_double operator /=(int_double &lhs, const double rhs)
  {
    return(lhs=lhs/rhs);
  }

};

inline double left(const int_double x)
{
  return(((double *) &x.v)[0]);
}

inline double right(const int_double x)
{
  return(((double *) &x.v)[1]);
}

// we can do better
// 
inline int_double sqrt(const int_double x) // strictly increasing
{
  __m128d a=_mm_xor_pd(d_zero,x.v);
  int_double res=_mm_sqrt_pd(a);
  res.set_right(-nextafter(right(res)));
  return(res);
}

inline int_double sqr(const int_double x)
{
  __m128d x0=x.v;
  __m128d x1=x0;
  __m128d x2=x0;
  x1=_mm_shuffle_pd(x1,x1,1);
  x0=_mm_min_sd(x0,x1);
  x1=_mm_max_sd(x1,x2);
  x1=_mm_max_sd(x1,d_zero_zero);
  x1=_mm_unpacklo_pd(x1,x0);
  x0=x1;
  x1=_mm_xor_pd(d_zero,x1);
  return(_mm_mul_pd(x1,x0));
}

inline int contains_zero(const int_double x)
{
  return((left(x)<=0)&&(right(x)<=0));
}

// Useful constants
//
// a half
int_double d_half;
//
// Pi
int_double d_pi;
//
// Euler gamma
int_double d_gamma;
//
// Pi/2
int_double d_pi_2;
//
// 2*Pi
int_double d_two_pi;
//
// log(pi)
int_double d_ln_pi;
//
// log(2*pi)
int_double d_ln_two_pi;


void print_int_double(const int_double x)
{
  printf("[%20.18e,%20.18e]",left(x),-right(x));
}

void print_int_double_str(const char *str, const int_double x)
{
  printf("%s",str);
  print_int_double(x);
  printf("\n");
}

inline double width(const int_double x)
{
  return(-right(x)-left(x));
}

inline int operator >= (const int_double lhs, const int_double rhs)
{
  return(left(lhs)>=-right(rhs));
}

inline int operator > (const int_double lhs, const int_double rhs)
{
  return(left(lhs)>-right(rhs));
}

inline int operator < (const int_double lhs, const int_double rhs)
{
  return(-right(lhs)<left(rhs));
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
  return(int_double(exp_rd(left(x)),exp_ru(-right(x))));
}

// if right hand endpoint is one, then interval is
// widened unnecessarily. could trap it but......
inline int_double log (const int_double &x) // nicely increasing
{
  return(int_double(log_rd(left(x)),log_ru(-right(x))));
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
  if(b>0.0)
    mask|=0b0100;
  if(c<0.0)
    mask|=0b0010;
  if(d>0.0)
    mask|=0b0001;
  return(mask);
}

inline int make_mask(const int_double &x, const int_double &y)
{
  return(make_mask1(left(x), right(x), left(y), right(y)));
}

// move a double prec float towards +infinity by
// the smallest delta possible
int_double delta_int_double;
int_double delta_int_double_neg;

inline double nextafter(const double x)
{
  int_double xd(x);
  xd+=delta_int_double;
  return(-right(xd));
}

inline double nnextafter(const double x)
{
  int_double xd(x);
  xd+=delta_int_double;
  return(right(xd));
}



inline double nextbefore(double x)
{
  int_double xd(x);
  xd-=delta_int_double;
  return(left(xd));
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
  //printf("%30.28e\n%30.28e\n%30.28e\n%30.28e\nmask=%x\n",cos_left,cos_right,sin_left,sin_right,mask);
  switch (mask)
    {
    case 0: // all -ve, 3rd quadrant
      sin_x->set_left(sin_right);
      sin_x->set_right(nnextafter(sin_left));
      cos_x->set_left(cos_left);
      cos_x->set_right(nnextafter(cos_right));
      return;
    case 2:
      sin_x->set_left(nextbefore(sin_right));  // sin_right < 0 and too close to 0
      sin_x->set_right(nnextafter(sin_left));
      cos_x->set_left(-1.0);
      if(cos_left>=cos_right)
	cos_x->set_right(nnextafter(cos_left));
      else
	cos_x->set_right(nnextafter(cos_right));
      return;
    case 3:
      sin_x->set_left(sin_right);
      sin_x->set_right(nnextafter(sin_left));
      cos_x->set_left(cos_right);
      cos_x->set_right(nnextafter(cos_left));
      return;
    case 4:
      sin_x->set_left(-1.0);
      if(sin_right>=sin_left)
	sin_x->set_right(nnextafter(sin_right));
      else
	sin_x->set_right(nnextafter(sin_left));
      cos_x->set_left(nextbefore(cos_left));
      cos_x->set_right(nnextafter(cos_right));
      return;
    case 11:
      if(sin_left<=sin_right)
	sin_x->set_left(sin_left);
      else
	sin_x->set_left(sin_right);
      sin_x->set_right(-1.0);
      cos_x->set_left(cos_right);
      cos_x->set_right(nnextafter(cos_left));
      return;
    case 12:
      sin_x->set_left(sin_left);
      sin_x->set_right(nnextafter(sin_right));
      cos_x->set_left(cos_left);
      cos_x->set_right(nnextafter(cos_right));
      return;
    case 13:
      cos_x->set_right(-1.0);
      if(cos_left<=cos_right)
	cos_x->set_left(cos_left);
      else
	cos_x->set_left(cos_right);
      sin_x->set_left(sin_left);
      sin_x->set_right(nnextafter(sin_right));
      return;
    case 15:
      sin_x->set_left(sin_left);
      sin_x->set_right(nnextafter(sin_right));
      cos_x->set_left(cos_right);
      cos_x->set_right(nnextafter(cos_left));
      return;
    default:
      printf("Weird error in sin_cos, mask = %d Exiting.\n",(int) mask);
      printf("\ncos_left: %20.18e\ncos_right: %20.18e\n",cos_left,cos_right);
      printf("sin_left: %20.18e\nsin_right: %20.18e\n",sin_left,sin_right);
      exit(1);
    }
}

// perform sin(x) into sin_x and cos(x) into cos_x
inline void sin_cos(const int_double &x, int_double *sin_x, int_double *cos_x)
{
  double cos_left,cos_right,sin_left,sin_right;
  unsigned char mask=0;

  int_double xp=x/d_pi;
  if(right(xp)+left(xp)<=-0.5)
    {
      printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");
      exit(1);
      return;
    }

  cos_left=cospi_rd(left(xp));
  cos_right=cospi_rd(-right(xp)); // - probably a nop because cos is even
  sin_left=sinpi_rd(left(xp));
  sin_right=sinpi_rd(-right(xp));
  sin_cos1(cos_left,cos_right,sin_left,sin_right,sin_x,cos_x);
}

//
// sin (pi*x), cos(pi*x)
//
inline void sin_cospi(const int_double &x, int_double *sin_x, int_double *cos_x)
{
  double cos_left,cos_right,sin_left,sin_right;


  if(right(x)+left(x)<=-0.5)
    {
      printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");
      exit(1);
      return;
    }

  cos_left=cospi_rd(left(x));
  cos_right=cospi_rd(-right(x)); // probably a nop because cos is even
  sin_left=sinpi_rd(left(x));
  sin_right=sinpi_rd(-right(x));
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

// simple atan
inline int_double atan(const int_double &x)
{
  return(int_double(atan_rd(left(x)),atan_ru(-right(x))));
}

// returns theta in [-Pi,Pi]
inline int_double atan2(const int_double &y, const int_double &x)
{
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

inline int_double pow(const int_double x, const int_double y)
{
  return(exp(y*log(x)));
}

inline int_double pow(const int_double x, const double y)
{
  return(exp(log(x)*y));
}


class int_complex{
 public:
  int_double real;
  int_double imag;

  inline int_complex()
    {
      real=int_double(0.0);
      imag=real;
    };

  inline int_complex(const int_double re)
    {
      real=re;
      imag=int_double(0.0);
    };

  inline int_complex(const int_double re, const int_double im)
    {
      real=re;
      imag=im;
    };

  inline friend int_complex operator +(const int_complex lhs, const int_complex rhs)
  {
    return(int_complex(lhs.real+rhs.real,lhs.imag+rhs.imag));
  }

  inline friend int_complex operator +(const int_complex lhs, const int_double rhs)
  {
    return(int_complex(lhs.real+rhs,lhs.imag));
  }

  inline friend int_complex operator +(const int_complex lhs, const double rhs)
  {
    return(int_complex(lhs.real+rhs,lhs.imag));
  }

  inline friend int_complex operator +=(int_complex &lhs, const int_complex rhs)
  {
    return(lhs=lhs+rhs);
  }

  inline friend int_complex operator -(const int_complex lhs, const int_complex rhs)
  {
    return(int_complex(lhs.real-rhs.real,lhs.imag-rhs.imag));
  }

  inline friend int_complex operator -(const int_complex lhs, const int_double rhs)
  {
    return(int_complex(lhs.real-rhs,lhs.imag));
  }

  inline friend int_complex operator -(const int_complex lhs, const double rhs)
  {
    return(int_complex(lhs.real-rhs,lhs.imag));
  }

  inline friend int_complex operator -=(int_complex &lhs, const int_complex rhs)
  {
    return(lhs=lhs-rhs);
  }

  inline friend int_complex operator -(const int_complex lhs)
  {
    return(int_complex(-lhs.real,-lhs.imag));
  }

  inline friend int_complex operator *(const int_complex lhs, const int_complex rhs)
  {
    return(int_complex(lhs.real*rhs.real-lhs.imag*rhs.imag,
		       lhs.real*rhs.imag+lhs.imag*rhs.real));
  }

  inline friend int_complex operator *(const int_complex lhs, const int_double rhs)
  {
    return(int_complex(lhs.real*rhs,lhs.imag*rhs));
  }

  inline friend int_complex operator *(const int_complex lhs, const double rhs)
  {
    return(int_complex(lhs.real*rhs,lhs.imag*rhs));
  }

  inline friend int_complex operator *=(int_complex &lhs, const int_complex rhs)
  {
    return(lhs=lhs*rhs);
  }


  inline friend int_complex operator /(const int_complex, const int_complex);

  inline friend int_complex operator /(const int_complex lhs, const int_double rhs)
  {
    return(int_complex(lhs.real/rhs,lhs.imag/rhs));
  }

  inline friend int_complex operator /(const int_complex lhs, const double rhs)
  {
    return(int_complex(lhs.real/rhs,lhs.imag/rhs));
  }

  inline friend int_complex operator /=(int_complex &lhs, const int_complex rhs)
  {
    return(lhs=lhs/rhs);
  }

  inline friend int_complex operator /=(int_complex &lhs, const int_double rhs)
  {
    return(lhs=lhs/rhs);
  }

  inline friend int_complex operator /=(int_complex &lhs, const double rhs)
  {
    return(lhs=lhs/rhs);
  }

};

void print_int_complex(const int_complex z)
{
  print_int_double(z.real);
  printf("+i");
  print_int_double(z.imag);
}

void print_int_complex_str(const char *str, const int_complex z)
{
  printf("%s",str);
  print_int_complex(z);
  printf("\n");
}

inline int_complex conj(const int_complex lhs)
{
  return(int_complex(lhs.real,-lhs.imag));
}

inline int_double norm(const int_complex lhs)
{
  return(sqr(lhs.real)+sqr(lhs.imag));
}

inline int_complex operator /(const int_complex lhs, const int_complex rhs)
{
  return(lhs*conj(rhs)/norm(rhs));
}

inline int contains_zero(const int_complex z)
{
  return(contains_zero(z.real)&&contains_zero(z.imag));
}

inline int_complex pow (const int_double x, const int_complex s)
{
  int_double lnx;
  int_complex z;
  
  lnx=log(x);
  sin_cos(lnx*s.imag,&z.imag,&z.real);
  return(z*exp(lnx*s.real));
}

inline int_complex pow (const double x, const int_complex s)
{
  return(pow(int_double(x),s));
}

inline int_double mod(const int_complex z)
{
  return(sqrt(norm(z)));
}

inline int_double argument(const int_complex z)
{
  return(atan2(z.imag,z.real));
}

inline int_complex log(const int_complex z)
{
  return(int_complex(log(norm(z))/2,argument(z)));
}

inline int_complex exp(const int_complex z)
{
	int_double xs,xc;
	sin_cos(z.imag,&xs,&xc);
	return(int_complex(xc,xs)*exp(z.real));
}

inline int_complex sqrt(const int_complex z1) // returns sqrt with arg in [-pi/2,pi/2]
{
  //print_int_complex_str("In sqrt with z=",z1);
  int_double theta,mod_z;
  int_complex res,z;
  bool negated;
  if(contains_zero(z.imag)&&(right(z.real)>0.0)) // don't try to atan this
    {
      z=-z1;
      negated=true;
    }
  else
    {
      z=z1;
      negated=false;
    }

  mod_z=pow(norm(z),0.25);
  theta=atan2(z.imag,z.real)/2.0;
  sin_cos(theta,&res.imag,&res.real);
  res=res*mod_z;
  if(negated)
    return(res*int_complex(0,1));
  else
    return(res);
}

int_double d_one;

int_complex c_zero;
int_complex c_one;

void _sse_rndd() 
{
  // must leave fpu rounding alone for crlibm etc.
  //_fpu_getcw(old_cw);
  //new_cw=old_cw&0x3FF;
  //new_cw=new_cw|0x400;
  //_fpu_setcw(new_cw);

  __SSE_getcw(old__SSE_cw);
  new__SSE_cw=old__SSE_cw&0x1FBF; // zeros FTZ(15),RC(14,13) and DAZ(6) bits
  new__SSE_cw=new__SSE_cw|0x2000; // sets RC to Round Down (14=0, 13=1)
  __SSE_setcw(new__SSE_cw);
  //printf("Set SSE control word to %X\n",new__SSE_cw);

  //
  ((double *) &d_zero)[0] = 0.0;
  ((double *) &d_zero)[1] = -0.0;
  ((double *) &d_zero_zero)[0] = 0.0;
  ((double *) &d_zero_zero)[1] = 0.0;
  ((double *) &d_neg_zero)[0] = -0.0;
  ((double *) &d_neg_zero)[1] = -0.0;
  ((double *) &d_neg_neg_zero)[0] = 0.0;
  ((double *) &d_neg_neg_zero)[1] = 0.0;
// Useful constants
//
// a half
d_half=int_double(0.5,0.5);
//
// Pi
int _i_pi[2]={1413754136,1074340347};   // this is what double pi less a bit looks like
int _i_pi2[2]={1413754137,1074340347};  // this is double pi plus a bit
double *_d_pi=(double *)&_i_pi;
double *_d_pi2=(double *)&_i_pi2;
//
// Euler gamma
int _i_gamma[2]={0xfc6fb618,0x3fe2788c};
int _i_gamma2[2]={0xfc6fb619,0x3fe2788c};
double *_d_gamma=(double *)&_i_gamma;
double *_d_gamma2=(double *)&_i_gamma2;
d_gamma=int_double(_d_gamma[0],_d_gamma2[0]);
//
// Pi/2
int _i_pi_2[2]={0x54442d18,0x3ff921fb};   // this is what double pi/2 less a bit looks like
int _i_pi2_2[2]={0x54442d19,0x3ff921fb};  // this is double pi/2 plus a bit
double *_d_pi_2=(double *)&_i_pi_2;
double *_d_pi2_2=(double *)&_i_pi2_2;
d_pi_2=int_double(_d_pi_2[0],_d_pi2_2[0]);
//
// 2*Pi
int _i_2pi_2[2]={0x54442d18,0x401921fb};   // this is what double pi*2 less a bit looks like
int _i_2pi2_2[2]={0x54442d19,0x401921fb};  // this is double pi*2 plus a bit
double *_d_2pi_2=(double *)&_i_2pi_2;
double *_d_2pi2_2=(double *)&_i_2pi2_2;
d_two_pi=int_double(_d_2pi_2[0],_d_2pi2_2[0]);

 d_ln_pi=log(d_pi);
 d_ln_two_pi=log(d_two_pi);
 d_one=int_double(1.0);
 
 c_zero=int_complex(d_zero);
 c_one=int_complex(d_one);
 d_pi=int_double(_d_pi[0],_d_pi2[0]);

 delta_int_double=int_double(DBL_MIN);
 delta_int_double_neg=int_double(-DBL_MIN);

}


#endif
