// 22/3/9 v 8.0 added use of nextbefore
// 15/4/9 v 9.0 added SSE assembler for GCC
// 30/4/9 v 10.0 added SSE for GCC + Intel Linux
// gcc compile flags:- -fomit-frame-pointer -O3 -march=nocona -msse3 -frounding-math
// MS compile flags:- /Ox /Ob2 /Oi /Ot /Oy /GT /GL /D "WIN32" /D "NDEBUG" /D "_CONSOLE" 
//                    /D "_UNICODE" /D "UNICODE" /FD /EHsc /MT /Zp16 /GS- /arch:SSE2 /FAs
//                    /Fa"Release\\" /Fo"Release\\" /Fd"Release\vc90.pdb" /W3 /nologo /c
//                    /Zi /TP /errorReport:prompt
// Intel compile flags:- 
#ifndef INT_DOUBLE10
#define INT_DOUBLE10

// tell Visual Studio to define M_PI_2 etc.
// does no harm in GCC/Intel
#define _USE_MATH_DEFINES

// Running under LINUX
#define LINUX
// or
// Running under Windows
//#define WINDOWS

// Use SSE instructions
#define SSE
// or dont

// Using MS Studio C++ Compiler
//#define MS
// or
// Using GCC compiler (under LINUX/cygwin)
#define GCC
// or
// using INTEL compiler (Linux or Windows)
//#define INTEL

#if defined(SSE)
#if defined(LINUX)
// use ATT syntax for SSE assembler
#define SSE_ATT
#else
// use Intel syntax
#define SSE_INTEL
#endif
#endif

#define ATTRIBS(x) inline

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <malloc.h>

using namespace std;

#define debug printf("Reached line number %d\n",__LINE__)
#define print_int_complex_str(str,x) {printf(str);printf(" ");print_int_complex(x);printf("\n");}
#define print_int_double_str(str,x) {printf(str);printf(" ");print_int_double(x);printf("\n");}

inline double nextafter(double);
inline double nextbefore(double);

// int_double must be 16 byte aligned for SSE instructions
class int_double{
public:
#if defined(MS)
	__declspec(align(16))
#endif
	double
#if !defined(MS)
__attribute__ ((aligned(16)))
#endif
left;
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
        exit(0);
    };
    left=l;
    right=-r;
}

   friend int_double operator + (const int_double &lhs, const int_double &rhs);
   friend int_double operator + (const int_double &lhs, const double &rhs);
   friend int_double operator - (const int_double &lhs, const int_double &rhs);
   friend int_double operator - (const int_double &lhs, const double &rhs);
   friend int_double operator - (const int_double &lhs);
   friend int_double operator * (const int_double &lhs, const int_double &rhs);
   friend int_double operator * (const int_double &lhs, const double &rhs);
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
};

int_double times_pos(const int_double &, const int_double &);

int_double exp (const int_double &);

int_double log (const int_double &);

void sin_cos(const int_double &, int_double *, int_double *);

int_double sqr(const int_double &);

int_double pow(const int_double &, const int_double &);

int_double pow(const int_double &, const double &);

int_double atan(const int_double &);

int_double sqrt(const int_double &);

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
	int_double 
//#if defined(LINUX)&defined(SSE)
//__attribute__ ((aligned(16)))
//#endif
	real;
	int_double imag;
	inline int_complex ()
	{};
	inline int_complex (const int_double &re,const int_double &im)
	{
		real=re;
		imag=im;
	};

	 friend int_complex operator + (const int_complex &lhs, const int_complex &rhs);
	 friend int_complex operator * (const int_complex &lhs, const int_complex &rhs);
	 friend int_complex operator * (const int_complex &lhs, const int_double &rhs);
	 friend int_complex operator * (const int_complex &lhs, const double &rhs);
	 friend int_complex operator + (const int_complex &lhs, const int_double &rhs);
	 friend int_complex operator - (const int_complex &lhs, const int_complex &rhs);
	 friend int_complex operator - (const int_complex &lhs);
	 friend int_complex operator / (const int_complex &lhs, const int_complex &rhs);
	 friend int_complex operator + (const int_complex &lhs, const double &rhs);
	 friend int_complex operator - (const int_complex &lhs, const int_double &rhs);
	 friend int_complex operator - (const int_complex &lhs, const double &rhs);
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

void print_int_complex(const int_complex &);

int_complex pow (const int_double &,const int_complex &);

int_complex pow1(const double &, const double &, const double &);

int_complex pow1(const int_double &, const double &, const double &);

int_complex pow (const double, const int_complex &);

int_double norm(const int_complex &);

int_complex sqrt(const int_complex &);

void print_interval (const int_double &);

int_complex conj (const int_complex &);

int_double argument (const int_complex &);

int_complex lngamma (const int_complex &);

int_complex hurwitz(const int_complex &, const int_double &);

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

/*static*/ const int_double d_half=int_double(0.5,0.5);

/*static*/ int _i_pi[2]={1413754135,1074340347};   // this is what double pi less a bit looks like

/*static*/ int _i_pi2[2]={1413754137,1074340347};  // this is double pi plus a bit

/*static*/ double *_d_pi=(double *)&_i_pi;

/*static*/ double *_d_pi2=(double *)&_i_pi2;

/*static*/ int_double d_pi=int_double(_d_pi[0],_d_pi2[0]);

/*static*/ int_double d_two_pi;

/*static*/ int_double d_ln_two_pi;

/*static*/ int_double d_pi_2;

unsigned int __nze[2]={0,0x80000000}; // -0.0

unsigned int __delta[2]={1,0};       // very small

unsigned int __delta_minus[2]={1,0x80000000};

double _nze = *(double*)__nze;

double _delta= *(double*)__delta;

double _delta_minus= *(double*)__delta_minus;

int_double d_zero_zero;
/*static*/ int_double d_zero;  // +0.0 -0.0
/*static*/ int_double d_neg_zero; // -0.0 -0.0
/*static*/ int_double d_neg_neg_zero; // -0.0 +0.0
/*static*/ int_double d_one=int_double(1.0,1.0);
/*static*/ int_complex c_zero;
/*static*/ int_complex c_one;

// move a double prec float towards +infinity by
// the smallest delta possible
inline double nextafter (double x)
{
  unsigned int *i=(unsigned int *) &x;

  if(i[1]&0x80000000)   // -ve
    {
      if(x==_nze)               // -0.0
	{
	  return(_delta);
	}

      i[0]--;
      if(i[0]==0xffffffff)
	{
	  i[1]--;
	}
      return(x);
    }

  if((i[1]&0x7ff00000)==0x7ff00000) // nan or +/-inf
    {
      return(x);  // don't mess with it
    }


  i[0]++;
  if(i[0]==0)
    i[1]++;
  return(x);
}

// move a double prec float towards -infinity by
// the smallest delta possible
inline double nextbefore (double x)
{
  unsigned int *i=(unsigned int *) &x;

  if(i[1]&0x80000000)   // -ve
    {
      if(x==0.0)
	{
	  return(_delta_minus);
	}

      i[0]++;
      if(i[0]==0x00000000)
	{
	  i[1]++;
	}
      return(x);
    }

  if((i[1]&0x7ff00000)==0x7ff00000) // nan or +/-inf
    {
      return(x);  // don't mess with it
    }


  i[0]--;
  if(i[0]==0xffffffff)
    i[1]--;
  return(x);
}

void print_int_double(const int_double &x)
{
  printf("[%20.18e,%20.18e]",x.left,-x.right);
};

inline int rel_error(int_double &x)
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
			return((int)log10(rerr));
	}
}

inline int abs_error(int_double &x)
{
	double aerr;
	aerr=x.left+x.right;
	if(aerr==0.0)
		return(0);
	if(aerr<0.0)
		return((int)log10(-aerr));
	return((int)log10(aerr));
}


int_double op_temp;



ATTRIBS(3)
#if defined(SSE_ATT)
int_double operator + (const int_double &lhs, const int_double &rhs)
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
#elif defined(_SSE_INTEL)
int_double operator + (const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	__asm{
	movapd XMM0,edx
	addpd XMM0,ecx
	movapd [temp],XMM0
	};
	return(temp);
}
#else
int_double operator + (const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	temp.left=lhs.left+rhs.left;
	temp.right=lhs.right+rhs.right;
	return(temp);
}
#endif

ATTRIBS(3)
#if defined (SSE_ATT)
int_double operator + (const int_double &lhs, const double &rhs)
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
#elif defined(_SSE_INTEL)
int_double operator + (const int_double &lhs, const double &rhs)
{
	int_double temp;
	__asm{
			movddup XMM0,rhs
			xorpd XMM0,d_zero
			mov eax,dword ptr lhs
			addpd XMM0,[eax]
			movapd temp,XMM0
	}
	return(temp);
}
#else
int_double operator + (const int_double &lhs, const double &rhs)
{
	int_double temp;
	temp.left=lhs.left+rhs;
	temp.right=lhs.right-rhs;
	return(temp);
}
#endif

inline int_double operator += (int_double &lhs, const int_double &rhs)
{ return(lhs=lhs+rhs);}

inline int_double operator += (int_double &lhs, const double &rhs)
{ return(lhs=lhs+rhs);}

ATTRIBS(3)
#if defined(SSE_ATT)
int_double operator - (const int_double &lhs, const int_double &rhs)
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
#elif defined(_SSE_INTEL)
int_double operator - (const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	__asm{
	mov eax,dword ptr lhs
	movapd XMM0,[eax]
	mov eax,dword ptr rhs
	movapd XMM1,[eax]
	shufpd XMM1,XMM1,1
	addpd XMM0,XMM1
	movapd [temp],XMM0
	}
	return(temp);
}
#else
int_double operator - (const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	temp.left=lhs.left+rhs.right;
	temp.right=lhs.right+rhs.left;
	return(temp);
}
#endif

ATTRIBS(3)
#if defined(SSE_ATT)
int_double operator - (const int_double &lhs, const double &rhs)
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
#elif defined(_SSE_INTEL)
int_double operator - (const int_double &lhs, const double &rhs)
{
	__asm{
	movddup XMM0,rhs
	mov eax,dword ptr [lhs]
	movapd XMM1,[eax]
	addsubpd XMM1,XMM0
	movapd (op_temp),XMM1
	}
	return(op_temp);
}
#else
int_double operator - (const int_double &lhs, const double &rhs)
{
	int_double temp;
	temp.left=lhs.left-rhs;
	temp.right=lhs.right+rhs;
	return(temp);
}
#endif

ATTRIBS(2)
#if defined(SSE_ATT)
int_double operator - (const int_double &lhs)
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
#elif defined(_SSE_INTEL)
int_double operator - (const int_double &lhs)
{
	  __asm{
		  mov eax,dword ptr [lhs]
		  movapd XMM0,[eax]
		  shufpd XMM0,XMM0,1
		  movapd [op_temp],XMM0
	  }
	  return(op_temp);
}
#else
int_double operator - (const int_double &lhs)
{
	  int_double temp;
      temp.left=lhs.right;
      temp.right=lhs.left;
      return(temp);
}
#endif

inline int_double operator -= (int_double &lhs, const int_double &rhs)
{
	return(lhs=lhs-rhs);
}

inline int_double operator -= (int_double &lhs, const double &rhs)
{
	return(lhs=lhs-rhs);
}

// Lambov's version
ATTRIBS(3)
#if defined(SSE_ATT)
int_double operator * (const int_double &lhs, const int_double &rhs)
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

#elif defined (_SSE_INTEL)
int_double operator * (const int_double &lhs, const int_double &rhs)
{
	double t1,t2,rl,rr;
	if(lhs.left>=0) // mask = 11xx
	{
		if(rhs.left>=0) // mask=1111
		{
			__asm{
			mov eax,dword ptr [lhs]
			movapd XMM0,[eax]
			mov eax,dword ptr [rhs]
			movapd XMM1,[eax]
			movapd XMM2,d_zero
			xorpd XMM2,XMM1
			mulpd XMM2,XMM0
			movapd [op_temp],XMM2
			}
			return(op_temp);
		}
		else // mask = 11xx
		{
			if(rhs.right<0) // mask=1101
			{
				__asm{
				mov eax,dword ptr [rhs]
				movapd XMM0,[eax]
				movapd XMM1,[d_neg_zero]
				xorpd XMM0,XMM1
				mov eax,dword ptr [lhs]
				movddup XMM1,8[eax]
				mulpd XMM0,XMM1
				movapd [op_temp],XMM0
				}
				return(op_temp);
			}
			else // mask=1100
			{
				__asm{
				mov eax,dword ptr [lhs]
				movapd XMM0,[eax]
				shufpd XMM0,XMM0,1
				mov eax,dword ptr [rhs]
				movapd XMM1,[eax]
				movapd XMM2,d_neg_neg_zero
				xorpd XMM1,XMM2
				mulpd XMM0,XMM1
				movapd (op_temp),XMM0
				}
				return(op_temp);
			}
		}
	}
	else // mask = 0xxx
	{
		if(lhs.right>=0)  // mask=00xx
		{
			if(rhs.left>=0)  // mask=0011
			{
				__asm{
				mov eax,dword ptr [rhs]
				movapd XMM0,[eax]
				shufpd XMM0,XMM0,1
				mov eax,dword ptr [lhs]
				movapd XMM1,[eax]
				movapd XMM2,d_neg_neg_zero
				xorpd XMM1,XMM2
				mulpd XMM0,XMM1
				movapd (op_temp),XMM0
				}
				return(op_temp);
			}
			else // mask=000x
			{
				if(rhs.right<0) // mask=0001
				{
					__asm{
					mov eax,dword ptr [lhs]
					movddup XMM0,[eax]
					mov eax,dword ptr [rhs]
					movapd XMM1,[eax]
					movapd XMM2,d_neg_zero
					xorpd XMM1,XMM2
					shufpd XMM1,XMM1,1
					mulpd XMM0,XMM1
					movapd (op_temp),XMM0
					}
					return(op_temp);
				}
				else // mask=0000
				{
					__asm{
					mov eax,dword ptr [lhs]
					mov ecx,dword ptr [rhs]
					movapd XMM0,[eax]
					movapd XMM1,[ecx]
					movapd XMM2,d_neg_neg_zero
					xorpd XMM1,XMM2
					mulpd XMM0,XMM1
					shufpd XMM0,XMM0,1
					movapd (op_temp),XMM0
					}
					return(op_temp);				
				}
			}
		}
		else // mask=01xx
		{
			if(rhs.left>=0) // mask= 0111
			{
				__asm{
				mov eax,dword ptr [lhs]
				movapd XMM0,[eax]
				movapd XMM1,[d_neg_zero]
				xorpd XMM0,XMM1
				mov eax,dword ptr [rhs]
				movddup XMM1,8[eax]
				mulpd XMM0,XMM1
				movapd [op_temp],XMM0
				}
				return(op_temp);
			}
			else  // mask=010x
			{
				if(rhs.right>=0) // mask = 0100
				{
					__asm{
					mov eax,dword ptr [rhs]
					movddup XMM0,[eax]
					mov eax,dword ptr [lhs]
					movapd XMM1,[eax]
					movapd XMM2,d_neg_zero
					xorpd XMM1,XMM2
					shufpd XMM1,XMM1,1
					mulpd XMM0,XMM1
					movapd (op_temp),XMM0
					}
					return(op_temp);
				} 
				else //mask = 0101
				{
					rl=-rhs.left;
					rr=-rhs.right;
					t1=lhs.left*rr;
					t2=lhs.right*rl;
					if(t1<=t2)
						op_temp.left=t1;
					else
						op_temp.left=t2;
					t1=lhs.left*rl;
					t2=lhs.right*rr;
					if(t1>=t2)
						op_temp.right=t2;
					else
						op_temp.right=t1;
					return(op_temp);
				}
			}
		}
	}
}
#else
int_double operator * (const int_double &lhs, const int_double &rhs)
{
	double t1,t2,rl,rr;
	int_double temp;
	if(lhs.left>=0) // mask = 11xx
	{
		if(rhs.left>=0) // mask=1111
		{
			temp.left=lhs.left*rhs.left;
			temp.right=lhs.right*(-rhs.right);
			return(temp);
		}
		else // mask = 11xx
		{
			if(rhs.right<0) // mask=1101
			{

				temp.left=lhs.right*(-rhs.left);
				temp.right=lhs.right*(-rhs.right);
				return(temp);
			}
			else // mask=1100
			{
				temp.left=lhs.right*(-rhs.left);
				temp.right=lhs.left*rhs.right;
				return(temp);
			}
		}
	}
	else // mask = 0xxx
	{
		if(lhs.right>=0)  // mask=00xx
		{
			if(rhs.left>=0)  // mask=0011
			{
				temp.left=lhs.left*(-rhs.right);
				temp.right=lhs.right*rhs.left;
				return(temp);
			}
			else // mask=000x
			{
				if(rhs.right<0) // mask=0001
				{
					temp.left=lhs.left*(-rhs.right);
					temp.right=lhs.left*(-rhs.left);
					return(temp);
				}
				else // mask=0000
				{
					temp.left=lhs.right*rhs.right;
					temp.right=lhs.left*(-rhs.left);
					return(temp);
				}
			}
		}
		else // mask=01xx
		{
			if(rhs.left>=0) // mask= 0111
			{
				temp.right=lhs.right*rhs.left;
				temp.left=lhs.left*(-rhs.right);
				return(temp);
			}
			else  // mask=010x
			{
				if(rhs.right>=0) // mask = 0100
				{
					temp.left=lhs.right*(-rhs.left);
					temp.right=lhs.left*(-rhs.left);
					return(temp);
				} 
				else //mask = 0101
				{
					rl=-rhs.left;
					rr=-rhs.right;
					t1=lhs.left*rr;
					t2=lhs.right*rl;
					if(t1<=t2)
						temp.left=t1;
					else
						temp.left=t2;
					t1=lhs.left*rl;
					t2=lhs.right*rr;
					if(t1>=t2)
						temp.right=t2;
					else
						temp.right=t1;
					return(temp);
				}
			}
		}
	}
}
#endif

// faster version of times if we know everything is +ve
ATTRIBS(3)
#if defined (SSE_ATT)
int_double times_pos(const int_double &lhs, const int_double &rhs)
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
#elif defined(_SSE_INTEL)
int_double times_pos(const int_double &lhs, const int_double &rhs)
{
			__asm{
			mov eax,dword ptr [lhs]
			movapd XMM0,[eax]
			mov eax,dword ptr [rhs]
			movapd XMM1,[eax]
			movapd XMM2,d_zero
			xorpd XMM2,XMM1
			mulpd XMM2,XMM0
			movapd [op_temp],XMM2
			}
			return(op_temp);
}
#else
int_double times_pos(const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	temp.left=lhs.left*rhs.left;
	temp.right=lhs.right*(-rhs.right);
	return(temp);
}
#endif

inline int_double operator *= (int_double &lhs, const int_double &rhs)
{
	return(lhs=lhs*rhs);
}

ATTRIBS(3)
#if defined(SSE_ATT)
int_double operator * (const int_double &lhs, const double &rhs)
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
/*
	int_double temp;
	if(rhs>=0.0)
	__asm("movddup %2,%%xmm0\n\t"
		"mulpd %1,%%xmm0\n\t"
		"movapd %%xmm0,%0"
		:"=m" (temp)
		:"m" (lhs), "m" (rhs)
		:"xmm0");
	else
		__asm("movddup %2,%%xmm0\n\t" // rhs rhs
		"xorpd %3,%%xmm0\n\t" // -rhs -rhs
		"mulpd %1,%%xmm0\n\t" // left*(-rhs) right*(-rhs)
		"shufpd $1,%%xmm0,%%xmm0\n\t" // right*(-rhs) left*(-rhs)
		"movapd %%xmm0,%0\n\t"
		:"=m" (temp)
		:"m" (lhs), "m" (rhs), "m" (d_neg_zero)
		:"xmm0");
	return(temp);
}
*/
#elif defined(_SSE_INTEL)
int_double operator * (const int_double &lhs, const double &rhs)
{
	if(rhs>=0.0)
		__asm{
			movddup xmm0,rhs
			mulpd xmm0,lhs
			movapd op_temp,xmm0
	};
	else
		__asm{
			movddup xmm0,rhs
			xorpd xmm0,d_neg_zero
			mulpd xmm0,lhs
			shufpd xmm0,xmm0,1
			movapd op_temp,xmm0};
	return(op_temp)
}
#else
int_double operator * (const int_double &lhs, const double &rhs)
{
	return(lhs*int_double(rhs));
}
#endif


inline int_double operator *= (int_double &lhs, const double &rhs)
{
	return(lhs=lhs*rhs);
}

ATTRIBS(3)
#if defined(SSE_ATT)
int_double operator / (const int_double &lhs, const int_double &rhs)
{
		int_double temp;
//		print_int_double_str("",rhs);
		if(contains_zero(rhs))
		{
			printf("Division by interval containing zero. Exiting\n");
			print_int_double_str("",rhs);
			exit(0);
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
#else
int_double operator / (const int_double &lhs, const int_double &rhs)
{
	int_double temp;
	if(lhs.left>=0.0)  // mask=11xx
	{
		if(rhs.left>=0.0)  // mask=1111
		{
			temp.left=lhs.left/(-rhs.right);
			temp.right=lhs.right/rhs.left;
			return(temp);
		}
		else // mask= 110x
		{
			if(rhs.right>=0) // mask = 1100
			{
				temp.left=lhs.right/rhs.right;
				temp.right=lhs.left/(-rhs.left);
				return(temp);
			}
			else // mask = 1101
			{
				printf("Division by interval containing zero.\n");
				exit(0);
			}
		}
	}
	else // mask = 0xxx
	{
		if(lhs.right>=0) // mask = 00xx
		{
			if(rhs.right>=0) // mask = 0000
			{
				temp.left=(-lhs.right)/rhs.left;
				temp.right=lhs.left/rhs.right;
				return(temp);
			}
			else // mask = 00x1
			{
				if(rhs.left>=0) // mask = 0011
				{
					temp.left=lhs.left/rhs.left;
					temp.right=-(lhs.right/rhs.right);
					return(temp);
				}
				else // mask = 0001
				{
					printf("Division by interval containing zero.\n");
					exit(0);
				}
			}
		}
		else // mask = 01xx
		{
			if(rhs.right>=0) // mask = 0100
			{
				temp.left=lhs.right/rhs.right;
				temp.right=lhs.left/(-rhs.right);  // bugged?
				return(temp);
			}
			else // mask = 01x1
			{
				if(rhs.left>=0) // mask = 0111
				{
					temp.left=lhs.left/rhs.left;
					temp.right=lhs.right/rhs.left;
					return(temp);
				}
				else // mask = 0101
				{
					printf("Division by interval containing zero.\n");
					exit(0);
				}
			}
		}
	}
}
#endif


inline int_double operator /= (int_double &lhs, const int_double &rhs)
{
	return(lhs=lhs/rhs);
}

ATTRIBS(3)
#if defined(SSE_ATT)
int_double operator / (const int_double &lhs, const double &rhs)
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
	exit(0);
}
#elif defined(_SSE_INTEL)
int_double operator / (const int_double &lhs, const double &rhs)
{
	exit(0);
}
#else
int_double operator / (const int_double &lhs, const double &rhs)
{

	double t1;
	int_double temp;
	if(rhs>0.0)
	{
		temp.left=lhs.left/rhs;
		temp.right=lhs.right/rhs;
		return(temp);
	}
	if(rhs<0.0)
	{
		t1=-rhs;
		temp.left=lhs.right/t1;
		temp.right=lhs.left/t1;
		return(temp);
	}
	printf("Division by zero in int_double / double. Exiting.\n");
	exit(0);
}
#endif


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


inline int_double exp (const int_double &x)  // nicely increasing
{
	int_double temp;
	temp.left=exp(x.left);
	if(x.right==0.0) // treat this specially because it happens a lot!
		temp.right=-1.0;
	else
	    temp.right=-nextafter(exp(-x.right));
	return(temp);
}

// if right hand endpoint is one, then interval is
// widened unnecessarily. could trap it but......
inline int_double log (const int_double &x) // nicely increasing
{
//	print_int_double_str("log called with",x);
//	print_int_double_str("returning",int_double(log(x.left),nextafter(log(-x.right))));
	return(int_double(log(x.left),nextafter(log(-x.right))));
}


inline void sin_cos(const int_double &x, int_double *sin_x, int_double *cos_x)
{
	double cos_left,cos_right,sin_left,sin_right;
	unsigned char mask=0;

	if(x.right+x.left<=-M_PI_2)
	{
		printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");

		exit(0);
		return;
	}

#ifdef LINUX
	__asm("fldl %2\n\t"
		"fsincos\n\t"
		"fstpl %0\n\t"
		"fstpl %1\n\t"
		:"=m" (cos_left), "=m" (sin_left)
		:"m" (x.left)
		:"st(6)", "st(7)");
	__asm("fldl %2\n\t"
		"fsincos\n\t"
		"fstpl %0\n\t"
		"fchs\n\t"
		"fstpl %1\n\t"
		:"=m" (cos_right), "=m" (sin_right)
		:"m" (x.right)
		:"st(6)", "st(7)");
#else
	{ double _x_=x.left;
	__asm{
		fld (_x_)
			fsincos
			fstp (cos_left)
			fstp (sin_left)
	}
	_x_=x.right;
	__asm{
		fld (_x_)
			fsincos
			fstp (cos_right)
			fchs
			fstp (sin_right)
	}
	}
#endif
/*
	cos_left=cos(x.left);
	cos_right=cos(x.right);
	sin_left=sin(x.left);
	sin_right=-sin(x.right);
*/
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
		printf("Weird error in sin_cos, mask = %d exiting.\n",(int) mask);
		printf("x: ");
		print_int_double(x);
		printf("\ncos_left: %20.18e\ncos_right: %20.18e\n",cos_left,cos_right);
		printf("sin_left: %20.18e\nsin_right: %20.18e\n",sin_left,sin_right);
		exit(0);
	}
}

// do me in asm
ATTRIBS(2)
#if defined(SSE_ATT)
int_double sqr(const int_double &x)
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
#else
int_double sqr(const int_double &x)
{
	int_double res;
	if(x.left<0)
	{
		if(x.right>0)
		{
			res.left=x.right*x.right;
			res.right=x.left*(-x.left);
		}
		else
		{
			res.left=0.0;
			if (x.left<x.right)
				res.right=x.left*(-x.left);
			else
				res.right=x.right*(-x.right);
		}
	}
	else
	{
		res.left=x.left*x.left;
		res.right=x.right*(-x.right);
	}
	return(res);
}
#endif

int_double pow_temp;
inline int_double pow(const int_double &x, const int_double &y)
{
	return(exp(log(x)*y));
}

inline int_double pow(const int_double &x, const double &y)
{
	return(exp(log(x)*y));
}

inline int_double atan2(const int_double &y, const int_double &x)
{
	int_double y_x,temp;
	// check x doesn't span a discontinuity
	//
	if(contains_zero(x))
	{
		print_int_double_str("y ",y);
		print_int_double_str("x ",x);
		printf("Error in atan2, x contains zero. Exiting.\n");
		printf("%20.18e\n%20.18e\n%20.18e\n%20.18e\n",atan2(y.left,x.left),atan2(y.left,-x.right),
			atan2(-y.right,x.left),atan2(-y.right,-x.right));
		print_int_double_str("Pi/2 ",d_pi_2);
		exit(0);
	}
	y_x=y/x;
	temp.left=atan(y_x.left);
	temp.right=nextbefore(atan(y_x.right)); // could be tighter? use nextbefore
	if(x.right>=0) // left half plane
	{
		if(temp.right>0) // theta <0
			temp=temp+d_pi;
		else
			temp=temp-d_pi;
	}
	return(temp);
}

inline int_double sqrt(const int_double &x) // strictly increasing
{
	int_double temp;

	temp.left=sqrt(x.left);
	temp.right=-nextafter(sqrt(-x.right));
	return(temp);
}

ATTRIBS(2)
#if defined(SSE_ATT)
int contains_zero(const int_double &x)
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
#else
int contains_zero(const int_double &x)
{
	return((x.left<=0.0)&&(x.right<=0.0));
}
#endif

inline int contains_zero(const int_complex &z)
{
	return(contains_zero(z.real)&&contains_zero(z.imag));
}

inline void print_interval (const int_double &x)
{
	print_int_double(x);
}

//int_complex c_op_temp,c_op_temp1,c_op_temp2,c_op_temp3;

ATTRIBS(2)
#if defined (SSE_ATT)
int_complex conj(const int_complex &rhs)
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
#else
int_complex conj(const int_complex &rhs)
{
	return(int_complex(rhs.real,-rhs.imag));
}
#endif


ATTRIBS(3)
#if defined (SSE_ATT)
int_complex operator + (const int_complex &lhs, const int_complex &rhs)
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
#else
int_complex operator + (const int_complex &lhs, const int_complex &rhs)
{
	return(int_complex(lhs.real+rhs.real,lhs.imag+rhs.imag));
}
#endif

ATTRIBS(3)
#if defined(SSE_ATT)
int_complex operator + (const int_complex &lhs, const int_double &rhs)
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
#else
int_complex operator + (const int_complex &lhs, const int_double &rhs)
{
	return(int_complex(lhs.real+rhs,lhs.imag));
}
#endif

ATTRIBS(3)
#if defined(SSE_ATT)
int_complex operator + (const int_complex &lhs, const double &rhs)
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
#else
int_complex operator + (const int_complex &lhs, const double &rhs)
{
	return(int_complex(lhs.real+rhs,lhs.imag));
}
#endif

inline int_complex operator += (int_complex &lhs, const int_complex &rhs)
{
	return(lhs=lhs+rhs);
}

inline int_complex operator += (int_complex &lhs, const int_double &rhs)
{
	return(lhs=lhs+rhs);
}

inline int_complex operator += (int_complex &lhs, const double &rhs)
{
	return(lhs=lhs+rhs);
}

ATTRIBS(2)
#if defined(SSE_ATT)
int_complex operator - (const int_complex &lhs)
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
#else
int_complex operator - (const int_complex &lhs)
{
	return(int_complex(-lhs.real,-lhs.imag));
}
#endif

ATTRIBS(3)
#if defined(SSE_ATT)
int_complex operator - (const int_complex &lhs, const int_complex &rhs)
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
#else
int_complex operator - (const int_complex &lhs, const int_complex &rhs)
{
	return(int_complex(lhs.real-rhs.real,lhs.imag-rhs.imag));
}
#endif

ATTRIBS(3)
#if defined(SSE_ATT)
int_complex operator - (const int_complex &lhs, const int_double &rhs)
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
#else
int_complex operator - (const int_complex &lhs, const int_double &rhs)
{
	return(int_complex(lhs.real-rhs,lhs.imag));
}
#endif

ATTRIBS(3)
#if defined(SSE_ATT)
int_complex operator - (const int_complex &lhs, const double &rhs)
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
#else
int_complex operator - (const int_complex &lhs, const double &rhs)
{
	return(int_complex(lhs.real-rhs,lhs.imag));
}
#endif

inline int_complex operator -= (int_complex &lhs, const int_complex &rhs)
{
	return(lhs=lhs-rhs);
}

inline int_complex operator -= (int_complex &lhs, const int_double &rhs)
{
	return(lhs=lhs-rhs);
}

inline int_complex operator -= (int_complex &lhs, const double &rhs)
{
	return(lhs=lhs-rhs);
}

ATTRIBS(3)
#if defined(SSE_ATT)
int_complex operator * (const int_complex &lhs, const int_complex &rhs)
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
#else
int_complex operator * (const int_complex &lhs, const int_complex &rhs)
{
	return(int_complex(lhs.real*rhs.real-lhs.imag*rhs.imag,
		lhs.real*rhs.imag+lhs.imag*rhs.real));
}
#endif

inline int_complex operator * (const int_complex &lhs, const int_double &rhs)
{
	return(int_complex(lhs.real*rhs,lhs.imag*rhs));
}

inline int_complex operator * (const int_complex &lhs, const double &rhs)
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

inline int_complex sqrt(const int_complex &z) // returns sqrt with arg in [-pi/2,pi/2]
{
	int_double theta,mod_z;
	int_complex res;
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


// these power routines work in long double prior to
// calculating the angle mod 2 Pi to improve accuracy
// when x*im(s) is large
const long double two_pi_left=6.283185307179586476L;
const long double two_pi_right=6.283185307179586477L;

inline int_complex pow1 (const double &x, const double &re, const double &im)
{
	return(pow1(int_double(x),re,im));
}

inline int_complex pow1 (const int_double &x, const double &re, const double &im)
{
//	if(sizeof(double)==sizeof(long double)) // long double is synonym for double
		return(pow(x,int_complex(int_double(re),int_double(im))));

	int_double tlnx,lnx=log(x);

	long double tlnx_left=lnx.left*im;
	long double tlnx_right=-lnx.right*im;

	long double n=floor(tlnx_left/two_pi_left); // how many multiples of 2 Pi
	if(n>=0.0)
	{
		tlnx.left=tlnx_left-n*two_pi_left;
		tlnx.right=-nextafter(tlnx_right-n*two_pi_right);
	}
	if(n<0.0)
	{
		tlnx.left=tlnx_left-n*two_pi_right;
		tlnx.right=-nextafter(tlnx_right-n*two_pi_left);
	}

	int_complex exp_itlnx;

	sin_cos(tlnx,&exp_itlnx.imag,&exp_itlnx.real);

	return(exp_itlnx*exp(re*lnx));
}


inline int_double norm(const int_complex &z)
{
	return(sqr(z.real)+sqr(z.imag));
}

inline int_double argument(const int_complex &z)
{
	return(atan2(z.imag,z.real));
}

inline int_complex log(const int_complex &z)
{
	return(int_complex(log(norm(z))/2,argument(z)));
}


#define MAX_BERNOULLI_K (7)

int_double bernoulli[MAX_BERNOULLI_K],h_bernoulli[MAX_BERNOULLI_K];


void set_h_bernoulli()
{
	bernoulli[0]=d_one/6;
	h_bernoulli[0]=bernoulli[0]/2;				// b(2)/2!
	bernoulli[1]=d_one/(-30);				// b(4)
	h_bernoulli[1]=bernoulli[1]/24;				// b(4)/4!
	bernoulli[2]=d_one/42;				// b(6)
	h_bernoulli[2]=bernoulli[2]/720;				// b(6)/6!
	bernoulli[3]=bernoulli[1];					// b(8)
	h_bernoulli[3]=bernoulli[3]/40320;
	bernoulli[4]=int_double(5)/66;				// b(10)
	h_bernoulli[4]=bernoulli[4]/3628800;
	bernoulli[5]=int_double(-691)/2730;			//b(12)
	h_bernoulli[5]=bernoulli[5]/479001600;                  //b(12)/12!
	bernoulli[6]=int_double(7)/6;				//b(14)
	h_bernoulli[6]=bernoulli[6]/479001600;					//b(14)/12!
	h_bernoulli[6]=h_bernoulli[6]/(13*14);                    //b(14)/14!

}

bool h_bernoulli_initialised=false;

#define DEFAULT_N (50)

int_complex hurwitz (const int_complex &s, const int_double &alpha)
{
	unsigned int i,k=MAX_BERNOULLI_K,n=DEFAULT_N;
	int_double err,n_alpha,s_mod=sqrt(norm(s));
	int_complex s1=s,res=c_zero,n_alpha_s,term;
	int_complex s_array[MAX_BERNOULLI_K];

	if(!h_bernoulli_initialised)
	{
		set_h_bernoulli();
		h_bernoulli_initialised=true;

	}

	if(n<-s_mod.right)
		n=(unsigned int) ceil(-s_mod.right);

	n_alpha=alpha+n;
	n_alpha_s=pow(n_alpha,-s+1);

	for(i=0;i<n;i++)
		res=res+pow(alpha+i,-s);

	res=res+n_alpha_s/(s-1);

	n_alpha_s=n_alpha_s/n_alpha; // n_alpha_s=(n+alpha)^(-s)
		
	res=res+n_alpha_s/2;

	s_array[0]=s*n_alpha_s/n_alpha;

	for(i=1;i<k;i++)
	{
		s1=s1+1;
		s_array[i]=s_array[i-1]*s1/n_alpha;
		s1=s1+1;
		s_array[i]=s_array[i]*s1/n_alpha;
	}

	for(i=0;i<k;i++)
	{
		term=s_array[i]*h_bernoulli[i];
		res=res+term;
	}

	err=sqrt(norm(term))/(s.real+2*k-1)*(s_mod+2*k-1);
	if(err.left<=err.right)
		err.right=err.left;
	else
		err.left=err.right;

	res.real=res.real+err;
	res.imag=res.imag+err;
	return(res);
}

//
// by Hare 1997 prop 4.1 with m=7 and |Im(z)|>=10
// 
#define MAX_LN_GAMMA_ERR (1.2e-14)
int_complex ln_gamma_err=int_complex(int_double(-MAX_LN_GAMMA_ERR,MAX_LN_GAMMA_ERR),
									 int_double(-MAX_LN_GAMMA_ERR,MAX_LN_GAMMA_ERR));

//
// by Hare 1997 equation 4.1 with m=7, sec(theta/2)<=2^(0.5) ie Re(z)>0
// used when Im(z)<10
//
#define MAX_LN_GAMMA_ERR_2 (8.3e-14)
int_complex ln_gamma_err_2=int_complex(int_double(-MAX_LN_GAMMA_ERR_2,MAX_LN_GAMMA_ERR_2),
									   int_double(-MAX_LN_GAMMA_ERR_2,MAX_LN_GAMMA_ERR_2));


// assumes Re(z)>0 so that sec(arg(z)/2)<=sqrt(2)
inline int_complex lngamma(const int_complex &z1)
{
	int_complex res,z_2n_1,lns,z=z1;
	int_double err;
	unsigned int n;

	if(!h_bernoulli_initialised)
	{
		set_h_bernoulli();
		h_bernoulli_initialised=true;
	}

	if(z.imag.left>=10.0) // no need to use recurrence
	{
		res=(z-0.5)*log(z)-z+d_ln_two_pi*0.5;
		res=res+int_complex(bernoulli[0],d_zero)/z/2;
		z_2n_1=z;
		n=1;
		while(n<MAX_BERNOULLI_K)
		{
			z_2n_1=z_2n_1*z*z;
			res=res+int_complex(bernoulli[n],d_zero)/z_2n_1/((n+n+2)*(n+n+1));
			n++;
		}
		return(res+ln_gamma_err);
	}
	else
	{
		lns=c_zero;
		for(n=0;n<10;n++)
			lns=lns+log(z+n);
		z=z+10;
		res=(z-0.5)*log(z)-z+d_ln_two_pi*0.5;
		res=res+int_complex(bernoulli[0],d_zero)/z/2;
		z_2n_1=z;
		n=1;
		while(n<MAX_BERNOULLI_K)
		{
			z_2n_1=z_2n_1*z*z;
			res=res+int_complex(bernoulli[n],d_zero)/z_2n_1/((n+n+2)*(n+n+1));
			n++;
		}
		return(res-lns+ln_gamma_err_2);
	}
}

#ifdef LINUX
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
	_fpu_getcw(old_cw);
	 new_cw=old_cw&0x3FF;
	 new_cw=new_cw|0x400;
	_fpu_setcw(new_cw);
	__SSE_getcw(old__SSE_cw);
	new__SSE_cw=old__SSE_cw&0x1FBF; // zeros FTZ(15),RC(14,13) and DAZ(6) bits
	new__SSE_cw=new__SSE_cw|0x2000; // sets RC to Round Down (14=0, 13=1)
	__SSE_setcw(new__SSE_cw);
	printf("Set FPU control word to %X and SSE control word to %X\n",new_cw,new__SSE_cw);
	d_zero_zero.left=0.0;
	d_zero_zero.right=0.0;
	d_zero.left=0.0;
	d_zero.right=_nze;
	d_neg_zero.left=_nze;
	d_neg_zero.right=_nze;
	d_neg_neg_zero.left=_nze;
	d_neg_neg_zero.right=0.0;
	d_two_pi=d_pi*2;
	d_pi_2=d_pi/2;
	d_ln_two_pi=log(d_two_pi);
	c_zero.real=d_zero;
	c_zero.imag=d_zero;
	c_one.real=int_double(1.0);
	c_one.imag=d_zero;
}
#else
// _control87 sets both FPU and SSE control words
// ref:- http://msdn.microsoft.com/en-us/library/e9b52ceh(VS.80).aspx
void _fpu_rndd ()
{
	_control87(_RC_DOWN,_MCW_RC);
	d_zero.left=0.0;
	d_zero.right=_nze;
	d_neg_zero.left=_nze;
	d_neg_zero.right=_nze;
	d_neg_neg_zero.left=_nze;
	d_neg_neg_zero.right=0.0;
	d_two_pi=d_pi*2;
	d_pi_2=d_pi/2;
	d_ln_two_pi=log(d_two_pi);
	c_zero.real=d_zero;
	c_zero.imag=d_zero;
	c_one.real=int_double(1.0);
	c_one.imag=d_zero;
}
#endif


#endif
