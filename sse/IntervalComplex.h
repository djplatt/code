/*

File: interval_complex.h

Created: 2008

Version: <v> = 1.0

Last Modified: 3rd september 2008

Dialect: C++

Requires: see build instructions

Implementation notes: Requires Intel FPU chip + SSE

Build instructions:
 g++ -o <target> <target.cpp> interval_complex.cpp IntervalsSSE2.cpp -I. -msse3


By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

//#define VERSION "1.0"

#ifndef INTERVALCOMPLEX
#define INTERVALCOMPLEX

#include "IntervalsSSE2.h"
#include "IntervalSuggested.h"
#include <iostream>
#include <iomanip>
#include <typeinfo>


using namespace RealLib;
using namespace std;

typedef MachineEstimate (*BinaryFuncME)(const MachineEstimate &a, const MachineEstimate &b);
typedef __m128d (*BinaryFuncM128D)(__m128d a, __m128d b);


bool AreSame(const MachineEstimate &a, const MachineEstimate &b)
{
	return _mm_movemask_pd(_mm_cmpeq_pd(a.GetInterval(), b.GetInterval())) == 3;
}

bool TestFunction(BinaryFuncME def, BinaryFuncM128D impl)
{
	const int cnt = 4;
	MachineEstimate a[cnt] =
	{		
		MachineEstimate(1, 2),
		MachineEstimate(3, 4),
		MachineEstimate(-6, -5),
		MachineEstimate(-7, 8)
	};

	for (int i=0;i<cnt;++i)
		for (int j=0;j<cnt;++j) {
			try {
				MachineEstimate v = def(a[i], a[j]);
				MachineEstimate w = impl(a[i].GetInterval(), a[j].GetInterval());
				if (!AreSame(v, w)) 
					return false;
			} catch (RealLibException &e) {
				try {
					MachineEstimate v = def(a[i], a[j]);
					return false;
				} catch (RealLibException &e1) {
					try {
						MachineEstimate w = impl(a[i].GetInterval(), a[j].GetInterval());
						return false;
					} catch (RealLibException &e2) {
						if (typeid(e1) != typeid(e2)) 
							return false;
					}
				}
			}
		}

	return true;
	
}


unsigned long __nze[2]={0,0x80000000}; // -0.0

unsigned long __delta[2]={1,0};       // very small

double _nze = *(double*)__nze;

double _delta= *(double*)__delta;


/*
inline double nextafter (double x)
{
  unsigned long int *i=(unsigned long int *) &x;

//  printf("nextafter called with %20.18e\n",x);
  if(i[1]&0x80000000)   // -ve
    {
      if(x==_nze)               // -0.0
	{
//	  printf("nextafter returning   %20.18e\n",_delta);
	  return(_delta);
	};

      i[0]--;
      if(i[0]==0xffffffff)
	{
	  i[1]--;
	};
//      printf("nextafter returning   %20.18e\n",x);
      return(x);
    };

  if((i[1]&0x7ff00000)==0x7ff00000) // nan or +/-inf
    {
//      printf("nextafter returning   %20.18e\n",x);
      return(x);
    };


  i[0]++;
  if(i[0]==0)
    i[1]++;
//  printf("nextafter returning   %20.18e\n",x);
  return(x);
};
*/

MachineEstimate exp (MachineEstimate x)  // nicely increasing
{
  double temp[2];
  unsigned long int *i;

  x.MachineEstimate::GetInterval(temp);

  // _mm_storeu_pd(temp,x.MachineEstimate::interval);

  temp[1]=exp(temp[1]);
  temp[0]=exp(temp[0]);   // exp(x) >= 0 for all double x

  /*
  i=(unsigned long int *) &temp;  // inline nextafter
  i[0]++;
  if(i[0]==0)
    i[1]++;
  */

  temp[0]=nextafter(temp[0],DBL_MAX);

  return(MachineEstimate(temp[1],temp[0]));
};

MachineEstimate log (MachineEstimate x) // nicely increasing
{
  double temp[2];
  unsigned long int *i;

  //  _mm_storeu_pd(temp,x.MachineEstimate::interval);

  x.MachineEstimate::GetInterval(temp);

  temp[1]=log(temp[1]);
  temp[0]=log(temp[0]);

  temp[0]=nextafter(temp[0],DBL_MAX);
  
  return(MachineEstimate(temp[1],temp[0]));
};


void sin_cos (const MachineEstimate x1, MachineEstimate *sin_x,
MachineEstimate *cos_x)
{
  double x[2],cos_left,cos_right,sin_left,sin_right;
  unsigned char mask;

  x1.MachineEstimate::GetInterval(x);

  //  _mm_storeu_pd(x,x1.MachineEstimate::interval);

  //printf("in sin_cos [%20.18e,%20.18e]\n",x[1],x[0]);

  if(x[0]+x[1]<=-M_PI_2)
    {
      printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");
      exit(0);
    };

  cos_left=cos(x[1]);
  cos_right=cos(x[0]);
  sin_left=sin(x[1]);
  sin_right=sin(x[0]);
  //printf("%24.22e %24.22e %24.22e %24.22e\n",cos_left,cos_right,sin_left,sin_right);
  //printf("%24.22e\n",sin(3.125));
  if(cos_left>=0.0)
    mask=8;
  else
    mask=0;
  if(cos_right>=0.0)
    mask=mask|0x4;
  if(sin_left>=0.0)
    mask=mask|0x2;
  if(sin_right>=0.0)
    mask++;
//  printf("%d\n",(int) mask);
  switch (mask)
    {
    case 0: // all -ve, 3rd quadrant
      sin_x->interval=_mm_set_pd(sin_right,-nextafter(sin_left,DBL_MAX));
      cos_x->interval=_mm_set_pd(cos_left,-nextafter(cos_left,DBL_MAX));
      //sin_x->left=sin_right;
      //sin_x->right=-nextafter(sin_left);
      //cos_x->left=cos_left;
      //cos_x->right=-nextafter(cos_right);
      return;
    case 2:
	printf("2: pathological sin_cos. exiting");
	exit(0);
	/*
      sin_x->left=-nextafter(-sin_right);
      sin_x->right=-nextafter(sin_left);
      cos_x->left=-1.0;
      if(cos_left>=cos_right)
	cos_x->right=-nextafter(cos_left);
      else
	cos_x->right=-nextafter(cos_right);
	*/
      return;
    case 3:
      sin_x->interval=_mm_set_pd(sin_right,-nextafter(sin_left,DBL_MAX));
      cos_x->interval=_mm_set_pd(cos_right,-nextafter(cos_left,DBL_MAX));
      /*
      sin_x->left=sin_right;
      sin_x->right=-nextafter(sin_left);
      cos_x->left=cos_right;
      cos_x->right=-nextafter(cos_left);
      */
      return;
    case 4:
	printf("4: pathological sin_cos. exiting");
	exit(0);
	/*
      sin_x->left=-1.0;
      if(sin_right>=sin_left)
	sin_x->right=-nextafter(sin_right);
      else
	sin_x->right=-nextafter(sin_left);
      cos_x->left=cos_left;
      cos_x->right=-nextafter(cos_right);
	*/
      return;
    case 11:
	printf("11: pathological sin_cos. exiting");
	exit(0);
	/*
      if(sin_left<=sin_right)
	sin_x->left=sin_left;
      else
	sin_x->left=sin_right;
      sin_x->right=-1.0;
      cos_x->left=cos_right;
      cos_x->right=-nextafter(cos_left);
	*/
      return;
    case 12:
      sin_x->interval=_mm_set_pd(sin_left,-nextafter(sin_right,DBL_MAX));
      cos_x->interval=_mm_set_pd(cos_left,-nextafter(cos_right,DBL_MAX));
      /*
      sin_x->left=sin_left;
      sin_x->right=-nextafter(sin_right);
      cos_x->left=cos_left;
      cos_x->right=-nextafter(cos_right);
      */
      return;
    case 13:
      if(cos_left<=cos_right)
	cos_x->interval=_mm_set_pd(cos_left,-1.0);
      else
	cos_x->interval=_mm_set_pd(cos_right,-1.0);
      sin_x->interval=_mm_set_pd(sin_left,-nextafter(sin_right,DBL_MAX));
      /*
      cos_x->right=-1.0;
      if(cos_left<=cos_right)
	cos_x->left=cos_left;
      else
	cos_x->left=cos_right;
      sin_x->left=sin_left;
      sin_x->right=-nextafter(sin_right);
      */
      return;
    case 15:
      sin_x->interval=_mm_set_pd(sin_left,-nextafter(sin_right,DBL_MAX));
      cos_x->interval=_mm_set_pd(cos_right,-nextafter(cos_left,DBL_MAX));
      /*
      sin_x->left=sin_left;
      sin_x->right=-nextafter(sin_right);
      cos_x->left=cos_right;
      cos_x->right=-nextafter(cos_left);
      */
      return;
    default:
      printf("Weird error in sin_cos, mask = %d exiting.\n",(int) mask);
      exit(0);
    };
};


/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Define int_complex */

class int_complex{
public:
  MachineEstimate real;
  MachineEstimate imag;
  int_complex ();                      // constructors
  int_complex (MachineEstimate,MachineEstimate);
  int_complex operator + (int_complex); // operators
  int_complex operator + (MachineEstimate);
  //  int_complex operator + (dcomplex);
  int_complex operator + (double);
  //int_complex operator + (int);
  int_complex operator += (int_complex); // operators
  int_complex operator += (MachineEstimate);
  //  int_complex operator += (dcomplex);
  int_complex operator += (double);
  //int_complex operator += (int);
  int_complex operator - ();
  int_complex operator - (int_complex);
  int_complex operator - (MachineEstimate);
  int_complex operator - (double);
  //int_complex operator - (int);
  int_complex operator -= (int_complex);
  int_complex operator -= (MachineEstimate);
  int_complex operator -= (double);
  //int_complex operator -= (int);
  int_complex operator * (int_complex);
  int_complex operator * (MachineEstimate);
  int_complex operator * (double);
  //int_complex operator * (int);
  int_complex operator *= (int_complex);
  int_complex operator *= (MachineEstimate);
  int_complex operator *= (double);
  //int_complex operator *= (int);
  int_complex operator / (int_complex);
  int_complex operator / (MachineEstimate);
  int_complex operator / (double);
  //int_complex operator / (int);
};

inline int_complex::int_complex()
{
  real=MachineEstimate(0.0,0.0);
  imag=MachineEstimate(0.0,0.0);
};

inline int_complex::int_complex(MachineEstimate re, MachineEstimate im)
{
  real=re;
  imag=im;
}; 

inline int_complex conj(int_complex rhs)
{
  return(int_complex(rhs.real,-rhs.imag));
};

inline int_complex int_complex::operator+ (const int_complex rhs)
{
  int_complex temp(real,imag);
  temp.real+=rhs.real;
  temp.imag+=rhs.imag;
  return(temp);
};

inline int_complex int_complex::operator+ (const MachineEstimate rhs)
{
  int_complex temp(real,imag);
  temp.real+=rhs;
  return(temp);
};

inline int_complex int_complex::operator+ (const double rhs)
{
  int_complex temp(real,imag);
  temp.real+=rhs;
  return(temp);
};
/*
inline int_complex int_complex::operator+ (const int rhs)
{
  int_complex temp(real,imag);
  temp.real+=rhs;
  return(temp);
};
*/
inline int_complex int_complex::operator+= (const int_complex rhs)
{
  //  int_complex temp(real,imag);
  real+=rhs.real;
  imag+=rhs.imag;
  return(*this);
};

inline int_complex int_complex::operator+= (const MachineEstimate rhs)
{
  //  int_complex temp(real,imag);
  real+=rhs;
  return(*this);
};

inline int_complex int_complex::operator+= (const double rhs)
{
  //  int_complex temp(real,imag);
  real+=rhs;
  return(*this);
};
/*
inline int_complex int_complex::operator+= (const int rhs)
{
  //  int_complex temp(real,imag);
  real+=rhs;
  return(*this);
};
*/

inline int_complex int_complex::operator- ()
{
  return(int_complex(-real,-imag));
};

inline int_complex int_complex::operator- (const int_complex rhs)
{
  return(int_complex(real-rhs.real,imag-rhs.imag));
};

inline int_complex int_complex::operator- (const MachineEstimate rhs)
{
  return(int_complex(real-rhs,imag));
};

inline int_complex int_complex::operator- (const double rhs)
{
  return(int_complex(real-rhs,imag));
};

/*
inline int_complex int_complex::operator- (const int rhs)
{
  return(int_complex(real-rhs,imag));
};
*/

inline int_complex int_complex::operator-= (const int_complex rhs)
{
  //  int_complex temp(real,imag);
  real-=rhs.real;
  imag-=rhs.imag;
  return(*this);
};

inline int_complex int_complex::operator-= (const MachineEstimate rhs)
{
  //  int_complex temp(real,imag);
  real-=rhs;
  return(*this);
};

inline int_complex int_complex::operator-= (const double rhs)
{
  //  int_complex temp(real,imag);
  real-=rhs;
  return(*this);
};

/*
inline int_complex int_complex::operator-= (const int rhs)
{
  //  int_complex temp(real,imag);
  real-=rhs;
  return(*this);
};
*/

inline int_complex int_complex::operator* (const int_complex rhs)
{
    MachineEstimate old_r;
//    printf("\nin int_complex * int_complex\n");
//    printf_int_complex(*this);printf(" * ");print_int_complex(rhs);
//    printf("\ngave );print_int_complex(res);printf("\n");
    old_r=real;
    return(int_complex(real*rhs.real-imag*rhs.imag,
		       old_r*rhs.imag+imag*rhs.real));
};

inline int_complex int_complex::operator* (const MachineEstimate rhs)
{
  return(int_complex(real*rhs,imag*rhs));
};

inline int_complex int_complex::operator* (const double rhs)
{
  MachineEstimate r(rhs,rhs);
  return(int_complex(real*r,imag*r));
};

/*
inline int_complex int_complex::operator* (const int rhs)
{
  return(int_complex(real*rhs,imag*rhs));
};
*/

inline int_complex int_complex::operator*= (const int_complex rhs)
{
    MachineEstimate old_real;
    old_real=real;
    real=real*rhs.real-imag*rhs.imag;
    imag=old_real*rhs.imag+imag*rhs.real;
  return(*this);
};

inline int_complex int_complex::operator*= (const MachineEstimate rhs)
{
  real*=rhs;
  imag*=rhs;
  return(*this);
};

inline int_complex int_complex::operator*= (const double rhs)
{
  MachineEstimate r(rhs,rhs);
    real*=r;
    imag*=r;
  return(*this);
};

/*
inline int_complex int_complex::operator*= (const int rhs)
{
  real*=rhs;
  imag*=rhs;
  return(*this);
};
*/

inline int_complex int_complex::operator/ (const int_complex rhs)
{
  return(*this*conj(rhs)/(sq(rhs.real)+sq(rhs.imag)));
};

inline int_complex int_complex::operator/ (const MachineEstimate rhs)
{
  return(int_complex(real/rhs,imag/rhs));
};

inline int_complex int_complex::operator/ (const double rhs)
{
  MachineEstimate r(rhs,rhs);
  return(int_complex(real/r,imag/r));
};

/*
inline int_complex int_complex::operator/ (const int rhs)
{
  return(int_complex(real/rhs,imag/rhs));
};
*/
#define real(z) z.real
#define imag(z) z.imag

/*
inline int_complex exp(int_complex z)
{
  MachineEstimate xs,xc;

  sin_cos(z.imag,&xs,&xc);
  return(int_complex(xc,xs)*exp(z.real));
};
*/


void print_interval(MachineEstimate z)
{
  double d[2];
  _mm_storeu_pd(d,z.MachineEstimate::interval);
  printf("[%20.18e,%20.18e]",d[1],-d[0]);
}

void print_int_complex(int_complex z)
{
  print_interval(z.real);
  printf("+");
  print_interval(z.imag);
  printf("i\n");
};

int_complex pow(MachineEstimate x,int_complex s)
  /* return x^s */
{
  MachineEstimate tlnx,x_to_sigma,tc,ts;

  tlnx=imag(s)*log(x);

  x_to_sigma=exp(real(s)*log(x));

  sin_cos(tlnx,&ts,&tc);

  return(int_complex(x_to_sigma*tc,x_to_sigma*ts));
};




#endif

