/*

File: interval_complex.cpp

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

#include "IntervalsSSE2.h"
#include "IntervalSuggested.h"
#include <iostream>
#include <iomanip>

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


int int_double::operator>= (const int_double rhs)
{
  return(left>=-rhs.right);
};

int int_double::operator> (const int_double rhs)
{
  return(left>(-rhs.right));
};


int_double exp (int_double x)  // nicely increasing
{
  int_double temp;
  unsigned long int *i;
  //if((x.left==0)&&(x.right==0))
  //return(int_double(1,1));
  temp.left=exp(x.left);
  temp.right=exp(-x.right);   // exp(x) >= 0 for all double x

  i=(unsigned long int *) &temp.right;
  i[0]++;
  if(i[0]==0)
    i[1]++;

  temp.right=-temp.right;

  return(temp);
};

int_double log (int_double x) // nicely increasing
{
  int_double temp;
  unsigned long int *i;

  temp.left=log(x.left);
  temp.right=log(-x.right);

  temp.right=nextafter(temp.right);

  temp.right=-temp.right;
  
  return(temp);
};

void sin_cos(const int_double x, int_double *sin_x, int_double *cos_x)
{
  double cos_left,cos_right,sin_left,sin_right;
  unsigned char mask=0;
  
  //if((x.left==0.0)&&(x.right==0.0))
  //{
  //sin_x->left=0.0;
  //sin_x->right=0.0;
  //cos_x->left=1.0;
  //cos_x->right=-1.0;
  //return;
  //};
  

  if(x.right+x.left<=-M_PI_2)
    {
      printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");
      _fpu_restore;
      exit(0);
      return;
    };

  cos_left=cos(x.left);
  cos_right=cos(x.right);
  sin_left=sin(x.left);
  sin_right=sin(-x.right);
//  printf("%24.22e %24.22e %24.22e %24.22e\n",cos_left,cos_right,sin_left,sin_right);
//  printf("%24.22e\n",sin(3.125));
  if(cos_left>=0.0)
    mask=8;
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
      sin_x->left=sin_right;
      sin_x->right=-nextafter(sin_left);
      cos_x->left=cos_left;
      cos_x->right=-nextafter(cos_right);
      return;
    case 2:
	printf("2: pathological sin_cos. exiting");
	exit(0);
      sin_x->left=-nextafter(-sin_right);
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
	printf("4: pathological sin_cos. exiting");
	exit(0);
      sin_x->left=-1.0;
      if(sin_right>=sin_left)
	sin_x->right=-nextafter(sin_right);
      else
	sin_x->right=-nextafter(sin_left);
      cos_x->left=cos_left;
      cos_x->right=-nextafter(cos_right);
      return;
    case 11:
	printf("11: pathological sin_cos. exiting");
	exit(0);
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
      exit(0);
    };
};

*/


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

int main()
{
  MachineEstimate::Computation comp;
  int_complex s;

  s=int_complex(MachineEstimate(0.5,0.5),MachineEstimate(-1,3));
  s=s*s/3;
  print_int_complex(s);

}

