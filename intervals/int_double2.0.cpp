/*

File: int_double2.0.cpp

Created: 28th June 2008

Version: <v> = 2.0

Last Modified: 26th August 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: Requires Intel FPU chip

Build instructions: g++ -oint_double<v> int_double<v>.cpp


By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

//#define VERSION "2.0"

//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  : "=m" (cw))

unsigned short int old_cw,new_cw;
/*
#define _fpu_with_rounding(id,iu) {_fpu_getcw(old_cw);\
 new_cw=old_cw&0x3FF;\
 new_cw=new_cw|0x400;\
 _fpu_setcw(new_cw);\
 id;\
 new_cw=new_cw^0xC00;\
 _fpu_setcw(new_cw);\
 iu;\
 _fpu_setcw(old_cw);\
}
*/

// set the fpu to round down (towards -infty) 
#define _fpu_rndd {_fpu_getcw(old_cw);\
 new_cw=old_cw&0x3FF;\
 new_cw=new_cw|0x400;\
 _fpu_setcw(new_cw);}

// restore the fpu control register from memory
#define _fpu_restore {_fpu_setcw(old_cw);}

using namespace std;

unsigned long __nze[2]={0,0x80000000}; // -0.0

unsigned long __delta[2]={1,0};       // very small

double _nze = *(double*)__nze;

double _delta= *(double*)__delta;



inline double nextafter (double x)
{
  unsigned long int *i=(unsigned long int *) &x;


  if(i[1]&0x80000000)   // -ve
    {

      if(x==_nze)               // -0.0
	return(_delta);

      i[0]--;
      //printf("i[0]=%x\n",i[0]);
      if(i[0]==0xffffffff)
	{
	  i[1]--;
	  //printf("i[1]=%x\n",i[1]);
	};
      return(x);
    };

  if((i[1]&0x7ff00000)==0x7ff00000) // nan or +/-inf
    return(x);


  i[0]++;
  if(i[0]==0)
    i[1]++;
  return(x);
};


//#define nextafter(x) x

/*
This was an attempt at being cute but was slower!

union ieee {double x;int i[2];};

inline void nextafterb (ieee *val)
{

 if((val->i[1]&0x7ff00000)==0x7ff00000) // nan or +/-inf
    return;


  if(val->x==_nze)
    {
      val->x=_delta;
      return;
    };

  if(val->i[1]&0x80000000)   // -ve
    {
      val->i[0]--;
      if(val->i[0]==0xffffffff)
	val->i[1]--;
      return;
    };

  val->i[0]++;
  if(val->i[0]==0)
    val->i[1]++;
  return;
}
*/

class int_double{
public:
  double left;
  double right;
  int_double ();                      // constructors
  int_double (double);
  int_double (double,double);
  int_double operator + (int_double); // operators
  int_double operator + (double);
  int_double operator + (int);
  int_double operator += (int_double);
  int_double operator += (double);
  int_double operator += (int);
  int_double operator - ();           // unary - int_double
  int_double operator - (int_double); // int_double - int_double
  int_double operator - (double);     
  int_double operator - (int);
  int_double operator -= (int_double);
  int_double operator -= (double);
  int_double operator -= (int);
  int_double operator * (int_double);
  int_double operator * (double);
  int_double operator * (int);
  int_double operator *= (int_double);
  int_double operator *= (double);
  int_double operator *= (int);
  int_double operator / (int_double);
  int_double operator / (double);
  int_double operator / (int);
  int_double operator /= (int_double);
  int_double operator /= (double);
  int_double operator /= (int);
  int operator >= (int_double);
};


void print_int_double(int_double x)
{
  printf("[%20.18e,%20.18e]",x.left,x.right);
};



inline int_double::int_double()  // left and right are garbage
{
  left=0.0;
  right=0.0;
};


inline int_double::int_double(double l) 
{
  left=l;
  right=l;
};


inline int_double::int_double(double l,double r)
{
    if(l>r)
      {
        printf("Error constructing int_double, right<left. Exiting.\n");
        exit(0);_fpu_restore;
      };
  left=l;
  right=r;
};


inline int_double int_double::operator+ (const int_double rhs)
{
  int_double temp;
  temp.left=left+rhs.left;
  temp.right=nextafter(right+rhs.right);
  return(temp);
};

inline int_double int_double::operator+ (const double rhs)
{
  int_double temp;
  temp.left=left+rhs;
  temp.right=nextafter(right+rhs);
  return(temp);
};

inline int_double int_double::operator+ (const int rhs)
{
  int_double temp;
  temp.left=left+rhs;
  temp.right=nextafter(right+rhs);
  return(temp);
};

inline int_double int_double::operator+= (const int_double rhs)
{
  left=left+rhs.left;
  right=nextafter(right+rhs.right);
  return(*this);
};

inline int_double int_double::operator+= (const double rhs)
{
  left=left+rhs;
  right=nextafter(right+rhs);
  return(*this);
};

inline int_double int_double::operator+= (const int rhs)
{
  left=left+rhs;
  right=nextafter(right+rhs);
  return(*this);
};

inline int_double int_double::operator- ()
{
  return(int_double(-right,-left));
};

inline int_double int_double::operator- (const int_double rhs)
{
  int_double temp;
  temp.left=left-rhs.right;
  temp.right=right-rhs.left;
  
  temp.right=nextafter(temp.right);
  
  return(temp);
};

inline int_double int_double::operator- (const double rhs)
{
  int_double temp;
  temp.left=left-rhs;
  temp.right=nextafter(right-rhs);
  return(temp);
};

inline int_double int_double::operator- (const int rhs)
{
  int_double temp;
  temp.left=left-rhs;
  temp.right=nextafter(right-rhs);
  return(temp);
};

inline int_double int_double::operator-= (const int_double rhs)
{
  left=left-rhs.right;
  right=nextafter(right-rhs.left);
  return(*this);
};

inline int_double int_double::operator-= (const double rhs)
{
  left-=rhs;
  right-=rhs;
  right=nextafter(right);
  return(*this);
};

inline int_double int_double::operator-= (const int rhs)
{
  left-=rhs;
  right-=rhs;
  right=nextafter(right);
  return(*this);
};

inline int_double int_double::operator* (const int_double rhs)
{
  int_double temp;
  double t1,t2;
  unsigned char mask=0;

  if(left>=0)
    mask=12;
  else
    if(right>0)
      mask=4;
  if(rhs.left>=0)
    mask+=3;
  else
    if(rhs.right>0)
      mask++;

  switch(mask)
    {

    case 0:  // all -ve
      temp.left=right*rhs.right;
      temp.right=nextafter(left*rhs.left);
      return(temp);
    case 1: 
      temp.left=left*rhs.right;
      temp.right=nextafter(left*rhs.left);
      return(temp);
    case 3: 
      temp.left=left*rhs.right;
      temp.right=nextafter(right*rhs.left);
      return(temp);
    case 4:
      temp.left=right*rhs.left;
      temp.right=nextafter(left*rhs.left);
      return(temp);
    case 5:
      t1=left*rhs.right;
      t2=right*rhs.left;
      if(t1<=t2)
	temp.left=t1;
      else
	temp.left=t2;
      t1=left*rhs.left;
      t2=right*rhs.right;
	if(t1>=t2)
	  temp.right=nextafter(t1);
	else
	  temp.right=nextafter(t2);
      return(temp);
    case 7:
      temp.left=left*rhs.right;
      temp.right=nextafter(right*rhs.right);
      return(temp);
    case 12:
      temp.left=right*rhs.left;
      temp.right=nextafter(left*rhs.right);
      return(temp);
    case 13:
      temp.left=right*rhs.left;
      temp.right=nextafter(right*rhs.right);
      return(temp);
    case 15:
      temp.left=left*rhs.left;
      temp.right=nextafter(right*rhs.right);
      return(temp);
    default: 
      printf("Bad interval passed to *. Mask was %x. Exiting",mask);
      exit(0);_fpu_restore;
      return(temp);
    };
};

inline int_double int_double::operator*= (const int_double rhs)
{

  double t1,t2,t3;
  unsigned char mask=0;


  if(left>=0)
    mask=12;
  else
    {
      if(right>0)
	mask=4;
    };


  if(rhs.left>=0)
    mask+=3;
  else
    {
      if(rhs.right>0)
	mask++;
    };

  switch(mask)
    {
    case 0:
      t1=left;
      left=right*rhs.right;
      right=nextafter(t1*rhs.left);
      return(*this);
    case 1: 
      t1=left;
      left*=rhs.right;
      right=nextafter(t1*rhs.left);
      return(*this);
    case 3: 
      left*=rhs.right;
      right*=rhs.left;
      right=nextafter(right);
      return(*this);
    case 4:
      right=nextafter(left*rhs.left);
      left=right*rhs.left;
      return(*this);
    case 5:
      t3=left;
      t1=left*rhs.right;
      t2=right*rhs.left;
      if(t1<=t2)
	left=t1;
      else
	left=t2;
      t1=t3*rhs.left;
      t2=right*rhs.right;
      if(t1>=t2)
	right=nextafter(t1);
      else
	right=nextafter(t2);
      return(*this);
    case 7:
      left*=rhs.right;
      right*=rhs.right;
      right=nextafter(right);
      return(*this);
    case 12:
      t1=left;
      left=right*rhs.left;
      right=nextafter(t1*rhs.right);
      return(*this);
    case 13:
      left=right*rhs.left;
      right*=rhs.right;
      right=nextafter(right);
      return(*this);
    case 15:
      left*=rhs.left;
      right*=rhs.right;
      right=nextafter(right);
      return(*this);
    default: 
      printf("Bad interval passed to *. Mask was %x. Exiting",mask);
      exit(0);_fpu_restore;
      return(*this);
    };
};


inline int_double int_double::operator* (const double rhs)
{
  double t1,t2;
  int_double temp;
  t1=left*rhs;
  t2=right*rhs;
  if(t1>=t2)
    {
      temp.left=t2;
      temp.right=nextafter(t1);
      return(temp);
    };
  temp.left=t1;
  temp.right=nextafter(t2);
  return(temp);
};

inline int_double int_double::operator* (const int rhs)
{
  double t1,t2;
  int_double temp;
  t1=left*rhs;
  t2=right*rhs;
  if(t1>=t2)
    {
      temp.left=t2;
      temp.right=nextafter(t1);
      return(temp);
    };
  temp.left=t1;
  temp.right=nextafter(t2);
  return(temp);
};

inline int_double int_double::operator*= (const double rhs)
{
  double t1,t2;
  t1=left*rhs;
  t2=right*rhs;
  if(t1>=t2)
    {
      left=t2;
      right=nextafter(t1);
      return(*this);
    };
  left=t1;
  right=nextafter(t2);
  return(*this);
};

inline int_double int_double::operator*= (const int rhs)
{
  double t1,t2;
  t1=left*rhs;
  t2=right*rhs;
  if(t1>=t2)
    {
      left=t2;
      right=nextafter(t1);
      return(*this);
    };
  left=t1;
  right=nextafter(t2);
  return(*this);
};


inline int_double int_double::operator/ (const int_double rhs)
  /* assumes that intervals are in correct order */
  /* doesnt catch interval divisor with 1/2 endpoints = 0 */
{
  int_double temp;
  unsigned char mask=0;


  if(left>=0)
    mask=12;
  else
    if(right>0)
      mask=4;
  if(rhs.left>=0)
    mask+=3;
  else
    if(rhs.right>0)
      mask++;
  /*
  if(left>=0)
    mask=mask|0x8;
  if(right>=0)
    mask=mask|0x4;
  if(rhs.left>=0)
    mask=mask|0x2;
  if(rhs.right>=0)
    mask=mask|0x1;
  */

  //  printf (" 0x%X ",mask);

  switch(mask)
    {

    case 0:
      temp.left=right/rhs.left;
      temp.right=nextafter(left/rhs.right);
      return(temp);
    case 1:
      printf("Division by interval containing zero.\n");
      exit(0);_fpu_restore;
      return(temp);
    case 3: 
      temp.left=left/rhs.left;
      temp.right=nextafter(right/rhs.right);
      return(temp);
    case 4:
      temp.left=right/rhs.right;
      temp.right=nextafter(left/rhs.right);
      return(temp);
    case 5:
      printf("Division by interval containing zero.\n");
      _fpu_restore;exit(0);
      return(temp);
    case 7:
      temp.left=left/rhs.left;
      temp.right=nextafter(right/rhs.left);
      return(temp);
    case 12:
      temp.left=right/rhs.right;
      temp.right=nextafter(left/rhs.left);
      return(temp);
    case 13:
      printf("Division by interval containing zero.\n");
      _fpu_restore;exit(0);
      return(temp);
    case 15:
      temp.left=left/rhs.right;
      temp.right=nextafter(right/rhs.left);
      return(temp);
    default: 
      printf("Bad interval passed to *. Mask was %x. Exiting",mask);
      exit(0);_fpu_restore;
      return(temp);
    };
};  


inline int_double int_double::operator/ (const double rhs)
{
  double t1,t2;
  int_double temp;
  t1=left/rhs;
  t2=right/rhs;
  if(t1>t2)
    {
      temp.left=t2;
      temp.right=nextafter(t1);
      return(temp);
    };
  temp.left=t1;
  temp.right=nextafter(t2);
  return(temp);
};

inline int_double int_double::operator/ (const int rhs)
{
  double t1,t2;
  int_double temp;
  t1=left/rhs;
  t2=right/rhs;
  if(t1>t2)
    {
      temp.left=t2;
      temp.right=nextafter(t1);
      return(temp);
    };
  temp.left=t1;
  temp.right=nextafter(t2);
  return(temp);
};

inline int_double int_double::operator/= (const int_double rhs)
  /* assumes that intervals are in correct order */
{
  unsigned char mask=0;
  double t1;

  if(left>=0)
    mask=12;
  else
    if(right>0)
      mask=4;
  if(rhs.left>=0)
    mask+=3;
  else
    if(rhs.right>0)
      mask++;
  /*
  if(left>=0)
    mask=mask|0x8;
  if(right>=0)
    mask=mask|0x4;
  if(rhs.left>=0)
    mask=mask|0x2;
  if(rhs.right>=0)
    mask=mask|0x1;
  */

  //  printf (" 0x%X ",mask);

  switch(mask)
    {

    case 0:
      t1=left;
      left=right/rhs.left;
      right=nextafter(t1/rhs.right);
      return(*this);
    case 1:
      printf("Division by interval containing zero.\n");
      exit(0);_fpu_restore;
      return(*this);
    case 3: 
      left/=rhs.left;
      right/=rhs.right;
      right=nextafter(right);
      return(*this);
    case 4:
      t1=left;
      left=right/rhs.right;
      right=nextafter(t1/rhs.right);
      return(*this);
    case 5:
      printf("Division by interval containing zero.\n");
      exit(0);_fpu_restore;
      return(*this);
    case 7:
      left/=rhs.left;
      right/=rhs.left;
      right=nextafter(right);
      return(*this);
    case 12:
      t1=left;
      left=right/rhs.right;
      right=nextafter(t1/rhs.left);
      return(*this);
    case 13:
      printf("Division by interval containing zero.\n");
      exit(0);_fpu_restore;
      return(*this);
    case 15:
      left/=rhs.right;
      right/=rhs.left;
      right=nextafter(right);
      return(*this);
    default: 
      printf("Bad interval passed to *. Mask was %x. Exiting",mask);
      exit(0);_fpu_restore;
      return(*this);
    };
};  

inline int_double int_double::operator/= (const double rhs)
{
  double t1,t2;
  t1=left/rhs;
  t2=right/rhs;
  if(t1>=t2)
    {
      left=t2;
      right=nextafter(t1);
      return(*this);
    };
  left=t1;
  right=nextafter(t2);
  return(*this);
};

inline int_double int_double::operator/= (const int rhs)
{
  double t1,t2;
  t1=left/rhs;
  t2=right/rhs;
  if(t1>=t2)
    {
      left=t2;
      right=nextafter(t1);
      return(*this);
    };
  left=t1;
  right=nextafter(t2);
  return(*this);
};

int int_double::operator>= (const int_double rhs)
{
  return(left>=rhs.right);
};

int_double exp (int_double x)  // nicely increasing
{
  int_double temp;
  temp.left=exp(x.left);
  temp.right=nextafter(exp(x.right));
  return(temp);
};

int_double log (int_double x) // nicely increasing
{
  int_double temp;

  temp.left=log(x.left);
  temp.right=nextafter(log(x.right));
  return(temp);
};

void sin_cos(int_double x, int_double *sin_x, int_double *cos_x)
{
  double cos_left,cos_right,sin_left,sin_right;
  unsigned char mask=0;

if(x.right-x.left>=M_PI_2)
  {
    printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");
    exit(0);_fpu_restore;
    return;
  };

  cos_left=cos(x.left);
  cos_right=cos(x.right);
  sin_left=sin(x.left);
  sin_right=sin(x.right);
  if(cos_left>=0.0)
    mask=mask|0x8;
  if(cos_right>=0.0)
    mask=mask|0x4;
  if(sin_left>=0.0)
    mask=mask|0x2;
  if(sin_right>=0.0)
    mask=mask|0x1;
  //printf("%d\n",(int) mask);
  switch (mask)
    {
    case 0:
      sin_x->left=sin_right;
      sin_x->right=nextafter(sin_left);
      cos_x->left=cos_left;
      cos_x->right=nextafter(cos_right);
      return;
    case 2:
      sin_x->left=sin_right;
      sin_x->right=nextafter(sin_left);
      cos_x->left=-1.0;
      if(cos_left>=cos_right)
	cos_x->right=nextafter(cos_left);
      else
	cos_x->right=nextafter(cos_right);
      return;
    case 3:
      sin_x->left=sin_right;
      sin_x->right=nextafter(sin_left);
      cos_x->left=cos_right;
      cos_x->right=nextafter(cos_left);
      return;
    case 4:
      sin_x->left=-1.0;
      if(sin_right>=sin_left)
	sin_x->right=nextafter(sin_right);
      else
	sin_x->right=nextafter(sin_left);
      cos_x->left=cos_left;
      cos_x->right=nextafter(cos_right);
      return;
    case 11:
      if(sin_left<=sin_right)
	sin_x->left=sin_left;
      else
	sin_x->left=sin_right;
      sin_x->right=1.0;
      cos_x->left=cos_right;
      cos_x->right=nextafter(cos_left);
      return;
    case 12:
      sin_x->left=sin_left;
      sin_x->right=nextafter(sin_right);
      cos_x->left=cos_left;
      cos_x->right=nextafter(cos_right);
      return;
    case 13:
      sin_x->right=1.0;
      if(sin_left>=sin_right)
	sin_x->left=sin_right;
      else
	sin_x->left=sin_left;
      cos_x->left=cos_left;
      cos_x->right=nextafter(cos_right);
      return;
    case 15:
      sin_x->left=sin_left;
      sin_x->right=nextafter(sin_right);
      cos_x->left=cos_right;
      cos_x->right=nextafter(cos_left);
      return;
    default:
      printf("Weird error in sin_cos, mask = %d exiting.\n",(int) mask);
      exit(0);_fpu_restore;
    };
};

int_double sqr(int_double x)
{
  int_double res;
  if(x.left<0)
    {
      if(x.right<0)
	{
	  res.left=x.right*x.right;
	  res.right=nextafter(x.left*x.left);
	}
      else
	{
	  res.left=0.0;
	  if (fabs(x.left)>x.right)
	    res.right=nextafter(x.left*x.left);
	  else
	    res.right=nextafter(x.right*x.right);
	}
    }
  else
    {
      res.left=x.left*x.left;
      res.right=nextafter(x.right*x.right);
    };
  return(res);
};

/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* Define int_complex */

class int_complex{
public:
  int_double real;
  int_double imag;
  int_complex ();                      // constructors
  int_complex (int_double,int_double);
  int_complex operator + (int_complex); // operators
  int_complex operator + (int_double);
  //  int_complex operator + (dcomplex);
  int_complex operator + (double);
  int_complex operator + (int);
  int_complex operator += (int_complex); // operators
  int_complex operator += (int_double);
  //  int_complex operator += (dcomplex);
  int_complex operator += (double);
  int_complex operator += (int);
  int_complex operator - ();
  int_complex operator - (int_complex);
  int_complex operator - (int_double);
  int_complex operator - (double);
  int_complex operator - (int);
  int_complex operator -= (int_complex);
  int_complex operator -= (int_double);
  int_complex operator -= (double);
  int_complex operator -= (int);
  int_complex operator * (int_complex);
  int_complex operator * (int_double);
  int_complex operator * (double);
  int_complex operator * (int);
  int_complex operator *= (int_complex);
  int_complex operator *= (int_double);
  int_complex operator *= (double);
  int_complex operator *= (int);
  int_complex operator / (int_complex);
  int_complex operator / (int_double);
  int_complex operator / (double);
  int_complex operator / (int);
};

inline int_complex::int_complex()
{
  real=int_double(0.0,0.0);
  imag=int_double(0.0,0.0);
};

inline int_complex::int_complex(int_double re, int_double im)
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

inline int_complex int_complex::operator+ (const int_double rhs)
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

inline int_complex int_complex::operator+ (const int rhs)
{
  int_complex temp(real,imag);
  temp.real+=rhs;
  return(temp);
};

inline int_complex int_complex::operator+= (const int_complex rhs)
{
  //  int_complex temp(real,imag);
  real+=rhs.real;
  imag+=rhs.imag;
  return(*this);
};

inline int_complex int_complex::operator+= (const int_double rhs)
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

inline int_complex int_complex::operator+= (const int rhs)
{
  //  int_complex temp(real,imag);
  real+=rhs;
  return(*this);
};

inline int_complex int_complex::operator- ()
{
  return(int_complex(-real,-imag));
};

inline int_complex int_complex::operator- (const int_complex rhs)
{
  return(int_complex(real-rhs.real,imag-rhs.imag));
};

inline int_complex int_complex::operator- (const int_double rhs)
{
  return(int_complex(real-rhs,imag));
};

inline int_complex int_complex::operator- (const double rhs)
{
  return(int_complex(real-rhs,imag));
};

inline int_complex int_complex::operator- (const int rhs)
{
  return(int_complex(real-rhs,imag));
};


inline int_complex int_complex::operator-= (const int_complex rhs)
{
  //  int_complex temp(real,imag);
  real-=rhs.real;
  imag-=rhs.imag;
  return(*this);
};

inline int_complex int_complex::operator-= (const int_double rhs)
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

inline int_complex int_complex::operator-= (const int rhs)
{
  //  int_complex temp(real,imag);
  real-=rhs;
  return(*this);
};

inline int_complex int_complex::operator* (const int_complex rhs)
{
  return(int_complex(real*rhs.real-imag*rhs.imag,
		     real*rhs.imag+imag*rhs.real));
};

inline int_complex int_complex::operator* (const int_double rhs)
{
  return(int_complex(real*rhs,imag*rhs));
};

inline int_complex int_complex::operator* (const double rhs)
{
  return(int_complex(real*rhs,imag*rhs));
};

inline int_complex int_complex::operator* (const int rhs)
{
  return(int_complex(real*rhs,imag*rhs));
};

inline int_complex int_complex::operator*= (const int_complex rhs)
{
  int_double r=real;
  real=r*rhs.real-imag*rhs.imag;
  imag=r*rhs.imag+imag*rhs.real;
  return(*this);
};

inline int_complex int_complex::operator*= (const int_double rhs)
{
  real*=rhs;
  imag*=rhs;
  return(*this);
};

inline int_complex int_complex::operator*= (const double rhs)
{
  real*=rhs;
  imag*=rhs;
  return(*this);
};

inline int_complex int_complex::operator*= (const int rhs)
{
  real*=rhs;
  imag*=rhs;
  return(*this);
};

inline int_complex int_complex::operator/ (const int_complex rhs)
{
  return(*this*conj(rhs)/(sqr(rhs.real)+sqr(rhs.imag)));
};

inline int_complex int_complex::operator/ (const int_double rhs)
{
  return(int_complex(real/rhs,imag/rhs));
};

inline int_complex int_complex::operator/ (const double rhs)
{
  return(int_complex(real/rhs,imag/rhs));
};

inline int_complex int_complex::operator/ (const int rhs)
{
  return(int_complex(real/rhs,imag/rhs));
};

#define real(z) z.real
#define imag(z) z.imag

inline int_complex exp(int_complex z)
{
  int_double xs,xc;

  sin_cos(z.imag,&xs,&xc);
  return(int_complex(xc,xs)*exp(z.real));
};

void print_int_complex(int_complex z)
{
  print_int_double(z.real);
  printf("+");
  print_int_double(z.imag);
  printf("i");
};

/*
int main()
{
  int_complex z(int_double(0.5,0.5),int_double(5,5));
  _fpu_rndd;

  print_int_complex(exp(z));

  _fpu_restore;
  return(0);
};
*/
