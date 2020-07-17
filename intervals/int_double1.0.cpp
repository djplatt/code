/*

File: int_double1.0.cpp

Created: 28th June 2008

Version: <v> = 1.0

Last Modified: 28th June 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: Requires Intel FPU chip

Build instructions: g++ -oint_double<v> int_double<v>.cpp -O2

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "1.0"

#include <iostream>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <ieeefp.h>

#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  : "=m" (cw))

unsigned short int old_cw,new_cw;
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

using namespace std;

class int_double{

public:
  double left;
  double right;
  int_double ();                      // constructors
  int_double (double);
  int_double (double,double);
  int_double operator + (int_double); // operators
  int_double operator * (int_double);
  int_double operator * (double);
  int_double operator *= (int_double);
  int_double operator * (int);
  int_double operator / (int_double);
  int_double operator / (double);
  int_double operator / (int);
};

int_double::int_double()  // left and right are garbage
{
    left=0.0;
    right=0.0;
};

int_double::int_double(double l)  // right is garbage
{
  left=l;
  right=0.0;
};

int_double::int_double(double l,double r)
{
  left=l;
  right=r;
};


int_double int_double::operator+ (const int_double rhs)
{
  int_double temp;
  _fpu_with_rounding(temp.left=left+rhs.left,temp.right=right+rhs.right);
  return(temp);
};

inline int_double int_double::operator* (const int_double rhs)
{
  int_double temp;
  double t1,t2;
  unsigned char mask=0;

  if(left>=0)
    mask=mask|0x8;
  if(right>=0)
    mask=mask|0x4;
  if(rhs.left>=0)
    mask=mask|0x2;
  if(rhs.right>=0)
    mask=mask|0x1;


  //  printf (" 0x%X ",mask);
  switch(mask)
    {

    case 0:
      _fpu_with_rounding(temp.left=right*rhs.right,temp.right=left*rhs.left);
      return(temp);
    case 1: 
      _fpu_with_rounding(temp.left=left*rhs.right,temp.right=left*rhs.left);
      return(temp);
    case 3: 
      _fpu_with_rounding(temp.left=left*rhs.right,temp.right=right*rhs.left);
      return(temp);
    case 4:
      _fpu_with_rounding(temp.left=right*rhs.left,temp.right=left*rhs.left);
      return(temp);
    case 5:
      _fpu_with_rounding(
      {
	t1=left*rhs.right;
	t2=right*rhs.left;
	if(t1<=t2)
	  temp.left=t1;
	else
	  temp.left=t2;},
      {
	t1=left*rhs.left;
	t2=right*rhs.right;
	if(t1>=t2)
	  temp.right=t1;
	else
	  temp.right=t2;});
      return(temp);
    case 7:
      _fpu_with_rounding(temp.left=left*rhs.right,temp.right=right*rhs.right);
      return(temp);
    case 12:
      _fpu_with_rounding(temp.left=right*rhs.left,temp.right=left*rhs.right);
      return(temp);
    case 13:
      _fpu_with_rounding(temp.left=right*rhs.left,temp.right=right*rhs.right);
      return(temp);
    case 15:
      _fpu_with_rounding(temp.left=left*rhs.left,temp.right=right*rhs.right);
      return(temp);
    default: 
      printf("Bad interval passed to *. Mask was %x. Exiting",mask);
      exit(0);
      return(temp);
    };
};
inline int_double int_double::operator*= (const int_double rhs)
{

  double t1,t2;
  unsigned char mask=0;

  if(left>=0)
    mask=mask|0x8;
  if(right>=0)
    mask=mask|0x4;
  if(rhs.left>=0)
    mask=mask|0x2;
  if(rhs.right>=0)
    mask=mask|0x1;


  //  printf (" 0x%X ",mask);
  switch(mask)
    {

    case 0:
      _fpu_with_rounding(left=right*rhs.right,right=left*rhs.left);
      return(*this);
    case 1: 
      _fpu_with_rounding(left=left*rhs.right,right=left*rhs.left);
      return(*this);
    case 3: 
      _fpu_with_rounding(left=left*rhs.right,right=right*rhs.left);
      return(*this);
    case 4:
      _fpu_with_rounding(left=right*rhs.left,right=left*rhs.left);
      return(*this);
    case 5:
      _fpu_with_rounding(
      {
	t1=left*rhs.right;
	t2=right*rhs.left;
	if(t1<=t2)
	  left=t1;
	else
	  left=t2;},
      {
	t1=left*rhs.left;
	t2=right*rhs.right;
	if(t1>=t2)
	  right=t1;
	else
	  right=t2;});
      return(*this);
    case 7:
      _fpu_with_rounding(left=left*rhs.right,right=right*rhs.right);
      return(*this);
    case 12:
      _fpu_with_rounding(left=right*rhs.left,right=left*rhs.right);
      return(*this);
    case 13:
      _fpu_with_rounding(left=right*rhs.left,right=right*rhs.right);
      return(*this);
    case 15:
      _fpu_with_rounding(left=left*rhs.left,right=right*rhs.right);
      return(*this);
    default: 
      printf("Bad interval passed to *. Mask was %x. Exiting",mask);
      exit(0);
      return(*this);
    };
};

int_double int_double::operator* (const double rhs)
{
  int_double temp;
  if(rhs>=0)
    _fpu_with_rounding(temp.left=left*rhs,temp.right=right*rhs)
  else
    _fpu_with_rounding(temp.left=right*rhs,temp.right=left*rhs);
  return(temp);
};

int_double int_double::operator* (const int rhs)
{
  int_double temp;
  if(rhs>=0)
    _fpu_with_rounding(temp.left=left*rhs,temp.right=right*rhs)
  else
    _fpu_with_rounding(temp.left=right*rhs,temp.right=left*rhs);
  return(temp);
};

int_double int_double::operator/ (const int_double rhs)
  /* assumes that intervals are in correct order */
{
  int_double temp;

  return(temp);
};  


int_double int_double::operator/ (const int rhs)
{
  int_double temp;
  if (rhs>0)
    _fpu_with_rounding(temp.left=left/rhs,temp.right=right/rhs)
  else
    if (rhs<0)
      _fpu_with_rounding(temp.left=right/rhs,temp.right=left/rhs)
    else
      {
	printf("Attempt to divide interval by integer 0. Exiting.");
	exit(0);
      };
    
  return(temp);
};

int_double int_double::operator/ (const double rhs)
{
  int_double temp;
  if (rhs>0.0)
    _fpu_with_rounding(temp.left=left/rhs,temp.right=right/rhs)
  else
    if (rhs<0.0)
      _fpu_with_rounding(temp.left=right/rhs,temp.right=left/rhs)
    else
      {
	printf("Attempt to divide interval by double 0.0. Exiting.");
	exit(0);
      };
    
  return(temp);
};



void print_int_double(int_double x)
{
  printf("[%20.18e,%20.18e]",x.left,x.right);
};

int main()
{
  int_double x[6];
  int_double y;

  unsigned int i,j;
  int_double test_vec[100000];

  cout << "Entered main " << endl;
  /*
  x[0]=int_double(-10,-2);
  x[1]=int_double(-5,0);
  x[2]=int_double(-2,2);
  x[3]=int_double(0,7);
  x[4]=int_double(3.2,7.8);

  for(i=0;i<5;i++)
    for(j=i;j<5;j++)
      {
	print_int_double(x[i]);
	cout << endl << " * " << endl;
	print_int_double(x[j]);
	cout << endl << " = " << endl;
	print_int_double(x[i]*x[j]);
	cout << endl << endl;
      };
  */
  y=int_double(0.999999999999999,0.999999999999999);
  /* clever?
  _fpu_with_rounding(
  {for(i=0;i<1000000000;i++)
    y.left=y.left*y.left;},
  {for(i=0;i<1000000000;i++)
    y.right=y.right*y.right;});
  */

  /* dumb? */

  for(i=0;i<100000000;i++)
    y*=y;
  /*
  for (i=0;i<1000;i++)
    for (j=0;j<100000;j++)
      test_vec[j]=y*y;

    print_int_double(test_vec[100]);
    cout << endl;
  */    


  print_int_double(y);

/* _fpgetround isn't in cygwin unix!
  fp_rnd round_mode;

  round_mode=fpgetround();
  cout << round_mode << endl;
  */

return(0);
};


