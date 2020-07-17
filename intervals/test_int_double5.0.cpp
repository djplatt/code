/*

File: test_int_double5.0.cpp

Created: 30th August 2008

Version: <v> = 4.0

Last Modified: 30th August 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: Requires Intel FPU chip

Build instructions: expects to be included
                    Dont optimise in GCC. Breaks at -O1


By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

//#define VERSION "5.0"

//#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  : "=m" (cw))

unsigned short int old_cw,new_cw;

// set the fpu to round down (towards -infty) 
void _fpu_rndd ()
{
    _fpu_getcw(old_cw);
    new_cw=old_cw&0x3FF;
    printf("rounding down\n");
    new_cw=new_cw|0x400;
    _fpu_setcw(new_cw);
};

// restore the fpu control register from memory
#define _fpu_restore {_fpu_setcw(old_cw);}


using namespace std;


unsigned long __nze[2]={0,0x80000000}; // -0.0

unsigned long __delta[2]={1,0};       // very small

double _nze = *(double*)__nze;

double _delta= *(double*)__delta;



double nextafter (double x)
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

class int_double{
public:
  double left;
  double right;
  int_double ();                      // constructors
  int_double (double);
  int_double (double,double);
  int_double operator + (int_double); // operators
  int_double operator + (double);
  //int_double operator + (int);
  int_double operator += (int_double);
  int_double operator += (double);
  //int_double operator += (int);
  int_double operator - ();           // unary - int_double
  int_double operator - (int_double); // int_double - int_double
  int_double operator - (double);     
  //int_double operator - (int);
  int_double operator -= (int_double);
  int_double operator -= (double);
  //int_double operator -= (int);
  int_double operator * (int_double);
  int_double operator * (double);
  //int_double operator * (int);
  int_double operator *= (int_double);
  int_double operator *= (double);
  //int_double operator *= (int);
  int_double operator / (int_double);
  int_double operator / (double);
  //int_double operator / (int);
  int_double operator /= (int_double);
  int_double operator /= (double);
  //int_double operator /= (int);
  int operator >= (int_double);
  int operator > (int_double);
};


void print_int_double(int_double x)
{
  printf("[%20.18e,%20.18e]",x.left,-x.right);
};


int_double::int_double()
{
  left=0.0;
  right=0.0;
};


int_double::int_double(double l) 
{
  left=l;
  right=-l;
};


int_double::int_double(double l,double r)
{
    if(l>r)
    {
	printf("Error constructing int_double, right<left. Exiting.\n");
        exit(0);
    };
    left=l;
    right=-r;
};


int_double int_double::operator+ (const int_double rhs)
{
  int_double temp;
  temp.left=left+rhs.left;
  temp.right=right+rhs.right;
  return(temp);
};

int_double int_double::operator+ (const double rhs)
{
  int_double temp;
  temp.left=left+rhs;
  temp.right=right-rhs;
  return(temp);
};

int_double int_double::operator+= (const double rhs)
{
  left+=rhs;
  right-=rhs;
  return(*this);

};

int_double int_double::operator- (const int_double rhs)
{
  int_double temp;
  temp.left=left+rhs.right;
  temp.right=right+rhs.left;
  return(temp);
};

int_double int_double::operator* (const int_double rhs)
{
  int_double temp;
  unsigned int *i;
  double t1,t2,rl,rr;
  unsigned char mask=0;
  /*
  if(left>=0)
    mask=12;
  else
    if(right<0)
      mask=4;
  if(rhs.left>=0)
    mask+=3;
  else
    if(rhs.right<0)
      mask++;
  */
  if(left>=0) // mask = 11xx
    {
      if(rhs.left>=0) // mask=1111
	{
	  printf("multiplying \n");
	  print_int_double(*this);
	  printf("\nand\n");
	  print_int_double(rhs);
	  temp.left=left*rhs.left;
	  temp.right=right*(-rhs.right);
	  printf("\nreturning ");
	  print_int_double(temp);
	  printf("\n");
	  return(temp);
	}
      else // mask = 11xx
	{
	  if(rhs.right<0) // mask=1101
	    {
	      temp.left=right*(-rhs.left);
	      temp.right=right*(-rhs.right);
	      return(temp);
	    }
	  else // mask=1100
	    {
	      temp.left=right*(-rhs.left);
	      temp.right=left*rhs.right;
	      return(temp);
	    }
	}
    }
  else // mask = 0xxx
    {
      if(right>=0)  // mask=00xx
	{
	  if(rhs.left>=0)  // mask=0011
	    {
	      temp.left=left*(-rhs.right);
	      temp.right=right*rhs.left;
	      return(temp);
	    }
	  else // mask=000x
	    {
	      if(rhs.right<0) // mask=0001
		{
		  temp.left=left*(-rhs.right);
		  temp.right=left*(-rhs.left);
		  return(temp);
		}
	      else // mask=0000
		{
		  temp.left=right*rhs.right;
		  temp.right=left*(-rhs.left);
		  return(temp);
		}
	    }
	}
      else // mask=01xx
	{
	  if(rhs.left>=0) // mask= 0111
	    {
	      temp.left=left*(-rhs.right);
	      temp.right=right*temp.left;
	      return(temp);
	    }
	  else  // mask=010x
	    {
	      if(rhs.right>=0) // mask = 0100
		{
		  temp.left=right*(-rhs.left);
		  temp.right=left*(-rhs.left);
		  return(temp);
		} 
	      else //mask = 0101
		{
		  rl=-rhs.left;
		  rr=-rhs.right;
		  t1=left*rr;
		  t2=right*rl;
		  if(t1<=t2)
		    temp.left=t1;
		  else
		    temp.left=t2;
		  t1=left*rl;
		  t2=right*rr;
		  if(t1>=t2)
		    temp.right=t2;
		  else
		    temp.right=t1;
		  return(temp);
		}
	    }
	}
    }
};


int_double int_double::operator/ (const double rhs)
{
  double t1;
  int_double temp;
  //  printf("in int_double / double\n");
  if(rhs>0.0)
    {
      temp.left=left/rhs;
      temp.right=right/rhs;
      return(temp);
    };
  if(rhs<0.0)
    {
      t1=-rhs;
      temp.left=right/t1;
      temp.right=left/t1;
      return(temp);
    };
  printf("Division by zero in int_double / double. Exiting.\n");
  exit(0);
  return(temp);
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
  //int_complex operator + (int);
  int_complex operator += (int_complex); // operators
  int_complex operator += (int_double);
  //  int_complex operator += (dcomplex);
  int_complex operator += (double);
  //int_complex operator += (int);
  int_complex operator - ();
  int_complex operator - (int_complex);
  int_complex operator - (int_double);
  int_complex operator - (double);
  //int_complex operator - (int);
  int_complex operator -= (int_complex);
  int_complex operator -= (int_double);
  int_complex operator -= (double);
  //int_complex operator -= (int);
  int_complex operator * (int_complex);
  int_complex operator * (int_double);
  int_complex operator * (double);
  //int_complex operator * (int);
  int_complex operator *= (int_complex);
  int_complex operator *= (int_double);
  int_complex operator *= (double);
  //int_complex operator *= (int);
  int_complex operator / (int_complex);
  int_complex operator / (int_double);
  int_complex operator / (double);
  //int_complex operator / (int);
};

void print_int_complex(int_complex z)
{
  print_int_double(z.real);
  printf("+");
  print_int_double(z.imag);
  printf("i\n");
};

int_complex::int_complex()
{
  real=int_double(0.0,0.0);
  imag=int_double(0.0,0.0);
};

int_complex::int_complex(int_double re, int_double im)
{
  real=re;
  imag=im;
}; 

int_complex int_complex::operator+ (const double rhs)
{
  int_complex temp(real,imag);
  temp.real+=rhs;
  return(temp);
};

int_complex int_complex::operator* (const int_complex rhs)
{
    int_double old_r;
//    printf("\nin int_complex * int_complex\n");
//    printf_int_complex(*this);printf(" * ");print_int_complex(rhs);
//    printf("\ngave );print_int_complex(res);printf("\n");
    old_r=real;
    return(int_complex(real*rhs.real-imag*rhs.imag,
		       old_r*rhs.imag+imag*rhs.real));
};

int_complex int_complex::operator/ (const double rhs)
{
  return(int_complex(real/rhs,imag/rhs));
};


#define real(z) z.real
#define imag(z) z.imag


