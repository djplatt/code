// 22/3/9 v 8.0 added use of nextbefore
#ifndef INT_DOUBLE8
#define INT_DOUBLE8
#define LINUX
#ifndef LINUX
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES 
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>


using namespace std;


#ifdef LINUX
unsigned short int old_cw,new_cw;
#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  : "=m" (cw))
#define _fpu_rndd {_fpu_getcw(old_cw);\
 new_cw=old_cw&0x3FF;\
 new_cw=new_cw|0x400;\
 printf("Setting FPU control word to %X\n",new_cw);\
 _fpu_setcw(new_cw);\
}
#define _aligned_malloc(a,b) malloc(a)

#else
// _control87 sets both FPU and SSE control words
// ref:- http://msdn.microsoft.com/en-us/library/e9b52ceh(VS.80).aspx
#define _fpu_rndd _control87(_RC_DOWN,_MCW_RC)
#define _fpu_rndu _control87(_RC_UP,_MCW_RC)
#endif

inline double nextafter(double);
inline double nextbefore(double);


class int_double{
public:
	double left;
	double right;

  inline int_double ()                      // constructors
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


//  inline int_double (const  double);
//  inline int_double (const  double,const  double);
  inline int_double operator + (const int_double ); // operators
  inline int_double operator + (const double);
  inline int_double operator += (const int_double );
  inline int_double operator += (const double);
  inline int_double operator - (const int_double ); // int_double - int_double
  inline int_double operator - (const double);
  inline int_double operator - ()
  {
	  int_double temp;
      temp.left=right;
      temp.right=left;
      return(temp);
  };
  inline int_double operator -= (const int_double );
  inline int_double operator -= (const double);
  inline int_double operator * (const int_double );
  inline int_double operator * (const double);
  inline int_double operator *= (const int_double );
  inline int_double operator *= (const double);
  inline int_double operator / (const int_double );
  inline int_double operator / (const double);
//  inline int_double operator % (int_double );
  inline int_double operator /= (const int_double );
  inline int_double operator /= (const double);
  inline int operator >= (const int_double );
  inline int operator > (const int_double );
  inline int operator < (const int_double);
};

extern inline int_double times_pos(const int_double, const int_double);

extern inline int_double exp (const int_double );

extern inline int_double log (const int_double );

extern inline void sin_cos(const int_double, int_double *, int_double *);

extern inline int_double sqr(const int_double );

extern inline int_double pow(const int_double , const int_double );

extern inline int_double pow(const int_double , const double);

extern inline int_double atan(const int_double);

extern inline int_double sqrt(const int_double );

extern inline int rel_error(int_double );

extern inline int abs_error(int_double );

extern inline bool contains_zero (const int_double );

extern int_double times_mod_2pi (const int_double, const double);


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
	inline int_complex ();                      // constructors
	inline int_complex (const int_double ,const int_double );
	inline int_complex operator + (const int_complex ); // operators
	inline int_complex operator + (const int_double );
	inline int_complex operator + (const double);
	inline int_complex operator += (const int_complex ); // operators
	inline int_complex operator += (const int_double );
	inline int_complex operator += (const double);
	inline int_complex operator - ()
	{
		int_complex temp;
		temp.real=-real;
		temp.imag=-imag;
		return(temp);
	}
	inline int_complex operator - (const int_complex );
	inline int_complex operator - (const int_double );
	inline int_complex operator - (const double);
	inline int_complex operator -= (const int_complex );
	inline int_complex operator -= (const int_double );
	inline int_complex operator -= (const double);
	inline int_complex operator * (const int_complex );
	inline int_complex operator * (const int_double );
	inline int_complex operator * (const double);
	inline int_complex operator *= (const int_complex );
	inline int_complex operator *= (const int_double );
	inline int_complex operator *= (const double);
	inline int_complex operator / (const int_complex );
	inline int_complex operator / (const int_double );
	inline int_complex operator / (const double);
	inline int_complex operator /= (const int_complex );
	inline int_complex operator /= (const int_double );
	inline int_complex operator /= (const double);
};

extern inline int_complex exp(const int_complex );

extern inline void print_int_complex(const int_complex );

extern inline int_complex pow (int_double, int_complex);

extern inline int_complex pow (const double, const int_complex );

extern inline int_double norm(const int_complex );

extern inline int_complex sqrt(const int_complex );

extern inline void print_interval (const int_double );

extern inline int_complex conj (const int_complex );

extern inline int_double argument (const int_complex);

extern inline int_complex lngamma (const int_complex);

extern inline int_complex hurwitz(const int_complex, const int_double);

extern inline int_complex log (const int_complex);

extern inline bool contains_zero (const int_complex );

extern inline int_double real (const int_complex z)
{
	return(z.real);
}

extern inline int_double imag (const int_complex z)
{
	return(z.imag);
}


static int_double d_zero=int_double(0.0,0.0);
static int_double d_one=int_double(1.0,1.0);
static int_complex c_zero=int_complex(d_zero,d_zero);
static int_complex c_one=int_complex(d_one,d_zero);


static const int_double d_half=int_double(0.5,0.5);

static int _i_pi[2]={1413754135,1074340347};   // this is what double pi less a bit looks like

static int _i_pi2[2]={1413754137,1074340347};  // this is double pi plus a bit

static double *_d_pi=(double *)&_i_pi;

static double *_d_pi2=(double *)&_i_pi2;

static int_double d_pi=int_double(_d_pi[0],_d_pi2[0]);

static int_double d_two_pi=d_pi*2;

static int_double d_ln_two_pi=log(d_two_pi);

static int_double d_pi_2=d_pi/2;

unsigned long __nze[2]={0,0x80000000}; // -0.0

unsigned long __delta[2]={1,0};       // very small

unsigned long __delta_minus[2]={1,0x80000000};

double _nze = *(double*)__nze;

double _delta= *(double*)__delta;

double _delta_minus= *(double*)__delta_minus;


// move a double prec float towards +infinity by
// the smallest delta possible
inline double nextafter (double x)
{
  unsigned long int *i=(unsigned long int *) &x;

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
  unsigned long int *i=(unsigned long int *) &x;

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

void print_int_double(const int_double x)
{
  printf("[%20.18e,%20.18e]",x.left,-x.right);
};



inline int rel_error(int_double x)
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

inline int abs_error(int_double x)
{
	double aerr;
	aerr=x.left+x.right;
	if(aerr==0.0)
		return(0);
	if(aerr<0.0)
		return((int)log10(-aerr));
	return((int)log10(aerr));
}




inline int_double int_double::operator+ (const int_double rhs)
{
	int_double temp;
	temp.left=left+rhs.left;
	temp.right=right+rhs.right;
	return(temp);
}

inline int_double int_double::operator+ (const double rhs)
{
	int_double temp;
	temp.left=left+rhs;
	temp.right=right-rhs;
	return(temp);
}

inline int_double int_double::operator+= (const int_double rhs)
{
  left+=rhs.left;
  right+=rhs.right;
  return(*this);
}

inline int_double int_double::operator+= (const double rhs)
{
  left+=rhs;
  right-=rhs;
  return(*this);
}


inline int_double int_double::operator- (const int_double rhs)
{
  int_double temp;
  temp.left=left+rhs.right;
  temp.right=right+rhs.left;
  return(temp);
}

inline int_double int_double::operator- (const double rhs)
{
  int_double temp;
  temp.left=left-rhs;
  temp.right=right+rhs;
  return(temp);
}

inline int_double int_double::operator-= (const int_double rhs)
{
  left+=rhs.right;
  right+=rhs.left;
  return(*this);
}

inline int_double int_double::operator-= (const double rhs)
{
  left-=rhs;
  right+=rhs;
  return(*this);
}

inline int_double int_double::operator* (const int_double rhs)
{
	int_double temp;
	double t1,t2,rl,rr;
	if(left>=0) // mask = 11xx
	{
		if(rhs.left>=0) // mask=1111
		{
			temp.left=left*rhs.left;
			temp.right=right*(-rhs.right);
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

				temp.right=right*rhs.left;
				temp.left=left*(-rhs.right);
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
}

// faster version of times if we know everything is +ve
inline int_double times_pos(const int_double lhs, const int_double rhs)
{
	int_double res;
	res.left=lhs.left*rhs.left;
	res.right=lhs.right*(-rhs.right);
	return(res);
}

#define times_pos_equals(lhs,rhs){lhs.left*=rhs.left;lhs.right*=-rhs.right;}


inline int_double int_double::operator*= (const int_double rhs)
{
	double t1,t2,t3,rl,rr;
	if(left>=0) // mask = 11xx
	{
		if(rhs.left>=0) // mask=1111
		{
			left*=rhs.left;
			right*=(-rhs.right);
			return(*this);
		}
		else // mask = 11xx
		{
			if(rhs.right<0) // mask=1101
			{
				left=right*(-rhs.left);
				right*=(-rhs.right);
				return(*this);
			}
			else // mask=1100
			{
				t1=left;
				left=right*(-rhs.left);
				right=t1*rhs.right;
				return(*this);
			}
		}
	}
	else // mask = 0xxx
	{
		if(right>=0)  // mask=00xx
		{
			if(rhs.left>=0)  // mask=0011
			{
				left*=(-rhs.right);
				right*=rhs.left;
				return(*this);
			}
			else // mask=000x
			{
				if(rhs.right<0) // mask=0001
				{
					right=left*(-rhs.left);
					left*=(-rhs.right);

					return(*this);
				}
				else // mask=0000
				{
					t1=left;
					left=right*rhs.right;
					right=t1*(-rhs.left);
					return(*this);
				}
			}
		}
		else // mask=01xx
		{
			if(rhs.left>=0) // mask= 0111
			{

				right*=rhs.left;
				left*=(-rhs.right);
				return(*this);
			}
			else  // mask=010x
			{
				if(rhs.right>=0) // mask = 0100
				{
					t1=left;
					left=right*(-rhs.left);
					right=t1*(-rhs.left);
					return(*this);
				} 
				else //mask = 0101
				{
					rl=-rhs.left;
					rr=-rhs.right;
					t1=left*rr;
					t2=right*rl;
					t3=left;
					if(t1<=t2)
						left=t1;
					else
						left=t2;
					t1=t3*rl;
					t2=right*rr;
					if(t1>=t2)
						right=t2;
					else
						right=t1;
					return(*this);
				}
			}
		}
	}
}

inline int_double int_double::operator* (const double rhs)
{
	int_double temp;
	double t1;

	if(rhs>=0.0)
	{
		temp.left=left*rhs;
		temp.right=right*rhs;
		return(temp);
	};
	t1=-rhs;
	temp.left=right*t1;
	temp.right=left*t1;
	return(temp);
}

inline int_double int_double::operator*= (const double rhs)
{
	double t1,t2;
	if(rhs>=0.0)
	{
		left*=rhs;
		right*=rhs;
		return(*this);
	};
	t2=-rhs;
	t1=right*t2;
	right=left*t2;
	left=t1;
	return(*this);
}

inline int_double int_double::operator/ (const int_double rhs)
{
	int_double temp;

	if(left>=0.0)  // mask=11xx
	{
		if(rhs.left>0.0)  // mask=1111
		{
			temp.left=left/(-rhs.right);
			temp.right=right/rhs.left;
			return(temp);
		}
		else // mask= 110x
		{
			if(rhs.right>0) // mask = 1100
			{
				temp.left=right/rhs.right;
				temp.right=left/(-rhs.left);
				return(temp);
			}
			else // mask = 1101
			{
				printf("Division by interval containing zero.\n");
				exit(0);
				return(temp);
			}
		}
	}
	else // mask = 0xxx
	{
		if(right>=0) // mask = 00xx
		{
			if(rhs.right>0) // mask = 0000
			{
				temp.left=(-right)/rhs.left;
				temp.right=left/rhs.right;
				return(temp);
			}
			else // mask = 00x1
			{
				if(rhs.left>0) // mask = 0011
				{
					temp.left=left/rhs.left;
					temp.right=-(right/rhs.right);
					return(temp);
				}
				else // mask = 0001
				{
					printf("Division by interval containing zero.\n");
					exit(0);
					return(temp);
				}
			}
		}
		else // mask = 01xx
		{
			if(rhs.right>0) // mask = 0100
			{
				temp.left=right/rhs.right;
				temp.right=left/rhs.right;
				return(temp);
			}
			else // mask = 01x1
			{
				if(rhs.left>0) // mask = 0111
				{
					temp.left=left/rhs.left;
					temp.right=right/rhs.left;
					return(temp);
				}
				else // mask = 0101
				{
					printf("Division by interval containing zero.\n");
					exit(0);
					return(temp);
				}
			}
		}
	}
}


inline int_double int_double::operator/= (const int_double rhs)
{
	double t1;  
	if(left>=0.0)  // mask=11xx
	{
		if(rhs.left>0.0)  // mask=1111
		{
			left/=(-rhs.right);
			right/=rhs.left;
			return(*this);
		}
		else // mask= 110x
		{
			if(rhs.right>0) // mask = 1100
			{
				t1=left;
				left=right/rhs.right;
				right=t1/(-rhs.left);
				return(*this);
			}
			else // mask = 1101
			{
				printf("Division by interval containing zero.\n");
				exit(0);
				return(*this);
			}
		}
	}
	else // mask = 0xxx
	{
		if(right>=0) // mask = 00xx
		{
			if(rhs.right>0) // mask = 0000
			{
				t1=left;
				left=(-right)/rhs.left;
				right=t1/rhs.right;
				return(*this);
			}
			else // mask = 00x1
			{
				if(rhs.left>0) // mask = 0011
				{
					left/=rhs.left;
					right/=rhs.right;
					right=-right;
					return(*this);
				}
				else // mask = 0001
				{
					printf("Division by interval containing zero.\n");
					exit(0);
					return(*this);
				}
			}
		}
		else // mask = 01xx
		{
			if(rhs.right>0) // mask = 0100
			{
				t1=left;
				left=right/rhs.right;
				right=left/rhs.right;
				return(*this);
			}
			else // mask = 01x1
			{
				if(rhs.left>0) // mask = 0111
				{
					left/=rhs.left;
					right/=rhs.left;
					return(*this);
				}
				else // mask = 0101
				{
					printf("Division by interval containing zero.\n");
					exit(0);
					return(*this);
				}
			}
		}
	}
}  

inline int_double int_double::operator/ (const double rhs)
{
	double t1;
	int_double temp;
	//  printf("in int_double / double\n");
	if(rhs>0.0)
	{
		temp.left=left/rhs;
		temp.right=right/rhs;
		return(temp);
	}
	if(rhs<0.0)
	{
		t1=-rhs;
		temp.left=right/t1;
		temp.right=left/t1;
		return(temp);
	}
	printf("Division by zero in int_double / double. Exiting.\n");
	exit(0);
	return(temp);
}

inline int_double int_double::operator/= (const double rhs)
{
	double t1,t2;
	//  printf("in int_double /= double\n");
	if(rhs>0.0)
	{
		left/=rhs;
		right/=rhs;
		return(*this);
	}
	if(rhs<0.0)
	{
		t1=-rhs;
		t2=right/t1;
		right=left/t1;
		left=t2;
		return(*this);
	}
	printf("Division by zero in int_double /= double. Exiting.\n");
	exit(0);
	return(*this);
}

inline int int_double::operator>= (const int_double rhs)
{
	return(left>=(-rhs.right));
}

inline int int_double::operator> (const int_double rhs)
{
	return(left>(-rhs.right));
}

inline int int_double::operator< (const int_double rhs)
{
	return((-right)<rhs.left);
}


inline int_double exp (const int_double x)  // nicely increasing
{
	int_double temp;
//	unsigned long int *i;
	temp.left=exp(x.left);
	if(x.right==0.0) // treat this specially because it happens a lot!
		temp.right=-1.0;
	else
	{
/* doesn't work under -O3 with GCC
		temp.right=exp(-x.right);   // exp(x) >= 0 for all double x

		// very quick and dirty nextafter(+ve double)
		i=(unsigned long int *) &temp.right;
		i[0]++;
		if(i[0]==0)
			i[1]++;

		temp.right=-temp.right;
*/
	    temp.right=-nextafter(exp(-x.right));
	}
	return(temp);
}

// if right hand endpoint is one, then interval is
// widened unnecessarily. could trap it but......
inline int_double log (const int_double x) // nicely increasing
{
	return(int_double(log(x.left),nextafter(log(-x.right))));
}


inline extern void sin_cos(const int_double x, int_double *sin_x, int_double *cos_x)
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
	cos_left=cos(x.left);
	cos_right=cos(x.right);
	sin_left=sin(x.left);
	sin_right=-sin(x.right);
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
		printf("\ncos_left: %20.18e\ncos_right: %20.18e\n",cos_left,cos_right);		printf("sin_left: %20.18e\nsin_right: %20.18e\n",sin_left,sin_right);
		exit(0);
	}
}

inline extern int_double sqr(const int_double x)
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

inline extern int_double pow(const int_double x, const int_double y)
{
	return(exp(log(x)*y));
}

inline extern int_double pow(const int_double x, const double y)
{
	return(exp(log(x)*y));
}

inline extern int_double atan2(int_double y, int_double x)
{
	int_double y_x,temp;
	// check x doesn't span a discontinuity
	//
	if(contains_zero(x))
	{
		printf("Error in atan2, x contains zero. Exiting.\n");
		exit(0);
	}

	y_x=y/x;
	temp.left=atan(y_x.left);
	temp.right=nextbefore(atan(y_x.right)); // could be tighter? use nextbefore
	if(x.right>=0) // left half plane
	{
		if(temp.right>0) // theta <0
			temp+=d_pi;
		else
			temp-=d_pi;
	}
	return(temp);
}


extern inline int_double sqrt(const int_double x) // strictly increasing
{
	int_double temp;

	temp.left=sqrt(x.left);
	temp.right=-nextafter(sqrt(-x.right));
	return(temp);
}

extern inline bool contains_zero(int_double x)
{
	return((x.left<=0.0)&&(x.right<=0.0));
}

extern inline bool contains_zero(int_complex z)
{
	return(contains_zero(z.real)&&contains_zero(z.imag));
}

extern inline void print_interval (const int_double x)
{
	print_int_double(x);
}

// we don't initialise by default
inline int_complex::int_complex()
{
}

inline int_complex::int_complex(const int_double re, const int_double im)
{
	real=re;
	imag=im;
} 

extern inline int_complex conj(const int_complex rhs)
{
	int_double im=rhs.imag;
	return(int_complex(rhs.real,-im));
}

inline int_complex int_complex::operator+ (const int_complex rhs)
{
	return(int_complex(real+rhs.real,imag+rhs.imag));
}

inline int_complex int_complex::operator+ (const int_double rhs)
{
	return(int_complex(real+rhs,imag));
}

inline int_complex int_complex::operator+ (const double rhs)
{
	return(int_complex(real+rhs,imag));
}

inline int_complex int_complex::operator+= (const int_complex rhs)
{
	real+=rhs.real;
	imag+=rhs.imag;
	return(*this);
}

inline int_complex int_complex::operator+= (const int_double rhs)
{
	real+=rhs;
	return(*this);
}

inline int_complex int_complex::operator+= (const double rhs)
{
	real+=rhs;
	return(*this);
}

inline int_complex int_complex::operator- (const int_complex rhs)
{
	return(int_complex(real-rhs.real,imag-rhs.imag));
}

inline int_complex int_complex::operator- (const int_double rhs)
{
	return(int_complex(real-rhs,imag));
}

inline int_complex int_complex::operator- (const double rhs)
{
	return(int_complex(real-rhs,imag));
}

inline int_complex int_complex::operator-= (const int_complex rhs)
{
	real-=rhs.real;
	imag-=rhs.imag;
	return(*this);
}

inline int_complex int_complex::operator-= (const int_double rhs)
{
	real-=rhs;
	return(*this);
}

inline int_complex int_complex::operator-= (const double rhs)
{
	real-=rhs;
	return(*this);
}

inline int_complex int_complex::operator* (const int_complex rhs)
{
	return(int_complex(real*rhs.real-imag*rhs.imag,
		real*rhs.imag+imag*rhs.real));
}

inline int_complex int_complex::operator* (const int_double rhs)
{
	return(int_complex(real*rhs,imag*rhs));
}

inline int_complex int_complex::operator* (const double rhs)
{
	return(int_complex(real*rhs,imag*rhs));
}

inline int_complex int_complex::operator*= (const int_complex rhs)
{
	int_double old_real;
	old_real=real;
	real*=rhs.real;
	real-=imag*rhs.imag;
	imag=old_real*rhs.imag+imag*rhs.real;
	return(*this);
}

inline int_complex int_complex::operator*= (const int_double rhs)
{
	real*=rhs;
	imag*=rhs;
	return(*this);
}

inline int_complex int_complex::operator*= (const double rhs)
{
	real*=rhs;
	imag*=rhs;
	return(*this);
}

inline int_complex int_complex::operator/ (const int_complex rhs)
{
	return(*this*conj(rhs)/(sqr(rhs.real)+sqr(rhs.imag)));
}

inline int_complex int_complex::operator/ (const int_double rhs)
{
	return(int_complex(real/rhs,imag/rhs));
}

inline int_complex int_complex::operator/ (const double rhs)
{
	return(int_complex(real/rhs,imag/rhs));
}

inline int_complex int_complex::operator/= (const int_complex rhs)
{
	int_double den=sqr(rhs.real)+sqr(rhs.imag);

	*this*=conj(*this);
	*this/=den;
	return(*this);
}

inline int_complex int_complex::operator/= (const int_double rhs)
{
	real/=rhs;
	imag/=rhs;
	return(*this);
}

inline int_complex int_complex::operator/= (const double rhs)
{
	real/=rhs;
	imag/=rhs;
	return(*this);
}

extern inline int_complex exp(const int_complex z)
{
	int_double xs,xc;

	sin_cos(z.imag,&xs,&xc);
	return(int_complex(xc,xs)*exp(z.real));
}

extern inline int_complex sqrt(const int_complex z) // returns sqrt with arg in [-pi/2,pi/2]
{
	int_double theta,mod_z;
	int_complex res,z1;
	z1=z;
	mod_z=pow(norm(z1),0.25);
	theta=atan2(z1.imag,z1.real)/2.0;
	sin_cos(theta,&res.imag,&res.real);
	res*=mod_z;
	return(res);
}

extern inline void print_int_complex(const int_complex z)
{
	print_int_double(z.real);
	printf("+");
	print_int_double(z.imag);
	printf("i");
}

extern inline int_complex pow (const int_double x, const int_complex s)
{
	int_double lnx;
     	int_complex z;
//	int rts,rtc;
//	unsigned int n_2_pi;

	lnx=log(x);
	sin_cos(lnx*s.imag,&z.imag,&z.real);
	return(z*exp(lnx*s.real));
}

extern inline int_complex pow (const double x, const int_complex s)
{
	return(pow(int_double(x),s));
}

extern inline int_double norm(const int_complex z)
{
	return(sqr(z.real)+sqr(z.imag));
}

extern inline int_double argument(const int_complex z)
{
	return(atan2(z.imag,z.real));
}

extern inline int_complex log(int_complex z)
{
	int_double r,arg;

	return(int_complex(log(norm(z))/2,argument(z)));
}


#define MAX_BERNOULLI_K (7)

int_double bernoulli[MAX_BERNOULLI_K],h_bernoulli[MAX_BERNOULLI_K];


void set_h_bernoulli()
{
	bernoulli[0]=int_double(1)/6;
	h_bernoulli[0]=bernoulli[0]/2;				// b(2)/2!
	bernoulli[1]=int_double(-1)/30;				// b(4)
	h_bernoulli[1]=bernoulli[1]/24;				// b(4)/4!
	bernoulli[2]=int_double(1)/42;				// b(6)
	h_bernoulli[2]=bernoulli[2]/720;				// b(6)/6!
	bernoulli[3]=int_double(-1)/30;				// b(8)
	h_bernoulli[3]=bernoulli[3]/40320;
	bernoulli[4]=int_double(5)/66;				// b(10)
	h_bernoulli[4]=bernoulli[4]/3628800;
	bernoulli[5]=int_double(-691)/2730;			//b(12)
	h_bernoulli[5]=bernoulli[5]/479001600;                  //b(12)/12!
	bernoulli[6]=int_double(7)/6;				//b(14)
	h_bernoulli[6]=bernoulli[6]/479001600;					//b(14)/12!
	h_bernoulli[6]/=(13*14);                    //b(14)/14!

}

bool h_bernoulli_initialised=false;

#define DEFAULT_N (50)

int_complex hurwitz (int_complex s, int_double alpha)
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
		res+=pow(alpha+i,-s);

	res+=n_alpha_s/(s-1);

	n_alpha_s/=n_alpha; // n_alpha_s=(n+alpha)^(-s)
		
	res+=n_alpha_s/2;

	s_array[0]=s*n_alpha_s/n_alpha;

	for(i=1;i<k;i++)
	{
		s1=s1+1;
		s_array[i]=s_array[i-1]*s1/n_alpha;
		s1=s1+1;
		s_array[i]*=s1/n_alpha;
	}

	for(i=0;i<k;i++)
	{
		term=s_array[i]*h_bernoulli[i];
		res+=term;
	}

	err=sqrt(norm(term))/(s.real+2*k-1)*(s_mod+2*k-1);
	if(err.left<=err.right)
		err.right=err.left;
	else
		err.left=err.right;

	res.real+=err;
	res.imag+=err;
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
extern inline int_complex lngamma(int_complex z)
{
	int_complex res,z_2n_1,lns;
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
		res+=int_complex(bernoulli[0],d_zero)/z/2;
		z_2n_1=z;
		n=1;
		while(n<MAX_BERNOULLI_K)
		{
			z_2n_1*=z*z;
			res+=int_complex(bernoulli[n],d_zero)/z_2n_1/((n+n+2)*(n+n+1));
			n++;
		}
		return(res+ln_gamma_err);
	}
	else
	{
		lns=c_zero;
		for(n=0;n<10;n++)
			lns+=log(z+n);
		z+=10;
		res=(z-0.5)*log(z)-z+d_ln_two_pi*0.5;
		res+=int_complex(bernoulli[0],d_zero)/z/2;
		z_2n_1=z;
		n=1;
		while(n<MAX_BERNOULLI_K)
		{
			z_2n_1*=z*z;
			res+=int_complex(bernoulli[n],d_zero)/z_2n_1/((n+n+2)*(n+n+1));
			n++;
		}
		return(res-lns+ln_gamma_err_2);
	}
}


#define nextafter(x,y) nextafter(x)
#define nextbefore(x,y) nextbefore(x)

#endif
