#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>


using namespace std;


#define _fpu_rndd _control87(_RC_DOWN,_MCW_RC)

#define _fpu_rndu _control87(_RC_UP,_MCW_RC)

inline double nextafter(double);

void fatal_error(const char *);

class int_double{
public:
  double left;
  double right;
  inline int_double ();                      // constructors
  inline int_double (/* const */ double);
  inline int_double (/* const */ double,/* const */ double);
  inline int_double operator + (/* const */ int_double &); // operators
  inline int_double operator + (/* const */ double);
  inline int_double operator += (/* const */ int_double &);
  inline int_double operator += (/* const */ double);
  inline int_double operator - (/* const */ int_double &); // int_double - int_double
  inline int_double operator - (/* const */ double);
  inline int_double operator - ()
  {
	  int_double temp;
      temp.left=right;
      temp.right=left;
      return(temp);
  };
  inline int_double operator -= (/* const */ int_double &);
  inline int_double operator -= (/* const */ double);
  inline int_double operator * (/* const */ int_double &);
  inline int_double operator * (/* const */ double);
  inline int_double operator *= (/* const */ int_double &);
  inline int_double operator *= (/* const */ double);
  inline int_double operator / (/* const */ int_double &);
  inline int_double operator / (/* const */ double);
//  inline int_double operator % (int_double &);
  inline int_double operator /= (/* const */ int_double &);
  inline int_double operator /= (/* const */ double);
  inline int operator >= (/* const */ int_double &);
  inline int operator > (/* const */ int_double &);
};

extern inline int_double exp (/* const */ int_double &);

extern inline int_double log (/* const */ int_double &);

extern inline void sin_cos(/* const */ int_double, int_double *, int_double *);

extern inline int_double sqr(/* const */ int_double &);

extern inline int_double pow(/* const */ int_double &, /* const */ int_double &);

extern inline int_double pow(/* const */ int_double &, /* const */ double);

extern inline int_double atan(/* const */ int_double &);

extern inline int_double sqrt(/* const */ int_double &);

extern inline int rel_error(int_double &);

extern inline int abs_error(int_double &);

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
	inline int_complex (/* const */ int_double &,/* const */ int_double &);
	inline int_complex operator + (/* const */ int_complex &); // operators
	inline int_complex operator + (/* const */ int_double &);
	inline int_complex operator + (/* const */ double);
	inline int_complex operator += (/* const */ int_complex &); // operators
	inline int_complex operator += (/* const */ int_double &);
	inline int_complex operator += (/* const */ double);
	inline int_complex operator - ();
	inline int_complex operator - (/* const */ int_complex &);
	inline int_complex operator - (/* const */ int_double &);
	inline int_complex operator - (/* const */ double);
	inline int_complex operator -= (/* const */ int_complex &);
	inline int_complex operator -= (/* const */ int_double &);
	inline int_complex operator -= (/* const */ double);
	inline int_complex operator * (/* const */ int_complex &);
	inline int_complex operator * (/* const */ int_double &);
	inline int_complex operator * (/* const */ double);
	inline int_complex operator *= (/* const */ int_complex &);
	inline int_complex operator *= (/* const */ int_double &);
	inline int_complex operator *= (/* const */ double);
	inline int_complex operator / (/* const */ int_complex &);
	inline int_complex operator / (/* const */ int_double &);
	inline int_complex operator / (/* const */ double);
	inline int_complex operator /= (/* const */ int_complex &);
	inline int_complex operator /= (/* const */ int_double &);
	inline int_complex operator /= (/* const */ double);
};

extern inline int_complex exp(/* const */ int_complex &);

extern inline void print_int_complex(/* const */ int_complex &);

extern inline int_complex pow (/* const */ int_double &, /* const */ int_complex &);

extern inline int_complex pow (/* const */ double, /* const */ int_complex &);

extern inline int_double norm(/* const */ int_complex &);

extern inline int_complex sqrt(/* const */ int_complex &);

extern inline void print_interval (/* const */ int_double &);

extern inline int_complex conj (/* const */ int_complex &);

#define real(z) z.real
#define imag(z) z.imag

static /* const */ int_double d_zero=int_double(0.0,0.0);
static /* const */ int_double d_one=int_double(1.0,1.0);
static /* const */ int_complex c_zero=int_complex(d_zero,d_zero);
static /* const */ int_complex c_one=int_complex(d_one,d_zero);

static /* const */ int_double d_half=int_double(0.5,0.5);

static int _i_pi[2]={1413754136,1074340347};   // this is what double pi looks like

static int _i_pi2[2]={1413754137,1074340347};  // this is double pi plus a bit

static double *_d_pi=(double *)&_i_pi;

static double *_d_pi2=(double *)&_i_pi2;

static int_double d_pi=int_double(_d_pi[0],_d_pi2[0]);

static int_double d_two_pi=d_pi*2;

unsigned long __nze[2]={0,0x80000000}; // -0.0

unsigned long __delta[2]={1,0};       // very small

double _nze = *(double*)__nze;

double _delta= *(double*)__delta;


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


void print_int_double(/* const */ int_double &x)
{
  printf("[%20.18e,%20.18e]",x.left,-x.right);
};

// don't preset to anything sensible
inline int_double::int_double()
{
}

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


// does no rounding, see below
inline int_double::int_double(double l) 
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
inline int_double::int_double(double l,double r)
{
    if(l>r)
    {
	printf("Error constructing int_double, right<left. Exiting.\n");
        exit(0);
    };
    left=l;
    right=-r;
}


inline int_double int_double::operator+ (/* const */ int_double &rhs)
{
	int_double temp;
	temp.left=left+rhs.left;
	temp.right=right+rhs.right;
	return(temp);
}

inline int_double int_double::operator+ (/* const */ double rhs)
{
	int_double temp;
	temp.left=left+rhs;
	temp.right=right-rhs;
	return(temp);
}

inline int_double int_double::operator+= (/* const */ int_double &rhs)
{
  left+=rhs.left;
  right+=rhs.right;
  return(*this);
}

inline int_double int_double::operator+= (/* const */ double rhs)
{
  left+=rhs;
  right-=rhs;
  return(*this);
}


inline int_double int_double::operator- (/* const */ int_double &rhs)
{
  int_double temp;
  temp.left=left+rhs.right;
  temp.right=right+rhs.left;
  return(temp);
}

inline int_double int_double::operator- (/* const */ double rhs)
{
  int_double temp;
  temp.left=left-rhs;
  temp.right=right+rhs;
  return(temp);
}

inline int_double int_double::operator-= (/* const */ int_double &rhs)
{
  left+=rhs.right;
  right+=rhs.left;
  return(*this);
}

inline int_double int_double::operator-= (/* const */ double rhs)
{
  left-=rhs;
  right+=rhs;
  return(*this);
}

// very nasty, probably no quicker, so not used
//inline void myneg(unsigned int *i)
//{
//    i[1]^=0x80000000;
//}

inline int_double int_double::operator* (/* const */ int_double &rhs)
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
inline int_double times_pos(int_double &lhs, int_double &rhs)
{
	int_double res;
	res.left=lhs.left*rhs.left;
	res.right=lhs.right*(-rhs.right);
	return(res);
}

#define times_pos_equals(lhs,rhs){lhs.left*=rhs.left;lhs.right*=-rhs.right;}



inline int_double int_double::operator*= (/* const */ int_double &rhs)
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

inline int_double int_double::operator* (/* const */ double rhs)
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

inline int_double int_double::operator*= (/* const */ double rhs)
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

inline int_double int_double::operator/ (/* const */ int_double &rhs)
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


inline int_double int_double::operator/= (/* const */ int_double &rhs)
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

inline int_double int_double::operator/ (/* const */ double rhs)
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

inline int_double int_double::operator/= (/* const */ double rhs)
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

void minmax4(double a, double b, double c, double d, double &min, double &max)
{
	if(a<b)
	{
	min=a;
	max=b;
	}
	else
	{
		min=b;
		max=a;
	}
if(c<min)
		min=c;
	else
	{
		if(c>max)
			max=c;
	}
if(d<min)
		min=d;
	else
	{
		if(d>max)
			max=d;
	}
}
/*
// if needed this will have to be fixed.
inline int_double int_double::operator% (int_double %rhs)
{
	double a1,a2,a3,a4,r1,r2;
	if(rhs.left<=0.0)
	{
		printf("Remainder by negative interval or one containing zero. Exiting. \n");
		exit(0);
	}
	r1=-right;
	r2=-rhs.right;
	a1=left%rhs.left;
	a2=left%r2;
	a3=r1%rhs.left;
	a4=r1%r2;
	minmax4(a1,a2,a3,a4,temp.left,temp.right);
	temp.right=-nextafter(temp.right);
	return(temp);
}
*/




inline int int_double::operator>= (/* const */ int_double &rhs)
{
	return(left>=(-rhs.right));
}

inline int int_double::operator> (/* const */ int_double &rhs)
{
	return(left>(-rhs.right));
}

// trap rh limit =0.0
inline int_double exp (/* const */ int_double &x)  // nicely increasing
{
	int_double temp;
	unsigned long int *i;
	temp.left=exp(x.left);
	if(x.right==0.0)
		temp.right=-1.0;
	else
	{
		temp.right=exp(-x.right);   // exp(x) >= 0 for all double x

		// very quick and dirty nextafter(+ve double)
		i=(unsigned long int *) &temp.right;
		i[0]++;
		if(i[0]==0)
			i[1]++;

		temp.right=-temp.right;
	}
	return(temp);
}

// if right hand endpoint is one, then interval is
// widened unnecessarily. could trap it but......
inline int_double log (/* const */ int_double &x) // nicely increasing
{
	int_double temp;

	temp.left=log(x.left);
	temp.right=log(-x.right);

	temp.right=nextafter(temp.right);

	temp.right=-temp.right;

	return(temp);
}

#define M_PI_2 (double) 1.57079632679489662


inline extern void sin_cos(/* const */ int_double x, int_double *sin_x, int_double *cos_x)
{
	double cos_left,cos_right,sin_left,sin_right;
	unsigned char mask=0;

	if(x.right+x.left<=-M_PI_2)
	{
		printf("Endpoints more than Pi/2 apart in sin_cos. Exiting.\n");

		exit(0);
		return;
	}
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
		sin_x->left=-nextafter(-sin_right);  // sin_right < 0 and too close to 0
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
		cos_x->left=-nextafter(-cos_left);
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

		/*
		printf("13: pathological sin_cos. exiting");
		exit(0);
		sin_x->right=-1.0;
		if(sin_left>=sin_right)
		sin_x->left=sin_right;
		else
		sin_x->left=sin_left;
		cos_x->left=cos_left;
		cos_x->right=-nextafter(cos_right);
		return;
		*/
	case 15:
		sin_x->left=sin_left;
		sin_x->right=-nextafter(sin_right);
		cos_x->left=cos_right;
		cos_x->right=-nextafter(cos_left);
		return;
	default:
		printf("Weird error in sin_cos, mask = %d exiting.\n",(int) mask);
		exit(0);
	}
}

inline extern int_double sqr(/* const */ int_double &x)
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

inline extern int_double pow(/* const */ int_double &x, /* const */ int_double &y)
{
	return(exp(log(x)*y));
}

inline extern int_double pow(/* const */ int_double &x, /* const */ double y)
{
	return(exp(log(x)*y));
}

/*
// not debugged yet
inline extern int_double atan(int_double &x)
{
	int_double temp;
	// check x doesn't span a discontinuity
	//
	if(floor(x.left/pi)!=floor(-x.right/pi))
	{
		printf("Error in atan, x spans a discontinuity. Exiting.\n");
		exit(0);
	}

	temp.left=atan(x.left);
	temp.right=atan(x.right);  // is this right?
	return(temp);
}
*/

extern inline int_double sqrt(/* const */ int_double &x) // strictly increasing
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

extern inline void print_interval (/* const */ int_double &x)
{
	print_int_double(x);
}


// we don't initialise by default
inline int_complex::int_complex()
{
}

inline int_complex::int_complex(/* const */ int_double &re, /* const */ int_double &im)
{
	real=re;
	imag=im;
} 

extern inline int_complex conj(/* const */ int_complex &rhs)
{
	int_double temp=d_zero;
	return(int_complex(rhs.real,temp-rhs.imag));
}

inline int_complex int_complex::operator+ (/* const */ int_complex &rhs)
{
	return(int_complex(real+rhs.real,imag+rhs.imag));
}

inline int_complex int_complex::operator+ (/* const */ int_double &rhs)
{
	return(int_complex(real+rhs,imag));
}

inline int_complex int_complex::operator+ (/* const */ double rhs)
{
	return(int_complex(real+rhs,imag));
}

inline int_complex int_complex::operator+= (/* const */ int_complex &rhs)
{
	real+=rhs.real;
	imag+=rhs.imag;
	return(*this);
}

inline int_complex int_complex::operator+= (/* const */ int_double &rhs)
{
	real+=rhs;
	return(*this);
}

inline int_complex int_complex::operator+= (/* const */ double rhs)
{
	real+=rhs;
	return(*this);
}

inline int_complex int_complex::operator- ()
{
	return(int_complex(-real,-imag));
}

inline int_complex int_complex::operator- (/* const */ int_complex &rhs)
{
	return(int_complex(real-rhs.real,imag-rhs.imag));
}

inline int_complex int_complex::operator- (/* const */ int_double &rhs)
{
	return(int_complex(real-rhs,imag));
}

inline int_complex int_complex::operator- (/* const */ double rhs)
{
	return(int_complex(real-rhs,imag));
}

inline int_complex int_complex::operator-= (/* const */ int_complex &rhs)
{
	real-=rhs.real;
	imag-=rhs.imag;
	return(*this);
}

inline int_complex int_complex::operator-= (/* const */ int_double &rhs)
{
	real-=rhs;
	return(*this);
}

inline int_complex int_complex::operator-= (/* const */ double rhs)
{
	real-=rhs;
	return(*this);
}

inline int_complex int_complex::operator* (/* const */ int_complex &rhs)
{
	return(int_complex(real*rhs.real-imag*rhs.imag,
		real*rhs.imag+imag*rhs.real));
}

inline int_complex int_complex::operator* (/* const */ int_double &rhs)
{
	return(int_complex(real*rhs,imag*rhs));
}

inline int_complex int_complex::operator* (/* const */ double rhs)
{
	return(int_complex(real*rhs,imag*rhs));
}

inline int_complex int_complex::operator*= (/* const */ int_complex &rhs)
{
	int_double old_real;
	old_real=real;
	real*=rhs.real;
	real-=imag*rhs.imag;
	imag=old_real*rhs.imag+imag*rhs.real;
	return(*this);
}

inline int_complex int_complex::operator*= (/* const */ int_double &rhs)
{
	real*=rhs;
	imag*=rhs;
	return(*this);
}

inline int_complex int_complex::operator*= (/* const */ double rhs)
{
	real*=rhs;
	imag*=rhs;
	return(*this);
}

inline int_complex int_complex::operator/ (/* const */ int_complex &rhs)
{
	return(*this*conj(rhs)/(sqr(rhs.real)+sqr(rhs.imag)));
}

inline int_complex int_complex::operator/ (/* const */ int_double &rhs)
{
	return(int_complex(real/rhs,imag/rhs));
}

inline int_complex int_complex::operator/ (/* const */ double rhs)
{
	return(int_complex(real/rhs,imag/rhs));
}

inline int_complex int_complex::operator/= (/* const */ int_complex &rhs)
{
	int_double den=sqr(rhs.real)+sqr(rhs.imag);

	*this*=conj(*this);
	*this/=den;
	return(*this);
}

inline int_complex int_complex::operator/= (/* const */ int_double &rhs)
{
	real/=rhs;
	imag/=rhs;
	return(*this);
}

inline int_complex int_complex::operator/= (/* const */ double rhs)
{
	real/=rhs;
	imag/=rhs;
	return(*this);
}

extern inline int_complex exp(/* const */ int_complex &z)
{
	int_double xs,xc;

	sin_cos(z.imag,&xs,&xc);
	return(int_complex(xc,xs)*exp(z.real));
}

extern inline void print_int_complex(/* const */ int_complex &z)
{
	print_int_double(z.real);
	printf("+");
	print_int_double(z.imag);
	printf("i");
}

extern inline int_complex pow (/* const */ int_double &x, /* const */ int_complex &s)
{
	int_double tlnx,x_to_sigma,ts,tc;
//	int rts,rtc;
//	unsigned int n_2_pi;

	tlnx=log(x);
	x_to_sigma=exp(tlnx*s.real);
	tlnx*=s.imag;
	sin_cos(tlnx,&ts,&tc);
/*
	rts=rel_error(ts);
	rtc=rel_error(tc);
	if(((rts!=0)&&(rts>=-8))||((rtc!=0)&&(rtc>=-8)))
		printf(" ");
*/
	return(int_complex(x_to_sigma*tc,x_to_sigma*ts));
}

extern inline int_complex pow (/* const */ double x, /* const */ int_complex &s)
{
return(pow(int_double(x),s));
}

extern inline int_double norm(/* const */ int_complex &z)
{
	return(sqr(z.real)+sqr(z.imag));
}

/*
extern inline int_complex sqrt(int_complex &z)
{
	int_double r,theta;
	int_complex z1;

	r=pow(norm(z),0.25);

	theta=atan(z.imag/z.real);  // will fail if 0 in Re(z)

	if(z.real.left<0.0) // in quadrant 2 or 3
	{
		if(z.imag.left<=0.0)
		{
			if(z.imag.right>=0.0) // quadrant 3
				theta=(theta-d_pi)/2.0;
			else
				fatal_error("Imaginary part of theta contains zero in atan. Exiting");
		}
		else
		{
			if(z.imag.right<=0.0) // quadrant 2
				theta=(theta+d_pi)/2.0;
			else
				fatal_error("Imaginary part of theta contains zero in atan. Exiting");
		}
	}
	else
		theta=theta/2.0; // quadranr 1 or 4

	sin_cos(theta,&z1.imag,&z1.real);
	return(z1*r);

}
*/

#define nextafter(x,y) nextafter(x)

