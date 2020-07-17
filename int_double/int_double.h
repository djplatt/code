// a simple version of double interval arithmetic that avoids messing
// with rounding modes by doing nextbefore/nextafter after every operation
// may increase interval unnecessarily when results are exact
#ifndef INT_DOUBLE
#define INT_DOUBLE
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <algorithm>
#include "crlibm.h"

inline double min(const double a, const double b, const double c, const double d)
{
  return std::min(std::min(a,b),std::min(c,d));
}

inline double max(const double a, const double b, const double c, const double d)
{
  return std::max(std::max(a,b),std::max(c,d));
}

/*
double _nextbefore(const double);

inline double _nextafter (const double x)
{
  if(!isfinite(x)) return nan("");
  if(x==0.0) return DBL_MIN;
  double xx=x;
  uint64_t *i;
  i=(uint64_t *) &xx;
  i[0]+=(x>0.0 ? 1 : -1);
  return xx;
}

inline double _nextbefore (const double x)
{
  if(!isfinite(x)) return nan("");
  if(x==0.0) return -DBL_MIN;
  double xx=x;
  uint64_t *i;
  i=(uint64_t *) &xx;
  i[0]+=(x>0.0 ? -1 : 1);
  return xx;
}
*/


inline double _nextbefore(const double x)
{
  return nextafter(x,-INFINITY);
}

inline double _nextafter(const double x)
{
  return nextafter(x,INFINITY);
}
/*
#define _nextbefore(x) nextafter(x,-INFINITY)
#define _nextafter(x) nextafter(x,INFINITY)
*/

class int_double{
public:
  double left;
  double right;
	
inline  int_double ()                      // constructors
  {
  };
// does no rounding, see below
inline int_double(double l) 
{
  left=l;
  right=l;
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
    right=r;
}

   friend int_double operator + (const int_double &lhs, const int_double &rhs);
   friend int_double operator + (const int_double &lhs, const double rhs);
   friend int_double operator - (const int_double &lhs, const int_double &rhs);
   friend int_double operator - (const int_double &lhs, const double rhs);
   friend int_double operator - (const int_double &lhs);
   friend int_double operator * (const int_double &lhs, const int_double &rhs);
   friend int_double operator * (const int_double &lhs, const double rhs);
  friend int_double operator / (const int_double &lhs, const int_double &rhs);
  friend int_double operator / (const int_double &lhs, const double rhs);
  friend int_double operator += (int_double &lhs, const int_double &rhs);
  friend int_double operator += (int_double &lhs, const double rhs);
  friend int_double operator -= (int_double &lhs, const int_double &rhs);
  friend int_double operator -= (int_double &lhs, const double rhs);
  friend int_double operator *= (int_double &lhs, const int_double &rhs);
  friend int_double operator *= (int_double &lhs, const double rhs);
  friend int_double operator /= (int_double &lhs, const int_double &rhs);
  friend int_double operator /= (int_double &lhs, const double rhs);
  friend int operator >= (const int_double &lhs, const int_double &rhs);
  friend int operator <= (const int_double &lhs, const int_double &rhs);
  friend int operator > (const int_double &lhs, const int_double &rhs);
  friend int operator < (const int_double &lhs, const int_double &rhs);
}; // end class

void print_int_double(const int_double &x)
{
  printf("[%20.18e,%20.18e]",x.left,x.right);
};

void print_int_double_str(const char *s, const int_double &x)
{
  printf("%s [%20.18e,%20.18e]\n",s,x.left,x.right);
};

inline int_double operator + (const int_double &lhs, const int_double &rhs)
{
  return int_double(_nextbefore(lhs.left+rhs.left),_nextafter(lhs.right+rhs.right));
}

inline int_double operator + (const int_double &lhs, const double rhs)
{
  return int_double(_nextbefore(lhs.left+rhs),_nextafter(lhs.right+rhs));
}

inline int_double operator += (int_double &lhs, const int_double &rhs)
{
  lhs=lhs+rhs;
  return lhs;
}

inline int_double operator += (int_double &lhs, const double rhs)
{
  lhs=lhs+rhs;
  return lhs;
}

// unary minus
inline int_double operator - (const int_double &x)
{
  return int_double(-x.right,-x.left);
}

inline int_double operator -(const int_double &lhs, const int_double &rhs)
{
  return int_double(_nextbefore(lhs.left-rhs.right),_nextafter(lhs.right-rhs.left));
}

inline int_double operator -(const int_double &lhs, const double rhs)
{
  return int_double(_nextbefore(lhs.left-rhs),_nextafter(lhs.right-rhs));
}

inline int_double operator -=(int_double &lhs, const int_double &rhs)
{
  lhs=lhs-rhs;
  return lhs;
}

inline int_double operator -=(int_double &lhs, const double rhs)
{
  lhs=lhs-rhs;
  return lhs;
}
  
inline int_double operator * (const int_double &lhs, const int_double &rhs)
{
  double a=lhs.left*rhs.left;
  double b=lhs.left*rhs.right;
  double c=lhs.right*rhs.left;
  double d=lhs.right*rhs.right;
  double mx=max(a,b,c,d);
  double mn=min(a,b,c,d);
  return int_double(_nextbefore(mn),_nextafter(mx));
}  

inline int_double operator * (const int_double &lhs, const double rhs)
{
  if(rhs>=0.0)
    return int_double(_nextbefore(lhs.left*rhs),_nextafter(lhs.right*rhs));
  else
    return int_double(_nextbefore(lhs.right*rhs),_nextafter(lhs.left*rhs));
}


inline int_double operator *= (int_double &lhs, const int_double &rhs)
{
  lhs=lhs*rhs;
  return lhs;
}

inline int_double operator *= (int_double &lhs, const double rhs)
{
  lhs=lhs*rhs;
  return lhs;
}


bool contains_zero(const int_double &x)
{
  return (x.left<=0.0) && (x.right>=0.0);
}


inline int_double inv(const int_double &x)
{
  if(x.left<=0.0)
    {
      if(x.right>=0)
	return int_double(-nan(""),nan(""));
      else
	return int_double(_nextbefore(1.0/x.right),_nextafter(1.0/x.left));
    }
  return int_double(_nextbefore(1.0/x.left),_nextafter(1.0/x.right));
}


inline int_double operator / (const int_double &lhs, const int_double &rhs)
{ 
 if(contains_zero(rhs))
    return int_double(-nan(""),nan(""));
  double a=lhs.left/rhs.left;
  double b=lhs.left/rhs.right;
  double c=lhs.right/rhs.left;
  double d=lhs.right/rhs.right;
  double mx=max(a,b,c,d);
  double mn=min(a,b,c,d);
  return int_double(_nextbefore(mn),_nextafter(mx));
}

inline int_double operator / (const int_double &lhs, double rhs)
{
  if(rhs==0.0)
    return int_double(-nan(""),nan(""));
  if(rhs>0.0)
    return int_double(_nextbefore(lhs.left/rhs),_nextafter(lhs.right/rhs));
  return int_double(_nextbefore(lhs.right/rhs),_nextafter(lhs.left/rhs));
}

inline int_double operator /= (int_double &lhs, const int_double &rhs)
{
  lhs=lhs/rhs;
  return lhs;
}

inline int_double operator /= (int_double &lhs, const double rhs)
{
  lhs=lhs/rhs;
  return lhs;
}

inline int operator >= (const int_double &lhs, const int_double &rhs)
{
	return lhs.left>=rhs.right;
}

inline int operator > (const int_double &lhs, const int_double &rhs)
{
	return lhs.left>rhs.right;
}

inline int operator <= (const int_double &lhs, const int_double &rhs)
{
	return lhs.right<=rhs.left;
}

inline int operator < (const int_double &lhs, const int_double &rhs)
{
	return lhs.right<rhs.left;
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
//int_double d_ln_pi;
//
// log(2*pi)
//int_double d_ln_two_pi;

int_double d_one=int_double(1.0);



// use this if you know left!=right
inline int_double log_two_sided(const int_double &x)
{
  return int_double(log_rd(x.left),log_ru(x.right));
}

inline int_double log(const int_double &x)
{
  if(x.left==x.right) // faster?
    {
      double rd=log_rd(x.left);
      return int_double(rd,_nextafter(rd));
    }
  return log_two_sided(x);
}

inline int_double exp (const int_double &x)  // nicely increasing
{
  return int_double(exp_rd(x.left),exp_ru(x.right));
}

// log(1+x)
inline int_double log1p (const int_double &x) // nicely increasing
{
  return int_double(log1p_rd(x.left),log1p_ru(x.right));
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
  return mask;
}

inline int make_mask(const int_double &x, const int_double &y)
{
  return make_mask1(x.left, x.right, y.left, y.right);
}

inline void sin_cos1(const double cos_left, const double cos_right,
		     const double sin_left, const double sin_right,
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
      sin_x->right=_nextafter(sin_left);
      cos_x->left=cos_left;
      cos_x->right=_nextafter(cos_right);
      return;
    case 2:
      sin_x->left=_nextbefore(sin_right);  // sin_right < 0 and too close to 0
      sin_x->right=_nextafter(sin_left);
      cos_x->left=-1.0;
      if(cos_left>=cos_right)
	cos_x->right=_nextafter(cos_left);
      else
	cos_x->right=_nextafter(cos_right);
      return;
    case 3:
      sin_x->left=sin_right;
      sin_x->right=_nextafter(sin_left);
      cos_x->left=cos_right;
      cos_x->right=_nextafter(cos_left);
      return;
    case 4:
      sin_x->left=-1.0;
      if(sin_right>=sin_left)
	sin_x->right=_nextafter(sin_right);
      else
	sin_x->right=_nextafter(sin_left);
      cos_x->left=_nextbefore(cos_left);
      cos_x->right=_nextafter(cos_right);
      return;
    case 11:
      if(sin_left<=sin_right)
	sin_x->left=sin_left;
      else
	sin_x->left=sin_right;
      sin_x->right=1.0;
      cos_x->left=cos_right;
      cos_x->right=_nextafter(cos_left);
      return;
    case 12:
      sin_x->left=sin_left;
      sin_x->right=_nextafter(sin_right);
      cos_x->left=cos_left;
      cos_x->right=_nextafter(cos_right);
      return;
    case 13:
      cos_x->right=1.0;
      if(cos_left<=cos_right)
	cos_x->left=cos_left;
      else
	cos_x->left=cos_right;
      sin_x->left=sin_left;
      sin_x->right=_nextafter(sin_right);
      return;
    case 15:
      sin_x->left=sin_left;
      sin_x->right=_nextafter(sin_right);
      cos_x->left=cos_right;
      cos_x->right=_nextafter(cos_left);
      return;
    default:
      printf("Weird error in sin_cos, mask = %d Exiting.\n",(int) mask);
      printf("\ncos_left: %20.18e\ncos_right: %20.18e\n",cos_left,cos_right);
      printf("sin_left: %20.18e\nsin_right: %20.18e\n",sin_left,sin_right);
      exit(1);
    }
}


//
// sin (pi*x), cos(pi*x)
//
inline void sin_cospi(const int_double &x, int_double *sin_x, int_double *cos_x)
{
  double cos_left,cos_right,sin_left,sin_right;

  if(x.right-x.left>=0.5)
    {
      sin_x->left=-1;sin_x->right=1;
      cos_x->left=-1;cos_x->right=1;
      return;

    }

  cos_left=cospi_rd(x.left);
  cos_right=cospi_rd(x.right); // probably a nop because cos is even
  sin_left=sinpi_rd(x.left);
  sin_right=sinpi_rd(x.right);
  sin_cos1(cos_left,cos_right,sin_left,sin_right,sin_x,cos_x);
}

// perform sin(x) into sin_x and cos(x) into cos_x
inline void sin_cos(const int_double &x, int_double *sin_x, int_double *cos_x)
{
  sin_cospi(x/d_pi,sin_x,cos_x);
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
  return(int_double(atan_rd(x.left),atan_ru(x.right)));
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

inline int_double sqrt(const int_double &x)
{
  return int_double(_nextbefore(sqrt(x.left)),_nextafter(sqrt(x.right)));
}

inline int_double sqr(const int_double &x)
{
  if((x.left>=0.0)||(x.right<=0.0))
    return x*x;
  return int_double(0.0,_nextafter(std::max(-x.left,x.right)));
}

#endif


