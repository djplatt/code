#ifndef INT_COMPLEX
#define INT_COMPLEX

#include "int_double.h"

#define print_int_complex_str(str,x) {printf(str);printf(" ");print_int_complex(x);printf("\n");}
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

  inline int_complex ()
    {};
  
  /*
  // This is a step too far 
  inline int_complex(const double re)
  {
  real=int_double(re);
  imag=int_double(0.0);
  }
  */

  inline int_complex(const int_double &re)
    {
      real=re;
      imag=int_double(0.0);
    }

  inline int_complex (const int_double &re,const int_double &im)
    {
      real=re;
      imag=im;
    };


  friend int_complex operator + (const int_complex &lhs, const int_complex &rhs);
  friend int_complex operator + (const int_complex &lhs, const int_double &rhs);
  friend int_complex operator + (const int_double &lhs, const int_complex &rhs);
  friend int_complex operator + (const int_complex &lhs, const double rhs);
  friend int_complex operator + (const double lhs, const int_complex &rhs);

  friend int_complex operator * (const int_complex &lhs, const int_complex &rhs);
  friend int_complex operator * (const int_complex &lhs, const int_double &rhs);
  friend int_complex operator * (const int_double &lhs, const int_complex &rhs);
  friend int_complex operator * (const int_complex &lhs, const double rhs);
  friend int_complex operator * (const double lhs, const int_complex &rhs);

  friend int_complex operator - (const int_complex &lhs, const int_complex &rhs);
  friend int_complex operator - (const int_complex &lhs, const int_double &rhs);
  friend int_complex operator - (const int_complex &lhs, const double rhs);
  friend int_complex operator - (const int_complex &lhs);

  friend int_complex operator / (const int_complex &lhs, const int_complex &rhs);
  friend int_complex operator / (const int_complex &lhs, const int_double &rhs);
  friend int_complex operator / (const int_complex &lhs, const double rhs);

  friend int_complex operator += (int_complex &lhs, const int_complex &rhs);
  friend int_complex operator += (int_complex &lhs, const int_double &rhs);
  friend int_complex operator += (int_complex &lhs, const double rhs);
  friend int_complex operator -= (int_complex &lhs, const int_complex &rhs);
  friend int_complex operator -= (int_complex &lhs, const int_double &rhs);
  friend int_complex operator -= (int_complex &lhs, const double rhs);
  friend int_complex operator *= (int_complex &lhs, const int_complex &rhs);
  friend int_complex operator *= (int_complex &lhs, const int_double &rhs);
  friend int_complex operator *= (int_complex &lhs, const double rhs);
  friend int_complex operator /= (int_complex &lhs, const int_complex &rhs);
  friend int_complex operator /= (int_complex &lhs, const int_double &rhs);
  friend int_complex operator /= (int_complex &lhs, const double rhs);
};

int_complex exp(const int_complex &);

int_complex e(const int_double &);

void print_int_complex(const int_complex &);

int_complex pow (const int_double &,const int_complex &);

int_complex pow1(const double, const double, const double);

int_complex pow1(const int_double &, const double, const double);

int_complex pow (const double, const int_complex &);

int_double norm(const int_complex &);

int_complex sqrt(const int_complex &);

int_complex conj (const int_complex &);

int_double argument (const int_complex &);

int_complex lngamma (const int_complex &);

int_complex hurwitz(const int_complex &, const int_double &);

int_double hurwitz0(const int_double &);

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

int_complex c_zero=int_complex(int_double(0.0),int_double(0.0));
int_complex c_one=int_complex(d_one,int_double(0.0));
int_complex c_half=int_complex(int_double(0.5),int_double(0.0));

inline int contains_zero(const int_complex &z)
{
	return(contains_zero(z.real)&&contains_zero(z.imag));
}

inline int_complex conj(const int_complex &rhs)
{
  return int_complex(rhs.real,-rhs.imag);
}

inline int_complex operator + (const int_complex &lhs, const int_complex &rhs)
{
  return int_complex(lhs.real+rhs.real,lhs.imag+rhs.imag);
}

inline int_complex operator + (const int_complex &lhs, const int_double &rhs)
{
  return int_complex(lhs.real+rhs,lhs.imag);
}

inline int_complex operator + (const int_double &rhs, const int_complex &lhs)
{
  return int_complex(lhs.real+rhs,lhs.imag);
}

inline int_complex operator + (const int_complex &lhs, const double rhs)
{
  return int_complex(lhs.real+rhs,lhs.imag);

}

inline int_complex operator + (const double rhs, const int_complex &lhs)
{
  return int_complex(lhs.real+rhs,lhs.imag);
}

inline int_complex operator += (int_complex &lhs, const int_complex &rhs)
{
	lhs.real+=rhs.real;
	lhs.imag+=rhs.imag;
	return(lhs);
}

inline int_complex operator += (int_complex &lhs, const int_double &rhs)
{
	lhs.real+=rhs;
	return(lhs);
}

inline int_complex operator += (int_complex &lhs, const double rhs)
{
	lhs.real+=rhs;
	return(lhs);
}

inline int_complex operator - (const int_complex &lhs)
{
  return int_complex(-lhs.real,-lhs.imag);
}

inline int_complex operator - (const int_complex &lhs, const int_complex &rhs)
{
  return int_complex(lhs.real-rhs.real,lhs.imag-rhs.imag);
}

inline int_complex operator - (const int_complex &lhs, const int_double &rhs)
{
  return int_complex(lhs.real-rhs,lhs.imag);
}

inline int_complex operator - (const int_complex &lhs, const double rhs)
{
  return int_complex(lhs.real-rhs,lhs.imag);
}

inline int_complex operator -= (int_complex &lhs, const int_complex &rhs)
{
	lhs.real-=rhs.real;
	lhs.imag-=rhs.imag;
	return(lhs);
}

inline int_complex operator -= (int_complex &lhs, const int_double &rhs)
{
	lhs.real-=rhs;
	return(lhs);
}

inline int_complex operator -= (int_complex &lhs, const double rhs)
{
  lhs.real-=rhs;
  return(lhs);
}

inline int_complex operator * (const int_complex &lhs, const int_complex &rhs)
{
  return int_complex(lhs.real*rhs.real-lhs.imag*rhs.imag,
		     lhs.real*rhs.imag+lhs.imag*rhs.real);
}

inline int_complex operator * (const int_complex &lhs, const int_double &rhs)
{
  return int_complex(lhs.real*rhs,lhs.imag*rhs);
}

inline int_complex operator * (const int_double &rhs, const int_complex &lhs)
{
  return int_complex(lhs.real*rhs,lhs.imag*rhs);
}

inline int_complex operator * (const int_complex &lhs, const double rhs)
{
  return int_complex(lhs.real*rhs,lhs.imag*rhs);
}

inline int_complex operator * (const double rhs, const int_complex &lhs)
{
  return int_complex(lhs.real*rhs,lhs.imag*rhs);
}
inline int_complex operator *= (int_complex &lhs,const int_complex &rhs)
{
  lhs=lhs*rhs;
  return lhs;
}

inline int_complex operator *= (int_complex &lhs, const int_double &rhs)
{
  lhs=lhs*rhs;
  return lhs;
}

inline int_complex operator *= (int_complex &lhs, const double rhs)
{
  lhs=lhs*rhs;
  return lhs;
}

inline int_complex operator / (const int_complex &lhs, const int_complex &rhs)
{
  return lhs*conj(rhs)/norm(rhs);
}

inline int_complex operator / (const int_complex &lhs, const int_double &rhs)
{
  return int_complex(lhs.real/rhs,lhs.imag/rhs);
}

inline int_complex operator / (const int_complex &lhs, const double rhs)
{
  return int_complex(lhs.real/rhs,lhs.imag/rhs);
}

inline int_complex operator /= (int_complex &lhs, const int_complex &rhs)
{
  lhs=lhs/rhs;
  return lhs;
}

inline int_complex operator /= (int_complex &lhs, const int_double &rhs)
{
  lhs=lhs/rhs;
  return lhs;
}

inline int_complex operator /= (int_complex &lhs, const double rhs)
{
  lhs=lhs/rhs;
  return lhs;
}

inline int_complex exp(const int_complex &z)
{
  int_double xs,xc;
  sin_cos(z.imag,&xs,&xc);
  return int_complex(xc,xs)*exp(z.real);
}

// exp(2*Pi*i*x)
inline int_complex e(const int_double &x)
{
  int_double xs,xc;
  sin_cospi(x*2,&xs,&xc);
  return(int_complex(xc,xs));
}

inline int_complex sqrt(const int_complex &z) // returns sqrt with arg in [-pi/2,pi/2]
{
  int_double theta,mod_z;
  int_complex res;
  if(contains_zero(z.imag)&&(z.real.right>0.0)) // don't try to atan this
    return(sqrt(-z)*int_complex(0,1));
  mod_z=pow(norm(z),0.25);
  theta=atan2(z.imag,z.real)/2.0;
  sin_cos(theta,&res.imag,&res.real);
  res=res*mod_z;
  return res;
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
  return z*exp(lnx*s.real);
}

inline int_complex pow (const double x, const int_complex &s)
{
  return pow(int_double(x),s);
}

inline int_complex pow1(const int_double &x, const double re, const double im)
{
  return pow(x,int_complex(re,im));
}

inline int_complex pow1 (const double x, const double re, const double im)
{
  return pow(int_double(x),int_complex(re,im));
}

inline int_double norm(const int_complex &z)
{
  return(sqr(z.real)+sqr(z.imag));
}

inline int_double mod(const int_complex &z)
{
  return(sqrt(norm(z)));
}

inline int_double argument(const int_complex &z)
{
  return(atan2(z.imag,z.real));
}

inline int_complex log(const int_complex &z)
{
  return(int_complex(log(norm(z))/2,argument(z)));
}

inline int_double read_int_double(FILE *infile)
{
  int y;
  double x[2];
  y=fread(x,sizeof(double),2,infile);
  return(int_double(x[0],x[1]));
}

inline int_complex read_int_complex(FILE *infile)
{
  int y;
  double x[4];
  y=fread(x,sizeof(double),4,infile);
  return(int_complex(int_double(x[0],x[1]),int_double(x[2],x[3])));
}


#endif
