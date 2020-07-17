#ifndef HURWITZ
#define HURWITZ
#include "int_double12.0.h"

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
  h_bernoulli[6]=h_bernoulli[6]/(182);                    //b(14)/14!
}

bool h_bernoulli_initialised=false;

#define DEFAULT_HUR_N (50)

// zeta(0,alpha) alpha in (0,1]
int_double hurwitz0(const int_double &alpha)
{
  if(!h_bernoulli_initialised)
    {
      set_h_bernoulli();
      h_bernoulli_initialised=true;
    }
  long unsigned int as[MAX_BERNOULLI_K]={0,6,100,49*72,761*288,21257200,2972885760};
  long unsigned int bs[MAX_BERNOULLI_K]={2,4,48,1440,288*280,7257600,958003200};
#define THIS_K (5) // this gives smallest abs error with q=99,991
  int_double res=0;
  long int N=15;//floor(exp(as[THIS_K-1]/(double) bs[THIS_K-1])+0.5); // this gives a local min of last term
  for(int n=0;n<N;n++)
    res+=sqr(log(alpha+n));
  //print_int_double_str("res=",res);
  int_double na=alpha+N;
  int_double lna=log(na);
  int_double lna2=sqr(lna);
  res+=lna2-na*lna2+2*na*lna-2*N-2*alpha-0.5*lna2;
  na=d_one/na;
  int_double na2=sqr(na);
  //
  // k=1..K
  // terms are B_2k/(2k)! *(a-b*log(N+alpha))/(N+alpha)^(2k-1)
  // b=2*(2k-2)!
  // a=2*(2k-2)!*sum(1/n,n=1..2k-2)
  int_double last_term;
  for(int k=0;k<THIS_K;k++)
    {
      last_term=h_bernoulli[k]*(as[k]-bs[k]*lna)*na;
      res+=last_term;
      na*=na2;
    }
  // the error term
  if(last_term.left<=last_term.right)
    last_term.right=last_term.left;
  else
    last_term.left=last_term.right;

  //print_int_double_str("Error term=",last_term);
  res+=last_term;
  //printf("zeta''(0,");print_int_double(alpha);print_int_double_str(")=",res);
  return(res);
}


//simple Euler Mac. calculation of zeta(s,alpha)
//cache s_mod,sigma s_array*bernoulli etc.
int_complex hurwitz (const int_complex &s, const int_double &alpha)
{
  if((s.real.left+s.real.right!=0.0)||(s.imag.left+s.imag.right!=0.0)||(s.real.left<=0.0)||(s.imag.left<0.0))
    {
      print_int_complex_str("s must be exact and in first quadrant",s);
      printf("Exiting.\n");
      exit(1);
    }
  unsigned int i,k=MAX_BERNOULLI_K,n=DEFAULT_HUR_N;
  int_double err,n_alpha,s_mod=sqrt(norm(s));
  int_complex s1=s,res=c_zero,n_alpha_s,term;
  int_complex s_array[MAX_BERNOULLI_K];

  if(!h_bernoulli_initialised)
    {
      set_h_bernoulli();
      h_bernoulli_initialised=true;
    }
//debug;
  if(n<-s_mod.right)
    n=(unsigned int) ceil(-s_mod.right);
  //print_int_double_str("alpha=",alpha);
  n_alpha=alpha+n;
  //debug;
  //print_int_double_str("n+alpha=",n_alpha);
  //print_int_complex_str("-s+1=",-s+1);
  n_alpha_s=pow(n_alpha,-s+1);
  //debug;
  for(i=0;i<n;i++)
    {
      res=res+pow(alpha+i,-s);
      //      if((i%1000)==0)
      //	printf("relative error in res is %d %d\n",rel_error(res.real),rel_error(res.imag));
    }

  //debug;
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

  err=sqrt(norm(term))/((s.real+2*k-1)*(s_mod+2*k-1));
  if(err.left<=err.right)
    err.right=err.left;
  else
    err.left=err.right;
  //print_int_double_str("error term is ",err);
  res.real=res.real+err;
  res.imag=res.imag+err;
  return(res);
}

inline int_complex real_hur(double s, int_double alpha)
{
	return(hurwitz(int_complex(s,0),alpha));
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

int_double re_ln_gamma_err=int_double(0);
inline int_double lngamma(const int_double x)
{
  int_double res,z_2n_1,last_term;
  int n;
  if(!h_bernoulli_initialised)
    {
      //printf("Initialising Bernoulli Constants.\n");
      set_h_bernoulli();
      //printf("Initialised Bernoulli Constants.\n");
      h_bernoulli_initialised=true;
    }
  if(x.left<=0)
    {
      print_int_double_str("lngamma for real arguments only works for positive arguments. Exiting. ",x);
      exit(0);
    }
  if(x.left>=10.0) // big enough
    {
      //print_int_double_str("Doing lngamma on ",x);
      int_double x2=d_one/sqr(x);
      res=(x-0.5)*log(x)-x+d_ln_two_pi*0.5;
      res=res+bernoulli[0]/x/2;
      z_2n_1=d_one/x;
      n=1;
      while(n<MAX_BERNOULLI_K-1)
	{
	  z_2n_1=z_2n_1*x2;
	  last_term=bernoulli[n]*z_2n_1/((n+n+2)*(n+n+1));
	  res=res+last_term;
	  n++;
	}
      z_2n_1=z_2n_1*x2;
      last_term=bernoulli[n]*z_2n_1/((n+n+2)*(n+n+1));
      //print_int_double_str("lngamma called with ",x);
      //print_int_double_str("lngamma returning ",res);
      //
      // By AS, the error is < |first neglected term| and same sign
      
      if(last_term.left>=0)
	res.right+=last_term.right;
      else
	res.left+=last_term.left;
      //printf("lngamma(");print_int_double(x);printf(")=");print_int_double(res);printf("\n");
      return(res);
    }
  int_double new_x=x;
  res=0;
  while(new_x.left<10.0)
    {
      //print_int_double_str("new_x=",new_x);
      res-=log(new_x);
      new_x+=1;
    }
  return(res+lngamma(new_x));
}
     


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

inline int_double int_im1(double sigma,double t)
{
  int_double d_t=int_double(t);
  int_double d_s=int_double(sigma);
  int_double tan_x=atan2(d_t,d_s);
  int_double t_sqr=d_t*d_t;
  int_double s_sqr=d_s*d_s;
  int_double t_sqr_plus_s_sqr=t_sqr+s_sqr;
  int_double log_x=log(t_sqr_plus_s_sqr);
  int_double log_x1=log(t_sqr/s_sqr+1);
  return(tan_x*d_t*d_s-s_sqr*log_x1*0.5-tan_x*t*0.5+log_x1*sigma*0.25+log_x*t_sqr_plus_s_sqr*0.25-s_sqr*0.25-t_sqr*0.75);
}


// 0<sigma,t0<t1
inline int_double int_im_ln_gamma(double sigma, double t0, double t1)
{
  int_double err=int_double(t1-t0)/t0;
  err*=0.125;
  err.left=err.right;
  return(int_im1(sigma,t1)-int_im1(sigma,t0)+err);
}

inline int_double int_even(double t0, double t1)
{
  return(int_im_ln_gamma(0.25,t0/2.0,t1/2.0)*2.0);
}

inline int_double int_odd(double t0, double t1)
{
  return(int_im_ln_gamma(0.75,t0/2.0,t1/2.0)*2.0);
}


// cache values for efficiency when doing many Hurwitz computions using
// same s.
inline void hur_init(int_complex *s_array, double re_z, double im_z)
{
	int i;
	int_complex s=int_complex(int_double(re_z),int_double(im_z));

	s_array[0]=s;
	for(i=1;i<MAX_BERNOULLI_K;i++)
	{
		s_array[i]=s_array[i-1]*(s+1)*(s+2);
		s+=2;
	}

	for(i=0;i<MAX_BERNOULLI_K;i++)
		s_array[i]*=h_bernoulli[i];
}

// faster version of hurwitz using values cached by hur_init
inline int_complex hurwitz1(double re_z, double im_z, const int_double &alpha, int_complex *s_array)
{
	int N=(int) max(re_z+im_z,(double) DEFAULT_HUR_N),i;
	int_complex s=int_complex(int_double(re_z),int_double(im_z)),minus_s=-s;
	int_complex res=pow(alpha,minus_s),b_term;
	int_double err;

	for(i=1;i<N;i++)
		res+=pow(alpha+i,minus_s);

	res+=pow(alpha+N,minus_s+1)/(s-1)+pow(alpha+N,minus_s)/2;

	for(i=0;i<MAX_BERNOULLI_K;i++)
	{
		b_term=s_array[i]*pow(alpha+N,minus_s-i-i-1);
		res+=b_term;
	}

	err=mod(minus_s-2*MAX_BERNOULLI_K+1)*mod(b_term)/(re_z+2*MAX_BERNOULLI_K-1);
	err.left=err.right;
	res.real+=err;
	res.imag+=err;

	return(res);
}

// zeta(1/2,alpha)
inline int_double hurwitz_half(const int_double &alpha, int_complex *s_array)
{
  long int i;

  int_double res=d_one/sqrt(alpha),b_term;
  int_double err;

  for(i=1;i<DEFAULT_HUR_N;i++)
    res+=d_one/sqrt(alpha+i);

  b_term=sqrt(alpha+DEFAULT_HUR_N)*2.0;
  res-=b_term;
  res+=d_one/b_term;
  //res+=pow(alpha+N,minus_s+1)/(s-1)+pow(alpha+N,minus_s)/2;

  for(i=0;i<MAX_BERNOULLI_K;i++)
    {
      b_term=s_array[i].real*pow(alpha+DEFAULT_HUR_N,-0.5-i-i-1);
      res+=b_term;
    }

  err=(2*MAX_BERNOULLI_K+0.5)*abs(b_term)/(2*MAX_BERNOULLI_K-0.5);
  err.left=err.right;
  res+=err;
  return(res);
}

#endif
