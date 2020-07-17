#include "inttypes.h"
#include "../includes/int_double13.0.h"

// simple version of 
int_complex zeta(const int_complex &s, uint64_t N)
{
  //print_int_complex_str("In E-M with s=",s);
  //printf("and N=%lu\n",N);
  int_complex res(1);
  for(uint64_t n=2;n<=N;n++)
    res+=pow(n,-s);
  //print_int_complex_str("sum of powers=",res);
  res+=pow(N,-s+1)/(s-1)-pow(N,-s)/2;
  //print_int_complex_str("result before Bernoulli terms=",res);
  int_complex term=s*pow(N,-s-1)/12;
  res+=term;
  term=s*(s+1)*(s+2)*pow(N,-s-3)/720;
  res-=term;
  int_double factor=mod(s+1)/(s.real+1);
  int_double error=factor*mod(term);
  error.left=error.right;
  //print_int_double_str("error=",error);
  res.real+=error;
  res.imag+=error;
  //print_int_complex_str("Returning ",res);
  return(res);
} 
  
int_double theta(int_double t)
{
  int_double res=t*0.5*log(t/d_two_pi)-0.5*t-d_pi/8 +d_one/(t*48)+d_one*7/(5760*t*t*t);
  return(res);
}

//Riemann Siegel
int_double rs(int_double t, int_double *logs, int_double *sqrts)
{
  //print_int_double_str("in rs with t=",t);
  int_double(a)=sqrt(t/d_two_pi);
  uint64_t N=floor(a.left);
  if(floor(-a.right)!=N)
    {
      printf("resorting to E-M at N=%lu\n",N+1);
      int_complex s(int_double(0.5),t);
      return(mod(zeta(s,floor(-t.right))));
    }

  int_double thetat=theta(t);
  //print_int_double_str("theta(t)=",thetat);
  int_double res=0;
  for(uint64_t n=1;n<=N;n++)
    {
      int_double si,co;
      sin_cos(t*logs[n]-thetat,&si,&co);
      res+=co*sqrts[n];
    }
  //print_int_double_str("rs returning ",res);
  return(2*res);
}

int main(int argc, char **argv)
{
  _fpu_rndd();

  double full_delta=1.0/16.0;
  double ts=atof(argv[1]);
  double te=atof(argv[2]);
  uint64_t MAX_N=floor(sqrt(atof(argv[2])/6.283));
  int_double *logs=(int_double *)malloc(sizeof(int_double)*(MAX_N+1));
  int_double *sqrts=(int_double *)malloc(sizeof(int_double)*(MAX_N+1));
  for(uint64_t n=1;n<=MAX_N;n++)
    {
      logs[n]=log(int_double(n));
      sqrts[n]=d_one/sqrt(int_double(n));
    }

  while(ts<te)
    {
      double this_delta=full_delta;
      int_double t(ts);
      int_double delta(0,full_delta);
      t+=delta;
      int_double res=rs(t,logs,sqrts)-log(t)*pow(t,int_double(1.0)/6)*int_double(7)/10;
      uint64_t count=0;
      while(res.right<0) // too large
	{
	  printf("reducing delta.\n");
	  this_delta/=2.0;
	  t=int_double(ts);
	  delta=int_double(0,this_delta);
	  t+=delta;
	  print_int_double_str("calling rs with t=",t);
	  res=rs(t,logs,sqrts);
	  print_int_double_str("rs returned=",res);
	  res=res-log(t)*pow(t,int_double(1.0)/6)*int_double(7)/10;
	  if(++count>100)
	    exit(0);
	}
      ts=-t.right;
    }
  return(0);
}
