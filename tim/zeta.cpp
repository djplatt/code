//
// zetalim1.0.cpp
//
// DJ Platt Jan 2014
//
// Find upb on |zeta(1/2+it)|
//
#include "inttypes.h"
#include "../includes/int_double12.0.h"

//
// simple version of Euler Mac for Zeta(s) with sum of length N
//
int_complex zeta(const int_complex &s, uint64_t N)
{
  int_complex res(1);
  for(uint64_t n=2;n<=N;n++)
    res+=pow(n,-s);
  res+=pow(N,-s+1)/(s-1)-pow(N,-s)/2;
  int_complex term=s*pow(N,-s-1)/12;
  res+=term;
  term=s*(s+1)*(s+2)*pow(N,-s-3)/720;
  res-=term;
  int_double factor=mod(s+1)/(s.real+1);
  int_double error=factor*mod(term);
  error.left=error.right;
  res.real+=error;
  res.imag+=error;
  return(res);
} 

int_double theta_log_pi_2,theta_err,theta_pi8;;
//
// theta(t)=1/2(Im(lngamma(1/4+it/2)-lngamma(1/4-it/2))-tlog(pi))
// the error is simply |theta(200)-mytheta(200)|
// 
int_double theta(int_double t)
{
  //int_double res=t*0.5*log(t/d_two_pi)-0.5*t-d_pi/8 +d_one/(t*48)+d_one*7/(5760*t*t*t);
  int_double res=t*0.5*(log(sqrt((sqr(t)*4+1))*0.25)-1-d_ln_pi)-0.125*(d_pi-0.5/t);
  return(res+theta_err);
}

int_double c01(double p)
{
  int_double pd=int_double(p);
  int_double s,c,c1;
  sin_cospi(2*(sqr(pd)-pd-1.0/16.0),&s,&c);
  sin_cospi(2*pd,&s,&c1);
  return(c/c1);
}

int_double c0(int_double rho)
{
  int_double l=c01(rho.left);
  int_double r=c01(-rho.right);
  if((rho.left<=0.5)&&(rho.right<=-0.5)) // rho straddled 1/2
    return(int_double(0.382,max(-l.right,-r.right))); // c0(1/2)>0.382
  else
    return(int_double(min(l.left,r.left),max(-l.right,-r.right)));
}

int_double rs_err;
//
// Riemann Siegel, a=sqrt(t/2Pi)
// sqrts[n]=1/sqrt(n)
// logs(n)=log(n)
//
// The first and only C_n term C_0=cos(2pi(p^2-p-1/16))/cos(2pi p) is replaced
// by [0.382,0.984]/sqrt(a)
// 
int_double rs(int_double t, int_double a, int_double *logs, int_double *sqrts)
{
  uint64_t N=floor(a.left);
  int_double rho=a-N;
  int_double res=0,si,co;
  int_double tt=theta(t);
  //print_int_double_str("theta(t)=",tt);
  for(uint64_t n=1;n<=N;n++)
    {
      sin_cos(t*logs[n]-tt,&si,&co); // sin is thrown away
      res+=co*sqrts[n];
    }
  res+=res;
  if(N&1)
    res+=c0(rho)/sqrt(a);
  else
    res-=c0(rho)/sqrt(a);
  int_double tp=pow(t,-1.25); // Gabcke 0.053*t^(-5/4) for t>=200
  int_double err=rs_err*tp;
  err.left=err.right;
  res+=err;
  res=sqrt(sqr(res));
  return(res);
}

double THRESH;
bool too_wide (int_double x)
{
  return((x.left+x.right)/x.right>THRESH);
}

int main(int argc, char **argv)
{
  printf("Command Line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=2)
    {
      printf("Usage:- %s <t>\n",argv[0]);
      exit(0);
    }

  // set up the SSE system for rounding and initialise constants
  _fpu_rndd();

  int_double t(atof(argv[1]));
  int_double a=sqrt(t/(d_pi*2));
  uint64_t MAX_N=floor(sqrt(t.left/6.283));
  int_double *logs=(int_double *)malloc(sizeof(int_double)*(MAX_N+1));
  int_double *sqrts=(int_double *)malloc(sizeof(int_double)*(MAX_N+1));
  for(uint64_t n=1;n<=MAX_N;n++)
    {
      logs[n]=log(int_double(n));
      sqrts[n]=d_one/sqrt(int_double(n));
    }

  rs_err=53;
  rs_err/=1000.0;
  theta_log_pi_2=log(d_pi_2); // log(Pi/2)
  theta_pi8=d_pi*0.125; // Pi/8
  theta_err=int_double(-21,21);
  theta_err/=100000.0;
  print_int_double_str("rs returned  ",rs(t,a,logs,sqrts));
  print_int_double_str("E-M returned ",mod(zeta(int_complex(int_double(0.5),t),floor(t.left*10))));
  return(0);
}
