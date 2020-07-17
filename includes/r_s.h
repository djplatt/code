#ifndef R_S
#define R_S

#define pi_8 ((long double) 0.39269908169872415480783042290993786052L)
#define pi ((long double) 3.141592653589793238462643383279502884L)
#define two_pi ((long double) 6.2831853071795864769252867665590057684L)

#define MAX_T ((double) 1e10)
#define MAX_N ((long unsigned int) 100000)

long double logs[MAX_N],sqrts[MAX_N];

void init_r_s()
{
  long unsigned int i;
  for(i=1;i<MAX_N;i++)
    {
      logs[i]=logl(i+1.0);
      sqrts[i]=1.0L/sqrtl(i+1.0);
    }
}

long double theta(long double t)
/*************************************************************************
*                                                                              *
* Approximation to theta(t)=Im{log[Pi(it/2-3/4)]}-t/2*log(pi)                  *
*                                                                              *
*************************************************************************/
{
  long double t2=t/2.0L;
  return(t2*logl(t2/pi) - t2 - pi_8
	 + 1.0L/(48.0L*t) + 7.0L/(5760.0L*t*t*t));
}

long double Z (long double t)
{
  long double tt=theta(t);
  long double s=sqrtl(t/two_pi);
  long unsigned int N=s;
  long double p=s-N;
  long double C0=cosl(two_pi*(p*p-p-0.0625L))/cosl(two_pi*p)*expl(0.25L*log(two_pi/t));
  long double res=cosl(tt);
  unsigned long int i;

  //printf("t=%20.18e\ntt=%20.18e\ns=%20.18e\nN=%lu\np=%20.18e\nC0=%20.18e\n",t,tt,s,N,p,C0);
  //main sum
  for(i=1;i<N;i++)
    res+=cosl(tt-t*logs[i])*sqrts[i];

  res*=2.0L;

  if(N&1)
    return(res+C0);
  else
    return(res-C0);
}


#endif
