#include "../includes/int_double12.0.h"
#include "../primesieve/src/soe/PrimeSieve.h"

#define MAX_N (1000000000000000000LL) // 1e18
#define LI_TERMS (127)
#define LOG_TERMS (7) // 2^LOG_TERMS > LI_TERMS

int_double li_terms[LI_TERMS];
int_double logns[LOG_TERMS];

bool li_init=true;


int_double comp_bin(long unsigned int n)
{
  int_double res=1.0;
  //printf("comp_bin called with %lu\n",n);
  for(int i=0;i<LOG_TERMS;i++)
    {
      if(n&1)
	res*=logns[i];
      n>>=1;
    }
  //print_int_double_str("comp_bin returning ",res);
  return(res);
}

// 
int_double li_err;//=int_double(0.0,1.0e-21); // error is >=0

void init_li()
{
  li_terms[0]=int_double(1.0);
  for(int i=1;i<LI_TERMS;i++)
    li_terms[i]=li_terms[i-1]/(i+1);
  for(int i=1;i<LI_TERMS;i++)
    li_terms[i]/=(i+1);
  logns[0]=log(int_double(MAX_N));
  for(int i=1;i<LOG_TERMS;i++)
    logns[i]=sqr(logns[i-1]);
  li_err=comp_bin(LI_TERMS)*li_terms[LI_TERMS-1]*(LI_TERMS-1)/(LI_TERMS*LI_TERMS)/(int_double(1.0)-logns[0]/LI_TERMS);
  li_err.left=0.0;
  print_int_double_str("li error set to ",li_err);
}


int li_call_count=0;


//
// 5.1.10 of Abramowitz and Stegun
//
int_double li(unsigned long int n)
{
  li_call_count++;
  if(li_init)
    {
      init_li();
      li_init=false;
      //for(int i =0;i<10;i++) print_int_double_str("li_terms=",li_terms[i]);
    }
  //printf("Li called with %lu\n",n);
  int_double nn=n,res=int_double(0.5772,0.5773); // gamma
  logns[0]=log(nn);
  res+=log(logns[0]);
  for(int i=1;i<LOG_TERMS;i++)
    logns[i]=sqr(logns[i-1]);

  for(int i=0;i<LI_TERMS;i++)
    res+=comp_bin(i+1)*li_terms[i];

  return(res+li_err);
}

PrimeSieve ps;
long unsigned int pi_n;

int main(int argc,char **argv)
{
  _fpu_rndd();

  long unsigned int n=atol(argv[1]);
  double ustep;
  long unsigned int end=atol(argv[2]);
  if(end>MAX_N)
    {
      printf("maximum N supported=%lu. Exiting.\n",MAX_N);
      exit(0);
    }
  pi_n=atol(argv[3]);
  int_double li_n=li(n),diff,step,logn=log(int_double(n)),li_inc,tmp;

  ps.setSieveSize(64); // use 64Kbyte cache

  for(;n<end;)
    {
      //printf("Iterating with n=%lu, pi(n)=%lu, li(n)=",n,pi_n);
      //print_int_double(li_n);
      diff=li_n-pi_n;
      //print_int_double_str(" difference in ",diff);
      if(diff.left<=0.0) // this may be because our estimate for li is now miles wide
	{
	  printf("Recomputing li(n)\n");
	  li_n=li(n);
	  diff=li_n-pi_n;
	  if(diff.left<=0.0) // it isn't that so might be a crossover
	    {
	      printf("Apparent crossover at %lu\npi(n)=%lu\n",n,pi_n);
	      print_int_double_str("li(n)=",li_n);
	      exit(0);
	    }
	}
      step=diff*logn; // see Saouter and DeMichel 2010
      ustep=floor(step.left);
      if(ustep==0.0) ustep=1.0;
      tmp=ustep/logn;
      li_inc.right=tmp.right;
      pi_n+=ps.getPrimeCount(n+1,n+ustep);
      n+=ustep;
      logn=log(int_double(n));
      tmp=ustep/logn;
      li_inc.left=tmp.left;
      //print_int_double_str("Adding ",li_inc);
      li_n+=li_inc;
    }
  diff=li_n-pi_n;
  if(diff.left<=0.0)
    {
      printf("Recomputing li(n)\n");
      li_n=li(n);
      diff=li_n-pi_n;
      if(diff.left<=0.0)
	{
	  printf("Apparent crossover at %lu\npi(n)=%lu\n",n,pi_n);
	  print_int_double_str("li(n)=",li_n);
	  exit(0);
	}
    }

  if((n>MAX_N)&&(li_call_count>1)) // can't trust li_err anymore
    {
      printf("n exceeded MAX_N (%lu). Can't trust error term recomputing Li(n). Exiting.\n",MAX_N);
      exit(0);
    }

  printf("Terminating with n=%lu, pi(n)=%lu, li(n)=",n,pi_n);
  print_int_double_str("",li_n);
  printf("Had to recompute li %d times.\n",li_call_count-1);
  return(0);
}
