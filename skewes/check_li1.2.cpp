//
// File: check_li1.2.cpp
//
// Check that li(x)>pi(x)
//
// Input: start and end x
//        start and end pi(x)
//

#include "../includes/int_double12.0.h"
#include "../primesieve/src/soe/PrimeSieve.h"

#define MAX_N (1000000000000000000L) // 1e18
#define LOG_2_K (20)
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

// Global holding max possible truncation error from li
int_double li_err;

// setup li_terms, logns and li_err
// li_terms[i]=1/(i+1)/(i+1)!
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
// li(n)=gamma+sum_{i=1} log^{2i}(n)/n/n!
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
  // switch rounding mode for SSE
  // this makes int_doubles work
  _fpu_rndd();
  printf("Command line:-");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=6)
    {
      printf("Usage:- %s <x> <y> <pi(x)> <pi(y)> <cache size>.\n",argv[0]);
      exit(0);
    }
  long unsigned int x=atol(argv[1]);
  long unsigned int y=atol(argv[2]);
  long unsigned int pi_x=atol(argv[3]);
  long unsigned int pi_y=atol(argv[4]);

  if((x>=y)||(pi_x>=pi_y))
    {
      printf("Expect x<y and pi(x)<pi(y).\n");
      exit(0);
    }
  if(y>MAX_N)
    {
      printf("maximum y supported=%lu. Exiting.\n",MAX_N);
      exit(0);
    }
  if((x&1)||(y&1))
    {
      printf("specify x and y to be even to ensure they are composite.\n");
      exit(0);
    }

  int_double li_y=li(y),diff,log_y=log(int_double(y)),step,new_log_y;

  long unsigned int k=1<<LOG_2_K,ustep;
  int_double log_k=log(int_double(k-1)/k),li_inc,tmp;

  ps.setSieveSize(atol(argv[5])); // how many Kbyte cache

  for(;y>x;)
    {
      printf("Iterating with y=%lu, pi(y)=%lu, li(y)=",y,pi_y);
      print_int_double_str("",li_y);
      diff=li_y-pi_y;
      print_int_double_str(" difference in ",diff);
      if(diff.left<=0.0) // this may be because our estimate for li is now miles wide
	{
	  printf("Recomputing li(y)\n");
	  li_y=li(y);
	  diff=li_y-pi_y;
	  if(diff.left<=0.0) // it isn't that so might be a crossover
	    {
	      printf("Apparent crossover at %lu\npi(y)=%lu\n",y,pi_y);
	      print_int_double_str("li(y)=",li_y);
	      exit(0);
	    }
	}
      step=diff*(log_y+log_k); // if diff<y/k then 
      ustep=floor(step.left);
      ustep&=0xFFFFFFFFFFFFFFFEL; // round down to next even
      if((ustep<<LOG_2_K)>y)
	{
	  printf("Making too big a step for k. Scaling back.\n");
	  ustep=y>>LOG_2_K;
	  ustep&=0xFFFFFFFFFFFFFFFEL; // round down to next even
	}
      if(ustep==0) // no progress!
	{
	  printf("Step size is zero. Exiting.\n");
	  exit(0);
	}

      if(ustep>y-x)
	ustep=y-x;

      printf("ustep=%lu.\n",ustep);
      tmp=ustep/log_y;
      li_inc.left=tmp.left;
      pi_y-=ps.getPrimeCount(y-ustep,y);
      y-=ustep;
      log_y=log(int_double(y));
      tmp=ustep/log_y;
      li_inc.right=tmp.right;
      //print_int_double_str("Subtracting ",li_inc);
      li_y-=li_inc;
    }
  // now y=x, li_y is roughly li(y), pi_y should be pi(x)
  if(pi_x!=pi_y)
    {
      printf("Error in counting primes. Exiting.\n");
      exit(0);
    }
  diff=li_y-pi_y;
  if(diff.left<=0.0)
    {
      printf("Recomputing li(n)\n");
      li_y=li(y);
      diff=li_y-pi_y;
      if(diff.left<=0.0)
	{
	  printf("Apparent crossover at %lu\npi(y)=%lu\n",y,pi_y);
	  print_int_double_str("li(n)=",li_y);
	  exit(0);
	}
    }

  printf("Terminating with y=%lu, pi(y)=%lu, li(y)=",y,pi_y);
  print_int_double_str("",li_y);
  printf("Had to recompute li %d times.\n",li_call_count-1);
  return(0);
}
