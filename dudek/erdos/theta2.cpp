#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "math.h"
#include "primesieve.h"
#include "../includes/int_double12.0.h"


#define MAX_X (10000000000)
#define NQ (24)
#define MAX_Q (97)

uint64_t qs[NQ]={3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
uint64_t phi_qs[NQ];
uint64_t q2s[NQ];
double worsts[NQ];
uint64_t as[NQ];
uint64_t xs[NQ];

int_double thetas[NQ][MAX_Q*MAX_Q];

void do_prime(uint64_t p)
{
  int_double pd=p;
  int_double sp=sqrt(pd);
  int_double lp=log(pd);
  for(uint64_t qq=0;qq<NQ;qq++)
    {
      uint64_t a=p%q2s[qq];
      if(a%qs[qq]==0)
	continue;
      int_double th=thetas[qq][a];
      th-=pd/phi_qs[qq];
      if(th.left<0) th=-th;
      th/=sp;
      if(-th.right>worsts[qq])
	{
	  //printf("New record for a=%lu q=%lu just before p=%lu with theta=%20.18e\n",a,qs[qq],p,-th.right);
	  worsts[qq]=-th.right;
	  as[qq]=a;xs[qq]=p;
	}
      thetas[qq][a]+=lp;
      th=thetas[qq][a];
      //printf("after p=%lu thetas[%lu] now %20.18e\n",p,a,-th.right); 
      th-=pd/phi_qs[qq];
      if(th.left<0) th=-th;
      th/=sp;
      if(-th.right>worsts[qq])
	{
	  //printf("New record for a=%lu q=%lu just after  p=%lu with theta=%20.18e\n",a,qs[qq],p,-th.right);
	  worsts[qq]=-th.right;
	  as[qq]=a;xs[qq]=p;
	}
    }
}

int main(int argc, char **argv)
{

  _fpu_rndd();
  for(uint64_t qq=0;qq<NQ;qq++)
    {
      uint64_t q=qs[qq];
      q2s[qq]=q*q;
      phi_qs[qq]=q*(q-1);
      worsts[qq]=0.0;
      for(uint64_t i=1;i<q2s[qq];i++)
	thetas[qq][i]=0.0;
    }
  primesieve_set_sieve_size(atol(argv[1]));
  primesieve_callback_primes(2,(uint64_t) MAX_X,do_prime);

  for(uint64_t qq=0;qq<NQ;qq++)
    {
      uint64_t q=qs[qq];
      printf("worst seen for q=%lu was %20.18e at a=%lu x=%lu\n",q,worsts[qq],as[qq],xs[qq]);
    }

  return(0);
}
