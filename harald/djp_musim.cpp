#include <stdlib.h>
#include <stdio.h>
#include "inttypes.h"
#include "../includes/int_double12.0.h"
//#include "../int_double/int_double.h"
#include "primesieve.hpp"

int8_t *vec;
uint64_t st,en,wid;

inline uint64_t get_start(uint64_t n, uint64_t st)
{
  uint64_t r=st%n;
  if(r==0)
    return 0;
  else
    return n-r;
}

void callback(uint64_t p)
{
  uint64_t ptr=get_start(p,st);
  for(uint64_t ptr1=ptr;ptr1<wid;ptr1+=p)
    vec[ptr1]=-vec[ptr1];
  uint64_t p2=p*p;
  ptr=get_start(p2,st);
  for(uint64_t ptr1=ptr;ptr1<wid;ptr1+=p2)
    vec[ptr1]=0;
}

void sieve()
{
  for(uint64_t i=0;i<wid;i++)
    vec[i]=1;
  primesieve::callback_primes(2,en,callback);
}

int main(int argc, char **argv)
{

  if(argc!=5)
    return 0;

  _fpu_rndd();

  st=atol(argv[1]);
  wid=atol(argv[2]);
  en=st+wid;
  uint64_t its=atol(argv[3]);
  FILE *outfile=fopen(argv[4],"w");
  fwrite(&st,sizeof(uint64_t),1,outfile);
  fwrite(&wid,sizeof(uint64_t),1,outfile);
  fwrite(&its,sizeof(uint64_t),1,outfile);
  vec=(int8_t *)malloc(sizeof(int8_t)*wid);
  primesieve::set_sieve_size(32);

  for(uint64_t it=0;it<its;it++,st+=wid,en+=wid)
    {
      sieve();

      int_double max=0.0;
      int64_t max_n=0;
      int_double max_M=0;
      int_double min=0.0;
      int64_t min_n=0;
      int_double min_M=0;
      int_double M=0.0;

      for(int64_t i=0,n=st;i<wid;i++,n++)
	{
	  if(vec[i]==0) continue;
	  int_double dn=n;
	  M+=vec[i]/dn;
	  int_double del=M*sqrt(dn);
	  if(vec[i]==-1)
	    {
	      if(del.left<min.left)
		{min=del;min_n=n;min_M=M;}
	    }
	  else
	    if(del.right<max.right)
	      {max=del;max_n=n;max_M=M;}
	}
      
      fwrite(&min_n,sizeof(int64_t),1,outfile);
      fwrite(&min_M,sizeof(int_double),1,outfile);

      fwrite(&max_n,sizeof(int64_t),1,outfile);
      fwrite(&min_M,sizeof(int_double),1,outfile);

      fwrite(&M,sizeof(int_double),1,outfile);


      printf("Running from %lu to %lu we have \n",st,en);
      print_int_double_str("M = ",M);
      printf("min M*sqrt was to n=%lu",min_n);print_int_double_str(" ",min_M);
      printf("max M*sqrt was to n=%lu",max_n);print_int_double_str(" ",max_M);

      //for(uint64_t i=0;i<25;i++)	printf("mu(%lu)=%d\n",st+i,vec[i]);
    }

  fclose(outfile);
  return 0;
}
