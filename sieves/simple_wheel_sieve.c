#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"

#define debug printf("Reached line number %d.\n",__LINE__)
#define TRUE (1==1)
#define FALSE (1==0)

typedef long unsigned int ptype;
typedef __uint128_t bigint;

#define NUM_SPOKES (7)
#define PROD_SPOKES (2*3*5*7*11*13*17)
#define PHI_SPOKES (2*4*6*10*12*16)
ptype primes[NUM_SPOKES]={2,3,5,7,11,13,17};

ptype spokes[PHI_SPOKES];


inline ptype gcd (ptype a, ptype b)
// Euclid algorithm gcd
// best if a<=b ?
{
  unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};


inline void print_bigint(bigint i)
{
  if(i<10)
    printf("%1lu",(long unsigned int) i);
  else
    {
      print_bigint(i/10);
      printf("%lu",(long unsigned int) (i%10));
    }
}

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}

inline void do_spokes()
{
  ptype i,j;
  for(i=1,j=0;i<=PROD_SPOKES;i++)
    if(gcd(i,PROD_SPOKES)==1)
      spokes[j++]=i;
}


int main()
{
  unsigned char *target;
  ptype ptr,p,two_p;
  ptype i,j,target_len=40000;
  bigint target_start;
  ptype count=0,source_start,sqrt_end,offset;


  for(i=0,target_start=1;i<19;i++,target_start*=10);

  if(!(target=(unsigned char *) malloc(sizeof(unsigned char)*target_len)))
    fatal_error("Failed to allocate memory for target sieve.");

  for(i=0;i<target_len;i++)
    target[i]=TRUE;
  sqrt_end=sqrt(target_start+target_len);
  printf("sqrt_end set to %lu\n",sqrt_end);
  do_spokes();
 
  for(i=0;i<NUM_SPOKES;i++)
    {
      p=primes[i];
      two_p=p<<1;
      ptr=target_start%p;
      if(ptr==0)
	ptr=p;
      ptr=p-ptr+1;
      if(!(ptr&1))
	ptr+=p;
      for(;ptr<target_len;ptr+=two_p)
	target[ptr]=FALSE;	
    }

  for(i=1;i<PHI_SPOKES;i++) // skip p=1
    {
      p=spokes[i];
      //printf("clearing p=%lu\n",p);
      two_p=p<<1;
      ptr=target_start%p;
      if(ptr==0)
	ptr=p;
      ptr=p-ptr+1;
      if(!(ptr&1))
	ptr+=p;
      for(;ptr<target_len;ptr+=two_p)
	target[ptr]=FALSE;
    }
  for(offset=PROD_SPOKES;offset<=sqrt_end;offset+=PROD_SPOKES)
    {
      for(i=0;i<PHI_SPOKES;i++) // skip p=1
	{
	  p=spokes[i]+offset;
	  //printf("clearing p=%lu\n",p);
	  two_p=p<<1;
	  ptr=target_start%p;
	  if(ptr==0)
	    ptr=p;
	  ptr=p-ptr+1;
	  if(!(ptr&1))
	    ptr+=p;
	  for(;ptr<target_len;ptr+=two_p)
	    target[ptr]=FALSE;
	}
    }

  for(i=1;i<target_len;i+=2)
    if(target[i])
	count++;
  printf("Pi(");print_bigint(target_start+target_len-1);printf(")-Pi(");
  print_bigint(target_start);printf(")=%lu\n",count);
  return(0);
}
