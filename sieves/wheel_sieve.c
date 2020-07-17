#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"

#define debug printf("Reached line number %d.\n",__LINE__)

typedef long unsigned int ptype;
typedef __uint128_t bigint;

#define num_bits (sizeof(ptype)<<3)
#define log_num_bits (6)
#define all_ones ((num_bits<<1)-1)

ptype mask1[2*num_bits],mask[2*num_bits];

#define clear_prime(p,v) v[p>>(log_num_bits+1)]&=mask[p&all_ones]

#define primep(p,v) v[p>>(log_num_bits+1)]&mask1[(p&all_ones)]

#define NUM_SPOKES (7)
#define PROD_SPOKES (2*3*5*7*11*13*17)
#define PHI_SPOKES (2*4*6*10*12*16)
ptype primes[NUM_SPOKES-1]={3,5,7,11,13,17};

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


void num_prims(const ptype *source, const long unsigned int source_start, const long unsigned int source_len)
{
  long unsigned int i,count;
  for(i=1,count=0;i<=source_len;i++)
    if(primep(i,source))
      count++;
  printf("Pi(%lu)-Pi(%lu)=%lu\n",source_start+source_len-1,source_start,count);
}


inline void set_primes(ptype *vec, const ptype len)
{
  ptype i;
  for(i=0;i<=(len>>(log_num_bits+1));i++)
    {
      vec[i]=0;
      vec[i]=~vec[i];
    }
}

inline void do_masks()
{
  ptype i,j;
  for(i=1,j=1;i<2*num_bits;i+=2,j<<=1)
    {
      mask1[i]=j;
      mask[i]=0;
      mask[i]=~mask[i];
      mask[i]^=j;
      //printf("mask1[%lu]=%lX\nmask[%lu]=%lX\n",i,mask1[i],i,mask[i]);
    }
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
  ptype *target,ptr,p,two_p;
  ptype i,j,target_len=4000000;
  bigint target_start;
  ptype count=0,source_start,sqrt_end,offset;


  for(i=0,target_start=1;i<12;i++,target_start*=10);
  target_start++;
  if(!(target_start&1))
    fatal_error("target start must be odd.\n");

  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for target sieve.");

  set_primes(target,target_len);

  sqrt_end=sqrt(target_start+target_len);
  printf("sqrt_end set to %lu\n",sqrt_end);
  do_masks();
  do_spokes();
 
  //printf("target now %lX %lX\n",target[0],target[1]);

  // get rid of the spokes 3,5,7...
  for(i=0;i<NUM_SPOKES-1;i++)
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
	{
	  //printf("clearing ptr=%lu prime=%lu\n",ptr,target_start+ptr-1);
	  clear_prime(ptr,target);	
	  //printf("target now %lX %lX\n",target[0],target[1]);
	}
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
	{
	  //printf("clearing ptr=%lu prime=%lu\n",ptr,target_start+ptr-1);
	  clear_prime(ptr,target);	
	  //printf("target now %lX %lX\n",target[0],target[1]);
	}
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
	    {
	      //printf("clearing ptr=%lu prime=%lu\n",ptr,target_start+ptr-1);
	      clear_prime(ptr,target);	
	      //printf("target now %lX %lX\n",target[0],target[1]);
	    }
	}
    }

  for(i=1;i<target_len;i+=2)
    if(primep(i,target))
      {
	count++;
	//printf("%lu is prime\n",i+target_start-1);
      }
  printf("Pi(");print_bigint(target_start+target_len-1);printf(")-Pi(");
  print_bigint(target_start);printf(")=%lu\n",count);
  return(0);
}
