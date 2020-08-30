#include "malloc.h"
#include "stdio.h"
#include "inttypes.h"
#include <cstdlib>

#define N (18LL<<33) // should be a multiple of 18*8

uint8_t *ns; // our table of potential crocker numbers

// and with masks[i] to clear bit [i], i=0..7
uint64_t masks[]={255-1,255-2,255-4,255-8,255-16,255-32,255-64,255-128};

// cross out the number n
// if n!=0 mod(18) skip
void crossout(uint64_t n)
{
  uint64_t res=n%18;
  uint64_t quo=n/18;
  if(res==0) // its divisible by 18 so zero the corresponding bit
    ns[quo>>3]&=masks[quo&7];
}

void do_crocker()
{
  printf("Setting up table of potential Crocker Numbers.\n");
  for(uint64_t i=0,j=0;i<N;i+=18*8,j++)
    ns[j]=255;
  ns[0]=254; // can't do 1 as sum of 2 squares
  printf("Crossing out sums of squares plus 0,1,2 powers of 2.\n");
  for(uint64_t a=1;;a++)
    {
      //printf("a=%lu\n",a);
      uint64_t a2=a*a;
      if(a2>N) break;
      for(uint64_t b=a;;b++)
	{
	  uint64_t b2=b*b;
	  uint64_t a2b2=a2+b2;
	  if(a2b2>N) break;
	  crossout(a2b2);
	  for(uint64_t alpha=1;;alpha+=alpha)
	    {
	      uint64_t a2b2a=a2b2+alpha;
	      if(a2b2a>N) break;
	      crossout(a2b2a);
	      for(uint64_t beta=alpha+alpha;;beta+=beta)
		{
		  uint64_t a2b2ab=a2b2a+beta;
		  if(a2b2ab>N) break;
		  crossout(a2b2ab);
		}
	    }
	}
    }
  for(uint64_t i=0,j=0;i<N;i+=18*8,j++)
    if(ns[j]>0)
      printf("Found one %u at n=%lu-%lu\n",ns[j],i,i+7);
}

int main()
{
  printf("Checking for Crocker numbers <=%lu\n",N);
  ns=(uint8_t *)malloc(sizeof(uint8_t)*N/18/8);
  do_crocker();
  printf("Finished.\n");
  return(0);
}
