#include "inttypes.h"
#include "stdio.h"
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <assert.h>
#include <math.h>

// 128 bit unsigned built in to GCC (efficiency?)
typedef __int128_t bigint;

#define H (1000000)
#define W (30)

uint64_t cubes[H];

int main()
{
  for(int i=1;i<H;i++)
    {
      cubes[i]=i;cubes[i]*=i;cubes[i]*=i;
    }
  for(uint64_t x=1;x<H;x++)
    {
      for(uint64_t y=x;y<H;y++)
	{
	  uint64_t target=cubes[x]+cubes[y]-W;
	  double zd=exp(log((double)target)/3.0);
	  uint64_t z=floor(zd+0.5);
	  if(z>=H)
	    {
	      y=H;
	      break;
	    }
	  if(cubes[z]==target)
	    {
	      bigint x3=x;x3*=x;x3*=x;
	      bigint y3=y;y3*=y;y3*=y;
	      bigint z3=z;z3*=z;z3*=z;
	      if(x3+y3==z3+W)
		{
		  printf("Solution found %lu^3+%lu^3=%lu^3+%ld.\n",x,y,z,W);
		  exit(0);
		}
	      else
		printf("Solution found mod 2^64.\n");
	    }
	}
    }
  return(0);
}
