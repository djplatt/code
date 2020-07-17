#include "inttypes.h"
#include "stdio.h"
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <assert.h>

// 128 bit unsigned built in to GCC (efficiency?)
typedef __int128_t bigint;

#define H (1000000)
#define W (28)

bigint cubes[H];

int main()
{
  for(int i=1;i<H;i++)
    {
      cubes[i]=i;cubes[i]*=i;cubes[i]*=i;
    }
  for(uint64_t x=1;x<H;x++)
    {
      bigint x3=cubes[x];
      uint64_t z=(uint64_t)((double)x*(double) 1.259921049894873);
      if(z>=H)
	exit(0);
      bigint z3=cubes[z];
      for(uint64_t y=x;y<H;y++)
	{
	  bigint y3=cubes[y];
	  while(z3-W<y3+x3)
	    {
	      z++;
	      if(z>=H)
		{
		  y=H;
		  break;
		}
	      z3=cubes[z];
	    }
	  if(z3+W==y3+x3)
	    {
	      printf("Solution found %lu^3+%lu^3=%lu^3+%lu.\n",x,y,z,W);
	      exit(0);
	    }
	}
    }
  return(0);
}

	  
