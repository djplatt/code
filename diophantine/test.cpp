#include "stdio.h"
#include "math.h"
#define N (29)
#define MINX (4)
#define MAXX (100000)
#define CUBE_ROOT_2 ((double) 1.25992104989487)

inline int test_cube(long long int tgt)
{
  double approx=exp(log(tgt)/3.0)*2.0;
  long long int z=(long long int) approx;
  if(z&1)
    z=(z>>1)+1;
    else
      z>>=1;
  return(z*z*z==tgt);
}

int main()
{
  long long int x,y,z,x3,y3,z3;
  printf("using integers of length %ld\n",sizeof(x));
  printf("searching x^3+y^3+z^3=%d for solutions |x|<=|y|<=|z|and |x| in[%d,%d]\n",N,MINX,MAXX);
  for (x=MINX;x<=MAXX;x++)
    {
      //if((x&8191)==0) printf("x=%lld\n",x);
      x3=x*x*x;
      for(y=(long long int) ceil(((double) x)/CUBE_ROOT_2);y<x;y++)
	{
	  y3=y*y*y;
	  //z3=x3+y3-N; // assume x,y>=0
	  //if(z3<0) z3=-z3;
	  //if(z3<=y3)
	  //if(test_cube(z3))
	  //printf("solution x=%lld y=%lld\n",x,y);
	  z3=x3-y3-N; // must be at least 37-N
	  //printf("%lld %lld %lld\n",x3,y3,z3);
	  if(z3<=y3)
	    if(test_cube(z3))
	      printf("solution x=%lld y=%lld\n",x,-y);
	  z3+=(N<<1);
	  //printf("%lld %lld %lld\n",x3,y3,z3);
	  if(z3<=y3)
	    if(test_cube(z3))
	      printf("solution x=%lld y=%lld\n",-x,y);
	}
    }
  printf("Completed.\n");
}
