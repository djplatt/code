#include "stdio.h"
#include "math.h"

void process_prime(long int c, long int p)
{
  printf("c=%ld p=%ld.\n",c,p);
}


long int gulich(long int c1, long int c2)
{
  long int L=0;
  for(long int c=c1;c<=c2;c++)
    {
      bool not_found=true;
      double sqrt_x=sqrt(1.0+6.0*(double)c);
      long int x1=(c+1)/5,x2=(long int) (sqrt_x+1.0)/6;
      //printf("c=%lu x1=%ld x2=%ld\n",c,x1,x2);
      //printf("Entering step 2(a)\n");
      for(long int x=-x1;x<=-x2;x++)
	{
	  //printf("x=%ld\n",x);
	  if((c-x)%(6*x+1)==0)
	    {
	      L++;
	      not_found=false;
	      break;
	    }
	}
      if(not_found) // do step 2(b)
	{
	  not_found=true;
	  for(long int x=x1;x<=(long int) (sqrt_x-1.0)/6;x+=1)
	    {
	      if((c-x)%(6*x+1)==0)
		{
		  L++;
		  not_found=false;
		  break;
		}
	    }
	  if(not_found)
	    process_prime(c,6*c+1);
	}
      // step 2(c)
      x1=(c+1)/7;
      not_found=true;
      for(long int x=-x1;x<=-1;x++)
	{
	  if((-c-x)%(6*x+1)==0)
	    {
	      L++;
	      not_found=false;
	      break;
	    }
	}
      if(not_found)
	process_prime(-c,-6*c+1);
    }
  return(L);
}

int main()
{
  long int c1=1000000000,c2=c1+100;
  gulich(c1,c2);
  return(0);
}
