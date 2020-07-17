#include "stdlib.h"
#include "stdio.h"

inline bool smaller(long unsigned int this_q,long unsigned int this_a,long unsigned int q, long unsigned int a)
{
  this_q=this_q+this_a/q;
  this_a=this_a%q;
  if(this_q<q)
    return(true);
  if(this_q>q)
    return(false);
  return(this_a<a);
}

inline bool test_it(long unsigned int q, long unsigned int a)
{
  long unsigned int target=q*q+a;
  long unsigned int q1=(q*3)/2;
  long unsigned int a1=(a*3+1)/2;
  while(true)
    {
      if(smaller(q1,a1,q,a))
	{
	  //printf("Got under q=%ld a=%ld with q1=%ld a1=%ld.\n",q,a,q1,a1);
	  return(true);
	}
      if(q1&1)
	return(false);
      if(a1&1)
	{
	  if((q1>6148914691236517205)||(a1>6148914691236517204))
	    {
	      printf("Overflow happened.\n");
	      exit(0);
	    }
	  q1=q1*3;
	  a1=a1*3+1;
	}
      else
	{
	  q1>>=1;
	  a1>>=1;
	}
    }
}
	  
       
int main()
{
  long unsigned int a,q;
  long unsigned int worst_a=0;
  for(q=4;q<=100000000000;q+=4)
    {
      for(a=3;a<q;a+=4)
	{
	  if(!test_it(q,a))
	    {
	      if(a>worst_a)
		{
		  worst_a=a;
		  printf("Modulus %ld failed at a=%ld.\n",q,a);
		}
	      break;
	    }
	}
      if(a>=q)
	{
	  printf("Good Modulus found, q=%ld\n",q);
	  exit(0);
	}
    }
  return(0);
}
