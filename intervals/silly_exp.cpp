#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#define epsilon ((double) 1e-6)

double myexp(double x)
{
  double res=1.0;
  long unsigned int n=1;
  double next_term=x;

  while(true)
    {
      if(fabs(next_term)<=epsilon)
	return(res);
      res+=next_term;
      n++;
      next_term*=x;
      next_term/=(double) n;
    }
}


int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("command line:- silly_exp <x>\n");
      exit(0);
    }

  double x=atof(argv[1]);

  printf("  exp(%f)=%20.18e\nmyexp(%f)=%20.18e\n",x,exp(x),x,myexp(x));
  printf("     diff=%20.18e\n",fabs(exp(x)-myexp(x)));
  return(0);
}
