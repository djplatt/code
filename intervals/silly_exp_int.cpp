#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "../includes/int_double12.0.h"

#define epsilon ((double) 1e-6)

int_double int_epsilon(-epsilon,epsilon);

int_double myexp(double x)
{
  int_double res=1.0;
  long unsigned int n=1;
  int_double next_term=x;

  while(true)
    {
      if(fabs(next_term.right)<=epsilon)
	return(res+int_epsilon);
      res+=next_term;
      n++;
      next_term*=x;
      next_term/=n;
    }
}


int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("command line:- silly_exp <x>\n");
      exit(0);
    }

  _fpu_rndd();
  double x=atof(argv[1]);

  printf("  exp(%f)=%20.18e\nmyexp(%f)=",x,exp1(x),x);print_int_double_str("",myexp(x));
  return(0);
}
