#include "../includes/int_double12.0.h"
#include "inttypes.h"

#define LIMIT ((double) 1e-100)
void measure (int_double (*eval_fn) (int_double), int_double interval, int_double limit, int_double *totals)
{
  int_double res=(*eval_fn)(interval);
  //printf("f(");print_int_double(interval);print_int_double_str(") returned ",res);
  if(res.left>=-limit.right) // result is entirely > limit
    {
      //printf("> limit\n");
      totals[0]+=int_double(-interval.right)-interval.left;
      return;
    }
  if(-res.right<=limit.left) // result is entirely < limit
    {
      //printf("< limit\n");
      totals[1]+=int_double(-interval.right)-interval.left;
      return;
    }
  //printf("<> limit\n");
  double mp=(interval.left-interval.right)/2.0;
  //printf("%30.28g %30.28g %30.28g\n",mp,interval.left,-interval.right);
  if((mp-interval.left<=LIMIT)||(-interval.right-mp<=LIMIT))
    {
      totals[2]+=int_double(-interval.right)-interval.left;
      return;
    }
  measure(eval_fn,int_double(interval.left,mp),limit,totals);
  measure(eval_fn,int_double(mp,-interval.right),limit,totals);
}

int_double x2(int_double x)
{
  return(sqr(x));
}

uint64_t L;

int_double GL(int_double x)
{
  int_complex sum=c_zero,term;
  int_double x2=x*2; // start with exp(2 pi i 2^0 x)

  for(uint64_t j=0;j<L;j++,x2*=2)
    {
      sin_cospi(x2,&term.imag,&term.real);
      sum+=term;
    }
  return(mod(sum));
}

int_double TL(int_double x)
{
  int_complex sum=c_zero,term;
  int_double x2=x*4; // start with exp(2 pi i 2^1 x)

  for(uint64_t j=0;j<L;j++,x2*=2)
    {
      sin_cospi(x2,&term.imag,&term.real);
      sum+=term;
    }
  return(mod(sum));
}


int main(int argc, char **argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <lambda> <L>.\n",argv[0]);
      return(1);
    }
  L=atol(argv[2]);
  _fpu_rndd();
  int_double ms[3];
  ms[0]=0;ms[1]=0;ms[2]=0;
  //measure(*GL,int_double(0,1.0),int_double(atof(argv[1]))*L,ms);
  measure(*TL,int_double(0,1.0),int_double(atof(argv[1]))*L,ms);
  printf("%f ",atof(argv[1]));
  //print_int_double_str("Definately larger  =    ",ms[0]);
  print_int_double_str("Not smaller        =    ",1-ms[1]);
  //print_int_double_str("Definately smaller =    ",ms[1]);
  //print_int_double_str("Not larger         =    ",1-ms[0]);
  //print_int_double_str("Indeterminate      =    ",ms[2]);
  //print_int_double_str("Total              =    ",ms[0]+ms[1]+ms[2]); 
  //print_int_double_str("N^(-109/154)       =    ",exp(log(int_double(12<<(L+1)))*int_double(-109)/154));
  return(0);
}
  
