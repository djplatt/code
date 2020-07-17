#include "stdio.h"
#include "math.h"
#include "assert.h"

// print the hex representation of a 64 bit float
void print_double(double x)
{
  // comment out next two lines if compiler does not
  // support %A,%a
  printf("%A\n",x);
  return;
  int *xi;
  xi=(int *) &x;
  printf("%X %X\n",xi[1],xi[0]);
}

int main()
{
  // check we have sensible data types
  assert(sizeof(long long)>=8);
  assert(sizeof(double)==8);
  assert(sizeof(int)==4);

  long long num=884279718950133;
  long long denom=1;
  double best_res;
  int *best_resi=(int *) &best_res;
  best_resi[1]=0x3dda15c1;
  best_resi[0]=0x1a626331;

  denom<<=49;
  double res,x=(double) num/(double) denom;
  // x is now quite close to Pi/2 
  printf("x=%60.58e\n =",x);print_double(x);

  // load the FPU with x, cos it, store to res
  __asm("fldl %1\n\t"
	"fcos\n\t"
	"fstpl %0"
	: "=m" (res)
	: "m" (x)
	:);

  printf("using asm\ncos(x)=%60.58e\n      =",res);print_double(res);

  // use c++ built in cos
  res=cos(x);
  printf("using c++ RTL\ncos(x)=%60.58e\n      =",res);print_double(res);

  printf("we wanted\ncos(x)=%60.58e\n      =",best_res);print_double(best_res);


  return(0);
}
