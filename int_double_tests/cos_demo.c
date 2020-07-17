#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "math.h"
#define TEST_FN mycos
#define TEST_FN_STR "fcos " 
#define MPFI_TEST_FN mpfi_cos
#define GCC_LIB_FN cos
#define GCC_LIB_FN_STR "cos  "
#define GCC_LIB_FN_L cosl
#define GCC_LIB_FN_STR_L "cosl "

#define TRUE (1)
#define FALSE (0)

#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  :: "m" (cw))

short int old_cw,new_cw;

void rndd()
{
  _fpu_getcw(old_cw);
  //printf("old cw=%x\n",old_cw);
  new_cw=old_cw&0x00FF;
  new_cw=new_cw|0x0700;
  //printf("rndd cw=%X\n",new_cw);
  _fpu_setcw(new_cw);
}

void rndu()
{
  _fpu_getcw(old_cw);
  //printf("old cw=%x\n",old_cw);
  new_cw=old_cw&0x00FF;
  new_cw=new_cw|0x0B00;
  //printf("rndu cw=%X\n",new_cw);
  _fpu_setcw(new_cw);
}

void rndn()
{
  _fpu_getcw(old_cw);
  //printf("old cw=%x\n",old_cw);
  new_cw=old_cw&0x00FF;
  //new_cw=new_cw|0x0300;
  //printf("rndn cw=%X\n",new_cw);
  _fpu_setcw(new_cw);
}


// print the hex representation of a 64 bit float
void print_double(double x)
{
  // comment out next two lines if compiler does not
  // support %A,%a
  //printf("%A\n",x);
  //return;
  int *xi;
  xi=(int *) &x;
  printf("%X %X\n",xi[1],xi[0]);
}

int ulp_diff (double x1, double x2)
{
  int *i1,*i2;
  i1=(int *)&x1;
  i2=(int *)&x2;
  if(i1[1]!=i2[1])
    {
      printf("Huge difference. Exiting.\n");
      exit(0);
    }
  return(i1[0]-i2[0]);
}

double mycos(double x)
{
  double res;
  // load the FPU with x, cos it, store to res
  __asm("fldl %1\n\t"
	"fcos\n\t"
	"fstpl %0"
	: "=m" (res)
	: "m" (x)
	: "st");
  return(res);
}

double mysin(double x)
{
  double res;
  // load the FPU with x, cos it, store to res
  __asm("fldl %1\n\t"
	"fsin\n\t"
	"fstpl %0"
	: "=m" (res)
	: "m" (x)
	: "st");
  return(res);
}



int test_val(double x)
{
  double resd,resu,resn;
  printf("testing argument x=%30.28e = ",x);print_double(x);
  double gcc=GCC_LIB_FN(x);
  printf("Gcc library   %s produced %30.28e = ",GCC_LIB_FN_STR,gcc);print_double(gcc);
  long double gccl=GCC_LIB_FN_L(x);
  gcc=gccl;
  printf("Gcc L library %s produced %30.28e = ",GCC_LIB_FN_STR_L,gcc);print_double(gcc);
  rndd();
  resd=TEST_FN(x);
  printf("rounding down %s produced %30.28e = ",TEST_FN_STR,resd);print_double(resd);
  rndu();
  resu=TEST_FN(x);
  printf("rounding up   %s produced %30.28e = ",TEST_FN_STR,resu);print_double(resu);
  rndn();
  resn=TEST_FN(x);
  printf("rounding near %s produced %30.28e = ",TEST_FN_STR,resn);print_double(resn);
  return(0);
}

int main(int argc, char ** argv)
{
  double x=(double) 884279718950133.0;
  int i;
  for(i=0;i<49;i++)
	  x/=2.0;
  test_val(x);
    
  return(0);
}
