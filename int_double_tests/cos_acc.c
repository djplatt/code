#include "stdio.h"
#include "stdlib.h"
#include "assert.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "math.h"
#include "crlibm.h"
#define PREC (100)
#define NUM_STEPS (100000)
#define MIN_X ((double) -M_PI)
#define MAX_X ((double) M_PI)
#define DELTA ((double) (MAX_X-MIN_X)/(double) NUM_STEPS)
#define TEST_FN cos_rn
#define TEST_FN_STR "cos_rn " 
#define MPFI_TEST_FN mpfi_cos
#define TRUE (1)
#define FALSE (0)

#define _fpu_getcw(cw) asm ("fnstcw %0" : "=m" (cw))
#define _fpu_setcw(cw) asm ("fldcw %0"  :: "m" (cw))
#define __SSE_getcw(cw) asm ("stmxcsr %0" : "=m" (cw))
#define __SSE_setcw(cw) asm ("ldmxcsr %0" :: "m" (cw))

int old__SSE_cw,new__SSE_cw;

void rndd()
{
	__SSE_getcw(old__SSE_cw);
	new__SSE_cw=old__SSE_cw&0x1FBF; // zeros FTZ(15),RC(14,13) and DAZ(6) bits
	new__SSE_cw=new__SSE_cw|0x2000; // sets RC to Round Down (14=0, 13=1)
	__SSE_setcw(new__SSE_cw);
}

double my_sqrt(double x)
{
  double res;
  __asm("fldl %1\n\t"
	"fsqrt\n\t"
	"fstpl %0"
	: "=m" (res)
	: "m" (x));
  return(res);
}

void rndn()
{
  __SSE_setcw(old__SSE_cw);
}

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


  mpfi_t mx,mres;
  mpfr_t mleft,mright;


int test_val(double x)
{
  int diff1,diff2;
  double resn,left,right;
  resn=TEST_FN(x);
  mpfi_set_d(mx,x);
  MPFI_TEST_FN(mres,mx);
  mpfi_get_left(mleft,mres);
  mpfi_get_right(mright,mres);
  left=mpfr_get_d(mleft,GMP_RNDN);
  right=mpfr_get_d(mright,GMP_RNDN);
  if(left!=right)
    {
      printf("MPFI produced too wide a result.\n");
      exit(0);
    }
  if(resn!=left)
    {
      printf("error with x= %50.48e ",x);print_double(x);
      printf("mpfi returned %50.48e ",left);print_double(left);
      printf("rndn returned %50.48e ",resn);print_double(resn);
      return(1);
    }
  return(0);
}



int main(int argc, char ** argv)
{
  // check we have sensible data types
  assert(sizeof(long long)>=8);
  assert(sizeof(double)==8);
  assert(sizeof(int)==4);
  crlibm_init();
  rndd(); // set SSE rounding mode to down
  mpfr_set_default_prec(PREC);
  double left,right,x;
  mpfr_init(mleft);
  mpfr_init(mright);
  mpfi_init(mx);
  mpfi_init(mres); 
  unsigned int *xi,old;
  int count,err,err_count;
  x=(double) 884279718950133.0/(double) ((long long) 1<<49);
  printf("x=%A\n",x);
  test_val(x);
  double delta=1.0/(double) (1<<15);

  for(err_count=0,count=0,x=MIN_X;
      count<NUM_STEPS;
      x+=DELTA,count++,err_count+=test_val(x));
  
  printf("Total errors found=%d\n",err_count);
  
  return(0);
}
