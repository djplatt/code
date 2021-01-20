#include "inttypes.h"
#include "stdbool.h"
#include "arb.h"

void f(arb_t res, arb_t z, int64_t prec)
{
  static bool init=false;
  static arb_t lz,one_z,tmp1,tmp2,tmp3,tmp4;
  if(!init)
    {
      init=true;
      arb_init(lz);
      arb_init(one_z);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(tmp4);
    }
  arb_set_ui(tmp1,1);
  arb_sub(one_z,tmp1,z,prec); // (1-z)
  arb_log(lz,z,prec); // log z
  arb_mul_2exp_si(tmp1,one_z,2); // 4(1-z)
  arb_mul_ui(tmp2,z,5,prec);
  arb_add_ui(tmp3,tmp2,1,prec); // (1+5z)
  arb_mul(tmp2,tmp3,lz,prec);
  arb_div(tmp3,tmp2,tmp1,prec); // [(1+5z) log z]/4(1-z)
  arb_mul_ui(tmp2,one_z,3,prec); // 3(1-z)
  arb_add_ui(tmp1,z,1,prec);
  arb_mul(tmp4,tmp1,lz,prec); // (1+z) log z
  arb_add(tmp1,tmp4,tmp2,prec); // 3(1-z)+(1+z)log z
  arb_mul(tmp2,z,lz,prec); // z log z
  arb_div(tmp4,tmp2,tmp1,prec);
  arb_neg(tmp4,tmp4);
  arb_log(tmp1,tmp4,prec);
  arb_add(tmp4,tmp3,tmp1,prec);
  arb_set_d(tmp1,1.5);
  arb_add(res,tmp4,tmp1,prec);
  return;
}

#define L2_STEP (29)
#define PREC (100)

int main()
{
  arb_t err,z,res;
  arb_init(err);arb_init(z);arb_init(res);

  arb_set_d(z,0.5);
  f(res,z,PREC);
  printf("f(1/2) -> ");arb_printd(res,20);printf("\n");
  
  arb_set_d(err,1);
  arb_mul_2exp_si(err,err,-L2_STEP-1);
  
  uint64_t i;
  bool been_ok=false;
  for(i=0;i<(1LU<<L2_STEP);i++)
    {
      arb_set_d(z,0.5+(double) i);
      arb_mul_2exp_si(z,z,-L2_STEP);
      arb_add_error(z,err);
      f(res,z,PREC);
      if(!arb_is_positive(res))
	{
	  if(been_ok)
	    {
	      printf("Failed to be positive on interval %lu/%lu +/- 1/%lu",i+i+1,1LU<<(L2_STEP+1),(1LU<<L2_STEP+1));
	      printf(" result was ");
	      arb_printd(res,10);
	      printf("\n");
	      return 0;
	    }
	}
      else
	{
	  if(!been_ok) // first time it worked
	    {
	      printf("Successfully positive on interval %lu/%lu +/- 1/%lu",i+i+1,1LU<<(L2_STEP+1),(1LU<<L2_STEP+1));
	      printf(" result was ");
	      arb_printd(res,10);
	      printf("\n");
	      been_ok=true;
	    }
	}
	
    }

  arb_clear(z);arb_clear(err);arb_clear(res);
  return 0;
}
