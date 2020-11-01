#ifndef WIN_ZETA
#define WIN_ZETA

#include "inttypes.h"
//#define OP_ACC (101) // Output is correct to +/- 2^(-OP_ACC-1)

// defined so that POS&NEG == 0 but no other combination
typedef int sign_t;
#define POS (1)
#define NEG (2)
#define UNK (3)

typedef int dir_t;
#define UP (1)
#define DOWN (2)

inline sign_t sign(arb_t x)
{
  if(arb_is_negative(x))
    return(NEG);
  if(arb_is_positive(x))
    return(POS);
  return(UNK);
}

inline dir_t dir(arb_ptr left, arb_ptr right)
{
  if(arb_gt(left,right))
    return(DOWN);
  if(arb_lt(left,right))
    return(UP);
  return(UNK);
}

inline sign_t sign(acb_t x)
{
  if(arb_is_negative(acb_realref(x)))
    return(NEG);
  if(arb_is_positive(acb_realref(x)))
    return(POS);
  return(UNK);
}

inline dir_t dir(acb_ptr left, acb_ptr right)
{
  if(arb_gt(acb_realref(left),acb_realref(right)))
    return(DOWN);
  if(arb_lt(acb_realref(left),acb_realref(right)))
    return(UP);
  return(UNK);
}


double arb_get_d(arb_t x)
{
  static bool init=false;
  static arb_t tmp1;
  static mag_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      mag_init(tmp);
    }
  arb_get_mid_arb(tmp1,x);
  arb_get_mag(tmp,tmp1);
  if(arb_is_negative(tmp1))
    return(-mag_get_d(tmp));
  else
    return(mag_get_d(tmp));
}

void arb_mul_d(arb_t x, arb_t y, double d, int64_t prec)
{
  static bool init=false;
  static arb_t dd;
  if(!init)
    {
      arb_init(dd);
      init=true;
    }
  arb_set_d(dd,d);
  arb_mul(x,y,dd,prec);
}


void arb_div_d(arb_t x, arb_t y, double d, int64_t prec)
{
  static bool init=false;
  static arb_t dd;
  if(!init)
    {
      arb_init(dd);
      init=true;
    }
  arb_set_d(dd,d);
  arb_div(x,y,dd,prec);
}

void arb_add_d(arb_t x, arb_t y, double d, int64_t prec)
{
  static bool init=false;
  static arb_t dd;
  if(!init)
    {
      arb_init(dd);
      init=true;
    }
  arb_set_d(dd,d);
  arb_add(x,y,dd,prec);
}

void arb_sub_d(arb_t x, arb_t y, double d, int64_t prec)
{
  static bool init=false;
  static arb_t dd;
  if(!init)
    {
      arb_init(dd);
      init=true;
    }
  arb_set_d(dd,d);
  arb_sub(x,y,dd,prec);
}

int arb_cmp(arb_t x, arb_t y)
{
  if(arb_gt(x,y))
    return(1);
  if(arb_lt(x,y))
    return(-1);
  if(arb_eq(x,y))
    return(0);
  printf("fatal error in arb_cmp. Arguments overlap. Exiting.\n");
  exit(0);
}

#endif 
