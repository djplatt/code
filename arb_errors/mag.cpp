/*

investigate radius in ARB being wider than expected.

Here is what I expect

1 +/- 1.97521522..e-31 (2^-102)
1.0 +/- 0

Here is what I get

1 +/- 1.9722e-31
1.00000000186264514923095703125 +/- 0


Note that the radius is larger than necessary, even though it is an 
exact power of 2

*/

#include "arb.h"
#define SHIFT (102) // our delta is going to be 2^-SHIFT
int main()
{
  arb_t x,err;
  arb_init(x);
  arb_set_ui(x,1);
  arb_init(err);
  arb_set_ui(err,1);
  arb_mul_2exp_si(err,err,-SHIFT); // 2^-SHIFT 
  arb_add_error(x,err); // 1 +/- 2^-SHIFT
  arb_printd(x,200);printf("\n"); // print it
  arb_get_rad_arb(x,x); // should be (?) 2^-SHIFT
  arb_mul_2exp_si(x,x,SHIFT);  // should now be 2^0
  arb_printd(x,200);printf("\n"); // print it

  return 0;
}
