#include "arb.h"

int main()
{
  arb_t x,err;
  arb_init(x);
  arb_init(err);
  arb_set_d(x,(0.01-0.00004)/2);
  arb_set_d(err,(0.01+0.00004)/2);
  arb_add_error(x,err);
  arb_printd(x,20);printf("\n");
  arb_abs(x,x);
  arb_printd(x,20);printf("\n");
  return 0;
}
