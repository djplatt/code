#include "../includes/int_double13.0.h"

int main()
{
  _fpu_rndd();
  print_int_double_str("log(1-1e-10)=",log1p(int_double(-1.0)/10000000000.0));
  print_int_double_str("log(1-1e-10)=",log(1.0-int_double(1.0)/10000000000.0));
  return 0;
}
