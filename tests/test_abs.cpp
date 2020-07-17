#include "../includes/int_double14.0.h"

int main()
{
  _fpu_rndd();

  int_double x;
  x.left=-0.00004;
  x.right=-0.001;

  print_int_double_str("x=",x);
  print_int_double_str("|x|=",abs(x));
  return 0;
}
