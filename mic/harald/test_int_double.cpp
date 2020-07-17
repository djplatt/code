#include "int_double12.0.h"

int main()
{
  _fpu_rndd();
  int_double x=1;
  x/=3;
  print_int_double_str("1/3=",x);
  return(0);
}
