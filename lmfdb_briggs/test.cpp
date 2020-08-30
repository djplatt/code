#include "../includes/int_double14.0.h"
#include "inttypes.h"

int_double pow_ui(const int_double x, const uint64_t pow)
{
  uint64_t ppow=pow;
  int_double res=1;
  int_double sq=x;
  while(ppow>0)
    {
      if(ppow&1)
	res*=sq;
      ppow>>=1;
      sq*=sq;
    }
  return res;
}

int main()
{
  _fpu_rndd();
  int_double p=2;
  int_double lp=log(p);
  print_int_double_str("",log1p((p-1)/(pow(p,44)-p))/lp);
  print_int_double_str("",log1p((p-1)/(pow_ui(p,44)-p))/lp);

  return 0;
}
