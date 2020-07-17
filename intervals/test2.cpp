#include <stdio.h>
#include <math.h>

using namespace std;

int main()
{
  double foo,bar;
  unsigned int i;

  foo=0.999999999999999;

  for(i=0;i<10000000;i++)
    {
      bar=sin(foo);
      bar=cos(foo);
    };
    
  return(0);
};
