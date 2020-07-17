#include "stdio.h"
#include "mpfr.h"
#include "math.h"

#define LIM 1000000
#define PREC 53


int main()
{

  int i;
 double y,z;

 //  mpfr_set_default_prec(PREC);
 //  mpfr_init(y);


 z=1.0;
  for(i=0;i<LIM;i++)
    {
      //   mpfr_set_ui(y,i,GMP_RNDN);
      // mpfr_div_ui(y,y,LIM,GMP_RNDN);
      y=sin(z);
      z=y;
    };
  printf("%f\n",y);
  return (0);
  }


