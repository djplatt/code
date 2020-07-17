#include <iostream>
#include <math.h>
#include <float.h>
/* where the hell is this in Cygwin!!!!!
   #include <fpu_control.h>
*/


using namespace std;

union mydouble {double x; long long unsigned int i;} z;

int main()
{
  int i;
  double z1;
  z1=0.0;
  //  for(i=0;i<1000000;i++)
  //z1=nextafter(z1,1000000.0);
  
  printf("%20.18le %20.18le\n",z1,nextafter(z1,DBL_MAX));
return(1);
};
