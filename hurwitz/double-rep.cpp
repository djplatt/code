#include "stdio.h"
#include "stdlib.h"
#include "assert.h"

int main(int argc, char **argv)
{
  assert(argc==2);
  double x=atof(argv[1]);
  printf("Argument is %60.58e\n",x);
  return(0);
}