#include "flint/arb.h"
#include "stdbool.h"
#include "stdlib.h"

int main(int argc, char** argv)
{
  if(argc!=2)
    return 0;
  long int prec=atol(argv[1]);
  arb_t total,line;
  arb_init(total);arb_init(line);
  while(true)
    {
      if(arb_load_file(line,stdin)!=0) // read error
	{
	  printf("Error reading stdin. Exiting.\n");
	  arb_clear(total);arb_clear(line);
	  return 0;
	}
      if(arb_is_zero(line)) // end of data
	break;
      arb_add(total,total,line,prec);
    }
  printf("Sum = ");arb_printd(total,30);printf("\n");
  arb_clear(total);arb_clear(line);
  return 0;
}
