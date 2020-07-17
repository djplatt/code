#include "stdio.h"
#include "stdlib.h"

#define fname "all_primes.txt"

int main()
{
  FILE *infile;

  if(!(infile=fopen(fname,"r")))
    {
      printf("Error opening file %s for input. Exiting.\n",fname);
      exit(0);
    }
  unsigned long int res=0,i;
  while(1==1)
    {
      fscanf(infile,"%lu\n",&i);
      if(i==0)
	break;
      res+=i;
    }
  printf("Total = %lu\n",res);
}
