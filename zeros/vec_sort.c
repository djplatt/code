#include "stdio.h"

#define DTYPE double

#ifndef FAILURE
#define SUCCESS 1
#define FAILURE 0
#endif

int vec_sort(DTYPE *source_vec, DTYPE *target_vec, 
	     unsigned int gen, unsigned int len)
{
  
  /* does a sort of vec by primitive root */
  /* question, could I do this in place?  */

  unsigned int ptr,power;

  power=1;
  for(ptr=0;ptr<(len-2);ptr++)
  {
    target_vec[ptr]=source_vec[power];
    power=(power*gen)%len;
    /*    printf("%d ",power); */
    if(power==1)
      {
	printf("That was no primitive root!\n");
	return(FAILURE);
      };
  };
  target_vec[len-2]=source_vec[1];

  return(SUCCESS);
};

int main()
{
  
  DTYPE foo[10007];
  DTYPE bar[10007];

  if(!vec_sort(foo,bar,5,10007))
    printf("Catastrophic error.\n");

  return(SUCCESS);
};


