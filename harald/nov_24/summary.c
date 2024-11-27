#include <flint/arb.h>
#include <stdio.h>
#include <stdlib.h>

#define PREC (200)
#define load_file(a,b) arb_load_file(b,a)

int main()
{
  FILE *infile=fopen("data/results.txt","r");
  arb_t sum1,sum1a,sum1b,sum2,sum2a,sum2b,sum3,sum3a,sum3b,tmp;
  arb_init(sum1);
  arb_init(sum1a);
  arb_init(sum1b);
  arb_init(sum2);
  arb_init(sum2a);
  arb_init(sum2b);
  arb_init(sum3);
  arb_init(sum3a);
  arb_init(sum3b);
  arb_init(tmp);

  int count=0;
  while(load_file(infile,tmp)==0)
    {
      count++;
      arb_add(sum1,sum1,tmp,PREC);
      if(load_file(infile,tmp)!=0)
	{
	  printf("Error reading data.\n");
	  return 0;
	}
      arb_add(sum1a,sum1a,tmp,PREC);
      if(load_file(infile,tmp)!=0)
	{
	  printf("Error reading data.\n");
	  return 0;
	}
      arb_add(sum1b,sum1b,tmp,PREC);
      if(load_file(infile,tmp)!=0)
	{
	  printf("Error reading data.\n");
	  return 0;
	}
      arb_add(sum2,sum2,tmp,PREC);
      if(load_file(infile,tmp)!=0)
	{
	  printf("Error reading data.\n");
	  return 0;
	}
      arb_add(sum2a,sum2a,tmp,PREC);
      if(load_file(infile,tmp)!=0)
	{
	  printf("Error reading data.\n");
	  return 0;
	}
      arb_add(sum2b,sum2b,tmp,PREC);
      if(load_file(infile,tmp)!=0)
	{
	  printf("Error reading data.\n");
	  return 0;
	}
      arb_add(sum3,sum3,tmp,PREC);
      if(load_file(infile,tmp)!=0)
	{
	  printf("Error reading data.\n");
	  return 0;
	}
      arb_add(sum3a,sum3a,tmp,PREC);
      if(load_file(infile,tmp)!=0)
	{
	  printf("Error reading data.\n");
	  return 0;
	}
      arb_add(sum3b,sum3b,tmp,PREC);
    }
  printf("sum 1 ");arb_printd(sum1,20);
  printf("\nsum 1a ");arb_printd(sum1a,20);
  printf("\nsum 1b ");arb_printd(sum1b,20);
  printf("\nsum 2 ");arb_printd(sum2,20);
  printf("\nsum 2a ");arb_printd(sum2a,20);
  printf("\nsum 2b ");arb_printd(sum2b,20);
  printf("\nsum 3 ");arb_printd(sum3,20);
  printf("\nsum 3a ");arb_printd(sum3a,20);
  printf("\nsum 3b ");arb_printd(sum3b,20);
  printf("\n%d records processed.\n",count);
  return 0;
}
