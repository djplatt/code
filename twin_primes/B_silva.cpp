// compute an estimate of B from e Silva's data
//
// data in the form
// x p2(x)
// .....

#include "stdio.h"
#include "inttypes.h"
#include "arb.h"

int main(int argc, char **argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <infile> <prec>.\n",argv[0]);
      return 0;
    }

  FILE* infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Failed to open file %s for reading.\n",argv[1]);
      return(0);
    }

  int64_t prec=atol(argv[2]);

  arb_t sum_lo,sum_hi,tmp1,tmp2,tmp3,last_sum;
  arb_init(sum_lo);
  arb_init(sum_hi);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  arb_init(last_sum);

  uint64_t last_x,c,last_c,this_x,this_c;

  if(fscanf(infile,"%lu %lu\n",&last_x,&last_c)!=2)
    {
      printf("Error reading first line.\n");
      return 0;
    }
  while(fscanf(infile,"%lu %lu\n",&this_x,&this_c)==2)
    {
      uint64_t dc=this_c-last_c;

      arb_set_ui(tmp1,dc*2);

      arb_div_ui(tmp2,tmp1,last_x,prec);
      arb_add(sum_hi,sum_hi,tmp2,prec);

      arb_div_ui(tmp2,tmp1,this_x,prec);
      arb_add(sum_lo,sum_lo,tmp2,prec);
      
      last_c=this_c;
      last_x=this_x;
    }

  printf("sum ");arb_printd(sum_lo,20);printf(" ");arb_printd(sum_hi,20);printf("\n");
  return 0;
}
