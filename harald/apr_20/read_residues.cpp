#include "acb.h"
#include "stdio.h"
#include "stdlib.h"

#include "inttypes.h"

int main(int argc, char **argv)
{
  uint64_t i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Fatal error in %s: usage %s <prec> <residues file>. Exiting.\n",argv[0],argv[0]);
      exit(0);
    }
  int64_t prec=atol(argv[1]);
  printf("ARB working precision set to %ld\n",prec);
  FILE *infile=fopen(argv[2],"r");
  if(!infile)
    {
      printf("Fatal error in %s: failed to open zeros file %s for binary input. Exiting.\n",argv[0],argv[2]);
      exit(0);
    }
  acb_t res;
  arb_t max_mod,tmp;
  uint64_t zero_count=0;
  acb_init(res);
  arb_init(tmp);arb_init(max_mod);
  arb_set_ui(max_mod,0);
  while((arb_load_file(acb_realref(res),infile)==0)&&(arb_load_file(acb_imagref(res),infile)==0))
    {
      zero_count++;
      acb_abs(tmp,res,prec);
      arb_max(max_mod,max_mod,tmp,prec);
    }
  fclose(infile);
  printf("We processed %lu zeros.\nLargest absolute value found was ",zero_count);
  arb_printd(max_mod,20);
  printf("\n");

  acb_clear(res);
  arb_clear(tmp);
  arb_clear(max_mod);
  return 0;
}


