#include "../includes/int_double12.0.h"
#include "inttypes.h"

// sum records comprising the left and right end points of int doubles
// where the end points are given in decimal rep of bit pattern

int main(int argc, char **argv)
{

  if(argc!=2)
    {
      printf("Usage:- %s <records to sum>.\n",argv[0]);
      exit(0);
    }

  _fpu_rndd();

  FILE* zfile=fopen(argv[1],"r");
  if(!zfile)
    {
      printf("Error opening file %s for character input. Exiting.\n",argv[1]);
      exit(0);
    }

  uint64_t l,r;
  int_double sum=0.0,rec;

  while(fscanf(zfile,"%lu %lu\n",&rec.left,&rec.right)==2) // read as 64 bit unsigned int
                                                           // write as doubles
    sum+=rec;

  print_int_double_str("Total=",sum);
  fclose(zfile);
  return(0);
}
