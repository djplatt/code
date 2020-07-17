#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "primesieve.h"

int main(int argc, char **argv)
{
  if(argc!=4)
    {
      printf("Usage:- %s <start> <end> <step>.\n",argv[0]);
      return(0);
    }
  uint64_t st=atof(argv[1]);
  uint64_t en=atof(argv[2]);
  uint64_t step=atof(argv[3]);
  while(st<en)
    {
      uint64_t en1=st+step;
      uint64_t count=primesieve_parallel_count_twins(st,en1);
      printf("%lu %lu %lu\n",st,en1,count);fflush(stdout);
      st=en1;
    }
  return 0;
}
