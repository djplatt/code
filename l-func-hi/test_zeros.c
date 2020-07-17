#include "inttypes.h"
#include "fmpr.h"
#include "fmprb.h"
#include "fmpcb.h"
//#include <iostream>
//#include <cstdlib>
//#include "characters.h"

#define PREC (300)

//using namespace std;

uint64_t get_u64(FILE *infile)
{
  uint64_t res;
  if(fread(&res,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading u64. Exiting.\n");
      exit(0);
    }
  return(res);
}

int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>.\n",argv[0]);
      exit(0);
    }

  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Error opening %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  uint64_t q=get_u64(infile);

  //dft_init(PREC);
  //DirichletGroup G(q);

  return(0);
}
