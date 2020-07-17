//
// zero_counts.cpp
//
// reads output of merge_zeros.cpp and just outputs the zero counts.
//


#include <iostream>
#include <cstdlib>
#include "characters.h"

using namespace std;

void read_null(FILE *infile)
{
  uint64_t a;
  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading null zero a. Exiting.\n");
      exit(0);
    }
  if(a!=0)
    {
      printf("Error reading null zero a neq 0. Exiting.\n");
      exit(0);
    }
  uint32_t b;
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Error reading null zero b. Exiting.\n");
      exit(0);
    }
  if(b!=0)
    {
      printf("Error reading null zero b neq 0. Exiting.\n");
      exit(0);
    }
  uint8_t c;
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Error reading null zero c. Exiting.\n");
      exit(0);
    }
  if(c!=0)
    {
      printf("Error reading null zero c neq 0. Exiting.\n");
      exit(0);
    }
}

void skip_zero(FILE *infile)
{
  uint8_t buff[13];
  if(fread(buff,sizeof(uint8_t),13,infile)!=13)
    {
      printf("Fatal error reading zero. Exiting.\n");
      exit(0);
    }
}

int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>.\n",argv[0]);
      exit(0);
    }

  FILE *zfile=fopen(argv[1],"rb");

  uint64_t q,index,index_count=0,pair_count=0;

  if(fread(&q,sizeof(uint64_t),1,zfile)!=1)
    {
      printf("Error reading q. Exiting.\n");
      exit(0);
    }

  printf("q=%lu\n",q);
  while(fread(&index,sizeof(uint64_t),1,zfile)==1)
    {
      index_count++;
      pair_count++;
      printf("Index: %lu ",index);
      uint64_t num_zeros;
      if(fread(&num_zeros,sizeof(uint64_t),1,zfile)!=1)
	{
	  printf("Error reading num_zeros for index %lu. Exiting.\n");
	  exit(0);
	}
      printf("%lu",num_zeros);
      for(uint64_t i=0;i<num_zeros;i++)
	skip_zero(zfile);

      uint64_t index1=InvMod(index,q);
      if(index1==index)
	{
	  //printf("This character was real.\n");
	  printf("\n");
	  continue;
	}
      if(fread(&index,sizeof(uint64_t),1,zfile)!=1)
	{
	  printf("Error reading conjugate index. Exiting.\n");
	  exit(0);
	}
      index_count++;
      if(index!=index1)
	{
	  printf("Conjugate index read does not correspond to expected. Exiting.\n");
	  exit(0);
	}
      printf(" %lu ",index);
      if(fread(&num_zeros,sizeof(uint64_t),1,zfile)!=1)
	{
	  printf("Error reading num_zeros for index %lu. Exiting.\n",index);
	  exit(0);
	}
      printf("%lu\n",num_zeros);
      for(uint64_t i=0;i<num_zeros;i++)
	skip_zero(zfile);
    }
  printf("Total indices read: %lu.\n",index_count);
  printf("Total pairs read: %lu\n",pair_count);
  return(0);
}
