//
// merge_zeros.cpp
//
// merges output of find_zeros.cpp into a single file
// one per q

#include <iostream>
#include <cstdlib>
#include "characters.h"

#define NUM_FILES (16)
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

int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <specfile>.\n",argv[0]);
      exit(0);
    }

  FILE *specfile=fopen(argv[1],"r");

  if(!specfile)
    {
      printf("Failed to open %s for input. Exiting.\n",argv[1]);
      exit(0);
    }


  FILE *infiles[NUM_FILES];

  for(uint64_t f=0;f<NUM_FILES;f++)
    {
      char fname[1024];
      if(!fscanf(specfile,"%s\n",fname))
	{
	  printf("Error reading filename from spec file. Exiting.\n");
	  exit(0);
	}
      infiles[f]=fopen(fname,"rb");
      if(!infiles[f])
	{
	  printf("Error opening file %s for binary input. Exiting.\n",fname);
	  exit(0);
	}
    }
  uint64_t q,qc;

  if(fread(&q,sizeof(uint64_t),1,infiles[0])!=1)
    {
      printf("Error reading q from file 0. Exiting.\n");
      exit(0);
    }

  char outname[1024];
  sprintf(outname,"%lu/zeros_%lu.dat",q,q);
  FILE *outfile=fopen(outname,"wb");
  if(!outfile)
    {
      printf("Failed to open %s for output. Exiting.\n",outname);
      exit(0);
    }

  fwrite(&q,sizeof(uint64_t),1,outfile);
  for(uint64_t f=1,qc;f<NUM_FILES;f++)
    if((fread(&qc,sizeof(uint64_t),1,infiles[f])!=1)||qc!=q)
    {
      printf("Error reading q from file %lu. Exiting.\n",f);
      exit(0);
    }



  uint64_t f=NUM_FILES-1;
  while(true)
    {
      f++;
      if(f==NUM_FILES)
	f=0;
      uint64_t index;
      if(fread(&index,sizeof(uint64_t),1,infiles[f])!=1) // end of characters
	return(0);
      //printf("Doing index %lu\n",index);
      fwrite(&index,sizeof(uint64_t),1,outfile);

      uint64_t num_zeros;
      if(fread(&num_zeros,sizeof(uint64_t),1,infiles[f])!=1)
	{
	  printf("Error reading num_zeros for index %lu file %lu. Exiting.\n",index,f);
	  exit(0);
	}
      //printf("Processing %lu zeros\n",num_zeros);
      fpos_t pos;
      if(fgetpos(infiles[f],&pos)!=0)
	{
	  printf("Fatal error doing fgetpos. Exiting.\n");
	  exit(0);
	}
      uint64_t fzeros=0;
      for(;fzeros<num_zeros;fzeros++)
	{
	  uint64_t a;
	  if(fread(&a,sizeof(uint64_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",fzeros,index,f);
	      exit(0);
	    }
	  uint32_t b;
	  if(fread(&b,sizeof(uint32_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",fzeros,index,f);
	      exit(0);
	    }
	  uint8_t c;
	  if(fread(&c,sizeof(uint8_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",fzeros,index,f);
	      exit(0);
	    }
	  if((a==0)&&(b==0)&&(c==0))
	    break;
	}
      //printf("fzeros=%lu\n",fzeros);
      fwrite(&fzeros,sizeof(uint64_t),1,outfile);
      if(fsetpos(infiles[f],&pos)!=0)
	{
	  printf("Fatal error in fsetpos. Exiting.\n");
	  exit(0);
	}
      for(uint64_t z=0;z<fzeros;z++)
	{
	  uint64_t a;
	  if(fread(&a,sizeof(uint64_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",z,index,f);
	      exit(0);
	    }
	  uint32_t b;
	  if(fread(&b,sizeof(uint32_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",z,index,f);
	      exit(0);
	    }
	  uint8_t c;
	  if(fread(&c,sizeof(uint8_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",z,index,f);
	      exit(0);
	    }
	  fwrite(&a,sizeof(uint64_t),1,outfile);
	  fwrite(&b,sizeof(uint32_t),1,outfile);
	  fwrite(&c,sizeof(uint8_t),1,outfile);
	}

      uint64_t index1=InvMod(index,q);
      if(index1==index)
	{
	  //printf("This character was real.\n");
	  continue;
	}
      //printf("Processing conjugate index %lu\n",index1);
      read_null(infiles[f]);
      fwrite(&index1,sizeof(uint64_t),1,outfile);
      fzeros=num_zeros-fzeros;
      fwrite(&fzeros,sizeof(uint64_t),1,outfile);

      for(uint64_t z=0;z<fzeros;z++)
	{
	  uint64_t a;
	  if(fread(&a,sizeof(uint64_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",z,index1,f);
	      exit(0);
	    }
	  fwrite(&a,sizeof(uint64_t),1,outfile);
	  uint32_t b;
	  if(fread(&b,sizeof(uint32_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",z,index,f);
	      exit(0);
	    }
	  fwrite(&b,sizeof(uint32_t),1,outfile);
	  uint8_t c;
	  if(fread(&c,sizeof(uint8_t),1,infiles[f])!=1)
	    {
	      printf("Error reading zeros %lu for index %lu file %lu. Exiting.\n",z,index1,f);
	      exit(0);
	    }
	  fwrite(&c,sizeof(uint8_t),1,outfile);
	}
    }
  return(0);
}
