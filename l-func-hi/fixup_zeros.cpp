//
// fixup_zeros.cpp
//
// fixes output of find_zeros_upsam
//

#include <iostream>
#include <cstdlib>
#include "characters.h"

using namespace std;

bool sensible(uint64_t q, uint64_t num_zeros)
{
  double target=200.0/(2*M_PI)*log(q*200.0/(2*M_PI*exp(1.0)));
  uint64_t low=target*0.9;
  uint64_t high=target*1.1;
  return((num_zeros>=low)&&(num_zeros<=high));
}

void get_64(uint64_t *a, FILE *infile)
{
  if(fread(a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading zeros (64 bit). Exiting.\n");
      exit(0);
    }
}

void get_32(uint32_t *a, FILE *infile)
{
  if(fread(a,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Error reading zeros (32 bit). Exiting.\n");
      exit(0);
    }
}

void get_8(uint8_t *a, FILE *infile)
{
  if(fread(a,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Error reading zeros (8 bit). Exiting.\n");
      exit(0);
    }
}

void read_zero(uint64_t *a, uint32_t *b, uint8_t *c, FILE *infile)
{
  get_64(a,infile);
  get_32(b,infile);
  get_8(c,infile);
}

void write_zero(uint64_t *a, uint32_t *b, uint8_t *c, FILE *outfile)
{
  fwrite(a,sizeof(uint64_t),1,outfile);
  fwrite(b,sizeof(uint32_t),1,outfile);
  fwrite(c,sizeof(uint8_t),1,outfile);
}

void copy_zeros (uint64_t fzeros, FILE* infile, FILE *outfile)
{
  for(uint64_t z=0;z<fzeros;z++)
    {
      uint64_t a[1];uint32_t b[1]; uint8_t c[1];
      read_zero(a,b,c,infile);
      write_zero(a,b,c,outfile);
    }
}

void write_null(FILE *outfile)
{
  uint64_t a[1];uint32_t b[1];uint8_t c[1];
  a[0]=0;b[0]=0;c[0]=0;
  write_zero(a,b,c,outfile);
}

void read_null(FILE *infile)
{
  uint64_t a[1];uint32_t b[1];uint8_t c[1];
  read_zero(a,b,c,infile);
  if((a[0]!=0)||(b[0]!=0)||(c[0]!=0))
    {
      printf("Error reading null zero. Exiting.\n");
      exit(0);
    }
}

int main(int argc, char **argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <infile> <outfile>.\n",argv[0]);
      exit(0);
    }

  FILE *infile=fopen(argv[1],"rb");

  if(!infile)
    {
      printf("Failed to open %s for input. Exiting.\n",argv[1]);
      exit(0);
    }

  FILE *outfile=fopen(argv[2],"wb");

  if(!outfile)
    {
      printf("Failed to open %s for output. Exiting.\n",argv[2]);
      exit(0);
    }


  uint64_t q;

  if(fread(&q,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading q. Exiting.\n");
      exit(0);
    }
  printf("Processing q=%lu\n",q);
  fwrite(&q,sizeof(uint64_t),1,outfile);

  uint64_t index=0;

  while(true)
    {
      if(index==0) // need to read the index
	if(fread(&index,sizeof(uint64_t),1,infile)!=1) // end of characters
	  return(0);
      if(index>q)
	{
	  printf("Bad index %lu read. Exiting.\n",index);
	  exit(0);
	}
      printf("Doing index %lu\n",index);
      fwrite(&index,sizeof(uint64_t),1,outfile);
      uint64_t index1=InvMod(index,q);
      fpos_t pos;

      if(fgetpos(infile,&pos)!=0)
	{
	  printf("Fatal error doing fgetpos. Exiting.\n");
	  exit(0);
	}

      uint64_t num_zeros;
      if(fread(&num_zeros,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading num_zeros for index %lu. Exiting.\n",index);
	  exit(0);
	}

      if(index!=index1)
	num_zeros>>=1;
      uint64_t fzeros;
      printf("num_zeros read=%lu\n",num_zeros);
      if(sensible(q,num_zeros)) // looks like it was a zero count
	{
	  printf("num_zeros looks sensible\n");
	  if(fgetpos(infile,&pos)!=0)
	    {
	      printf("Fatal error doing fgetpos. Exiting.\n");
	      exit(0);
	    }
	  fzeros=0;
	}
      else // we just read the 64 bit part of a zero
	{
	  fzeros=1;
	  uint32_t b[1];uint8_t c[1];
	  get_32(b,infile);
	  get_8(c,infile);
	}

      //printf("Processing %lu zeros\n",num_zeros);
      if(index!=index1) // complex so will end with null zero
	{
	  printf("Complex character.\n");
	  index=0; // default behaviour is read next index from file
	  for(;;fzeros++)
	    {
	      uint64_t a[1]; uint32_t b[1]; uint8_t c[1];
	      read_zero(a,b,c,infile);
	      if((a[1]==0)&&(b[1]==0)&&(c[1]==0))
		break;
	    }
	  printf("Found %lu zeros for first index.\n",fzeros);
	  fwrite(&fzeros,sizeof(uint64_t),1,outfile);
	  if(fsetpos(infile,&pos)!=0)
	    {
	      printf("Fatal error in fsetpos. Exiting.\n");
	      exit(0);
	    }
	  copy_zeros(fzeros,infile,outfile);
	  write_null(outfile);
	  fzeros=0;
	  if(fgetpos(infile,&pos)!=0)
	    {
	      printf("Fatal error doing fgetpos. Exiting.\n");
	      exit(0);
	    }
	  for(;;fzeros++)
	    {
	      uint64_t a[1];uint32_t b[1];uint8_t c[1];
	      int res=fread(a,sizeof(uint64_t),1,infile);
	      if(res==0) // end of file
		break;
	      if(a[1]<=q) // looks like the index of the next character
		{
		  printf("Run into next index %lu\n",a[1]);
		  index=a[1];
		  break;
		}
	      get_32(b,infile);
	      get_8(c,infile);
	    }
	  fwrite(&fzeros,sizeof(uint64_t),1,outfile);
	  if(fsetpos(infile,&pos)!=0)
	    {
	      printf("Fatal error in fsetpos. Exiting.\n");
	      exit(0);
	    }
	  copy_zeros(fzeros,infile,outfile);
	}
      else // a real character
	{
	  index=0;
	  fzeros=0;
	  if(fgetpos(infile,&pos)!=0)
	    {
	      printf("Fatal error doing fgetpos. Exiting.\n");
	      exit(0);
	    }
	  for(;;fzeros++)
	    {
	      uint64_t a[1];uint32_t b[1];uint8_t c[1];
	      int res=fread(a,sizeof(uint64_t),1,infile);
	      if(res==0) // end of file
		break;
	      if(a[1]<=q) // looks like the index of the next character
		{
		  index=a[1];
		  break;
		}
	      get_32(b,infile);
	      get_8(c,infile);
	    }
	  fwrite(&fzeros,sizeof(uint64_t),1,outfile);
	  if(fsetpos(infile,&pos)!=0)
	    {
	      printf("Fatal error in fsetpos. Exiting.\n");
	      exit(0);
	    }
	  copy_zeros(fzeros,infile,outfile);
	}
    }
  return(0);
}
