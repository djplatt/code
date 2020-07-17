#include "stdlib.h"
#include "stdio.h"

//#define DEBUG

void write_file (FILE *infile, FILE *outfile, long unsigned int num_blocks)
{
  long unsigned int zeros[2];
  double st[2];
  unsigned char zero[13];
  for(long int n=0;n<num_blocks;n++)
    {
      if(fread(st,sizeof(double),2,infile)!=2)
	{
	  printf("reading start and finish failed. Exiting.\n");
	  exit(0);
	}
      fwrite(st,sizeof(double),2,outfile);
      fread(zeros,sizeof(long int),2,infile);
      fwrite(zeros,sizeof(long int),2,outfile);
      for(long int i=zeros[0];i<zeros[1];i++)
	{
	  fread(zero,sizeof(unsigned char),13,infile);
	  fwrite(zero,sizeof(unsigned char),13,outfile);
	}
    }
}

void concat_file (FILE *infile1, FILE *infile2, FILE *outfile)
{
  long int num_blocks,num_blocks1,num_blocks2;
  fread(&num_blocks1,sizeof(long int),1,infile1);
  fread(&num_blocks2,sizeof(long int),1,infile2);
  num_blocks=num_blocks1+num_blocks2;
  fwrite(&num_blocks,sizeof(long int),1,outfile);
  write_file(infile1,outfile,num_blocks1);
  write_file(infile2,outfile,num_blocks2);
}

int main(int argc, char **argv)
{
  FILE *infile1,*infile2,*outfile;
  if(argc!=4)
    {
      printf("Usage:- concat_zeros <infilename1> <infilename2> <outfilename>\n");
      exit(0);
    }
  infile1=fopen(argv[1],"rb");
  if(!infile1)
    {
      printf("failed to open infile1 %s. exiting.\n",argv[1]);
      exit(0);
    }
  infile2=fopen(argv[2],"rb");
  if(!infile2)
    {
      printf("failed to open infile2 %s. exiting.\n",argv[2]);
      exit(0);
    }
  outfile=fopen(argv[3],"wb");
  if(!outfile)
    {
      printf("failed to open outfile %s. exiting.\n",argv[3]);
      exit(0);
    }
  printf("Processing files %s and %s into %s\n",argv[1],argv[2],argv[3]);
  concat_file(infile1,infile2,outfile);
  return(0);
}
