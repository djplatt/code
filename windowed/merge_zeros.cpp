#include "stdlib.h"
#include "stdio.h"

//#define DEBUG

//#define FIRST_FILE "zeros_3209246000.dat"

void upd_num_blocks(long int new_val, FILE *outfile)
{
  fpos_t pos;
  fgetpos(outfile,&pos);
  fseek(outfile,0,SEEK_SET); // start of file
  fwrite(&new_val,sizeof(long int),1,outfile);
  fsetpos(outfile,&pos);
}

void merge_file (FILE *infile, FILE *outfile, double next_start)
{
  long int num_blocks,tot_blocks,zeros[2],z;
  double st[2];
  unsigned char zero[13];
  fread(&num_blocks,sizeof(long int),1,infile);
  fwrite(&num_blocks,sizeof(long int),1,outfile); // this may change
  tot_blocks=num_blocks;
#ifdef DEBUG
  printf("num blocks=%ld\n",num_blocks);
#endif
  for(long int n=0;n<num_blocks;n++)
    {
      if(fread(st,sizeof(double),2,infile)!=2)
	{
	  printf("reading start and finish failed. Exiting.\n");
	  exit(0);
	}

      //printf("start=%12.0f finish=%12.0f\n",st[0],st[1]);
      if(st[0]==0.0) // missing block(s)
	{
	  char fname[1000];
	  sprintf(fname,"../zeros9/zeros_%ld.dat",(long int) next_start);
	  printf("Recursing to file %s \n",fname);
	  FILE *newfile=fopen(fname,"rb");
	  if(!newfile)
	    {
	      printf("Failed to open file %s\n",fname);
	      exit(0);
	    }
	  //printf("Recursing...\n");
	  long int new_blocks;
	  fread(&new_blocks,sizeof(long int),1,newfile);
	  //printf("Processing %ld block in recursion.\n",new_blocks);
	  if(new_blocks!=1)
	    {
	      tot_blocks+=new_blocks-1;
	      upd_num_blocks(tot_blocks,outfile);
	    }
	  for(long int n=0;n<new_blocks;n++)
	    {
	      fread(st,sizeof(double),2,newfile);
	      next_start=st[1];
	      //printf("start=%12.0f finish=%12.0f\n",st[0],st[1]);
	      fwrite(st,sizeof(double),2,outfile);
	      fread(zeros,sizeof(long int),2,newfile);
	      fwrite(zeros,sizeof(long int),2,outfile);
	      for(long int i=zeros[0];i<zeros[1];i++)
		{
		  fread(zero,sizeof(unsigned char),13,newfile);
		  fwrite(zero,sizeof(unsigned char),13,outfile);
		}
	    }
	  fclose(newfile);
	  fread(zeros,sizeof(long int),1,infile);
	}
      else
	{
	  next_start=st[1];
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
  return;
}

int main(int argc, char **argv)
{
  FILE *infile,*outfile;
  if(argc!=4)
    {
      printf("Usage:- merge_zeros <infilename> <outfilename> <start t>\n");
      exit(0);
    }
  infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("failed to open infile %s. exiting.\n",argv[1]);
      exit(0);
    }
  outfile=fopen(argv[2],"wb");
  if(!outfile)
    {
      printf("failed to open outfile %s. exiting.\n",argv[2]);
      exit(0);
    }
  printf("Processing file %s into %s\n",argv[1],argv[2]);
  merge_file(infile,outfile,atof(argv[3]));
  return(0);
}
