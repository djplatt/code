#include "stdlib.h"
#include "stdio.h"

//#define DEBUG

//#define FIRST_FILE "zeros_3209246000.dat"

long int last_zero;

double check_file (char *ifname, FILE *infile, double start)
{
  long int num_blocks,n,zeros[2],z;
  double next_start=start,st[2];
  unsigned char zero[13];
  //printf("Checking file %s from %12.0f\n",ifname,start);
  fread(&num_blocks,sizeof(long int),1,infile);
  //printf("num blocks=%ld\n",num_blocks);
  for(n=0;n<num_blocks;n++)
    {
      //printf("Block number %ld\n",n);
      fread(st,sizeof(double),2,infile);
      //printf("start=%12.0f finish=%12.0f\n",st[0],st[1]);
      if(st[0]==0)
	{
	  printf("Fatal error. This file %s has not been merge_zero'd. Exiting.\n",ifname);
	  exit(0);
	}
      if(next_start!=st[0])
	{
	  printf("Mismatch in starting points.\n");
	}
      
      next_start=st[1];
      fread(zeros,sizeof(long int),2,infile);
      if(last_zero!=zeros[0])
	printf("expecting to start at zero %ld, starting at %ld\n",last_zero,zeros[0]);
      last_zero=zeros[1];
      printf("File %s Block %ld - t in [%.0f,%.0f] Zeros from %ld to %ld\n",ifname,n,st[0],st[1],zeros[0]+1,zeros[1]);
      if(fseek(infile,13*sizeof(unsigned char)*(zeros[1]-zeros[0]),SEEK_CUR)!=0)
	{
	  printf("error reading zeros in file %s\n",ifname);
	  exit(0);
	}
	  /*
	  for(z=zeros[0];z<zeros[1];z++)
	    if(fread(zero,sizeof(unsigned char),13,infile)!=13)
	      {
		printf("error reading zeros in file %s\n",ifname);
		exit(0);
	      }
	  */
    }
  return(next_start);
}

int main(int argc, char **argv)
{
  double ans=14.0;
  char fname[100],*iname;
  FILE *infile;
  long int i;
  if(argc!=4)
    {
      printf("Usage:- read_zeros <filename> <ans> <last_zero>\n");
      exit(0);
    }
  for(i=0;argv[1][i]!=0;i++)
    fname[i]=argv[1][i];
  fname[i]=0;

  infile=fopen(fname,"rb");
  ans=atof(argv[2]);
  last_zero=atol(argv[3]);
  if(!infile)
    {
      printf("failed to open first file %s. exiting.\n",fname);
      exit(0);
    }
  while(true)
    {
      ans=check_file(fname,infile,ans);
#ifdef DEBUG
      printf("check_file returned %12.0f\n",ans);
#endif
      fclose(infile);
#ifdef DEBUG
      printf("next file zeros_%ld.dat\n",(long int) ans);
#endif
      sprintf(fname,"zeros_%ld.dat",(long int) ans);
      if(!(infile=fopen(fname,"rb")))
	{
	  printf("Failed to open file %s\n",fname);
	  printf("Finished at zero %ld\n",last_zero);
	  return(0);
	}
    }
  return(0);
}
