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
  printf("Checking file %s from %12.0f\n",ifname,start);
  fread(&num_blocks,sizeof(long int),1,infile);
#ifdef DEBUG
  printf("num blocks=%ld\n",num_blocks);
#endif
  for(n=0;n<num_blocks;n++)
    {
      fread(st,sizeof(double),2,infile);
#ifdef DEBUG
      printf("start=%12.0f finish=%12.0f\n",st[0],st[1]);
#endif
      if(st[0]==0)
	{
	  char fname[100];
	  sprintf(fname,"zeros_%ld.dat",(long int) next_start);
	  FILE *newfile=fopen(fname,"rb");
	  if(!newfile)
	    {
	      printf("Failed to open file %s\n",fname);
	      exit(0);
	    }
	  printf("Recursing...\n");
	  next_start=check_file(fname,newfile,next_start); // recursive call
#ifdef DEBUG
	  printf("check_file returned %ld\n",(long int) next_start);
#endif
	  fclose(newfile);
	  fread(zeros,sizeof(long int),1,infile);
	}
      else
	{
	  if(next_start!=st[0])
	    {
	      printf("Mismatch in starting points.\n");
	    }
	  
	  next_start=st[1];
	  fread(zeros,sizeof(long int),2,infile);
	  if(last_zero!=zeros[0])
	    printf("expecting to start at zero %ld, starting at %ld\n",last_zero,zeros[0]);
	  last_zero=zeros[1];
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
    }
  return(next_start);
}

int main(int argc, char **argv)
{
  long unsigned int inblock,block,inzs[2],f,en,st=1075646000,b_tot;
  char fname[100];
  unsigned char zero[13];
  double inst[2];
  long int i;
  FILE *infile,*outfile;

  sprintf(fname,"zeros_%ld.dat",(long int) st);
  printf("Opening infile %s\n",fname);
  if(!(infile=fopen(fname,"rb")))
    {
      printf("Failed to open file %s\n",fname);
      return(0);
    }
  fread(&inblock,sizeof(long unsigned int),1,infile);
  fread(inst,sizeof(double),2,infile);
  fread(inzs,sizeof(long unsigned int),2,infile);

  block=0;
  b_tot=0;
  for(f=0,st=1075646000;f<7;f++,st+=2100000)
    {
      en=st+2100000;
      sprintf(fname,"zeros_%ld.new",(long int) st);
      printf("Opening out file %s\n",fname);
      if(!(outfile=fopen(fname,"wb")))
	{
	  printf("Failed to open file %s\n",fname);
	  return(0);
	}
      printf("Writing dummy block count.\n");
      fwrite(&f,sizeof(long unsigned int),1,outfile); // will get overwritten
      while(true)
	{
	  //printf("In loop.\n");
	  if(inst[0]==(double) en) // end of this outfile
	    {
	      printf("End of outfile.\n");
	      fseek(outfile,0,SEEK_SET);
	      fwrite(&b_tot,sizeof(long unsigned int),1,outfile);
	      b_tot=0;
	      fclose(outfile);
	      break; // next f
	    }
	  //printf("Writing start/end and zero counts = %ld %ld.\n",inzs[0],inzs[1]);
	  fwrite(inst,sizeof(double),2,outfile);
	  fwrite(inzs,sizeof(long unsigned int),2,outfile);
	  for(long unsigned int z=inzs[0];z<inzs[1];z++)
	    {
	      fread(zero,sizeof(unsigned char),13,infile);
	      fwrite(zero,sizeof(unsigned char),13,outfile);
	    }
	  block++;
	  b_tot++;
	  //printf("Doing block %ld of infile %ld of outfile.\n",block,b_tot);
	  if(block==inblock) // end of this infile
	    {
	      fclose(infile);
	      sprintf(fname,"zeros_%ld.dat",(long int) inst[1]);
	      printf("Opening infile %s\n",fname);
	      if(!(infile=fopen(fname,"rb")))
		{
		  printf("Failed to open file %s\n",fname);
		  return(0);
		}

	      fread(&inblock,sizeof(long unsigned int),1,infile);
	      fread(inst,sizeof(double),2,infile);
	      fread(inzs,sizeof(long unsigned int),2,infile);
	      block=0;
	      continue;
	    }
	  else
	    {
	      fread(inst,sizeof(double),2,infile);
	      fread(inzs,sizeof(long unsigned int),2,infile);
	    }
	}

    }
  return(0);
}
