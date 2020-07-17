#include "stdlib.h"
#include "stdio.h"
#include "inttypes.h"
#include "math.h"
#include "../includes/r_s.h"

//#define DEBUG


long double read_del_t (FILE *ifile)
{
  long double res,tmp;
  uint64_t i;
  uint32_t j;
  uint8_t k;
  if(fread(&i,sizeof(uint64_t),1,ifile)!=1)
    {
      printf("Fatal error reading 64 bits of zero. Exiting.\n");
      exit(0);
    }
  res=i;
  //printf("res=%20.18e\n",res);
  if(fread(&j,sizeof(uint32_t),1,ifile)!=1)
    {
      printf("Fatal error reading 32 bits of zero. Exiting.\n");
      exit(0);
    }
  res+=ldexpl(j,64);
  //printf("res=%20.18e\n",res);
  if(fread(&k,sizeof(uint8_t),1,ifile)!=1)
    {
      printf("Fatal error reading 8 bits of zero. Exiting.\n");
      exit(0);
    }
  res+=ldexpl(k,96);
  //printf("res=%20.18e\n",res);
  return(ldexpl(res,-101));
}

long unsigned int first_zero;
long unsigned int last_zero;
long unsigned int zero_count=0;
long unsigned int err_count=0;
double worst_t,worst_Z=0.0;
#define MAX_ERRS (100)

void read_zeros(char *ifname, FILE *infile, long unsigned int num_zeros, double start, double end)
{
  long unsigned int z;
  long double t=start;
  long double del_t,Zt,Ztabs;
  for(z=0;z<num_zeros;z++)
    {
      del_t=read_del_t(infile);
      t+=del_t;
      Zt=Z(t);
      Ztabs=fabsl(Zt);
      if(Ztabs>worst_Z)
	{
	  worst_Z=Ztabs;
	  worst_t=t;
	}
    }
  if(t>=end)
    {
      printf("Error, zeros in block passed beyond end point. Exiting.\n");
      exit(0);
    }
  zero_count+=num_zeros;
}

double check_file (char *ifname, FILE *infile)
{
  long unsigned int num_blocks,n,zeros[2];
  double next_start,st[2];
  printf("Checking file %s\n",ifname);
  if(fread(&num_blocks,sizeof(long int),1,infile)!=1)
    {
      printf("Error reading num blocks. Exiting.\n");
      exit(0);
    }
#ifdef DEBUG
  printf("file %s: num blocks=%ld\n",ifname,num_blocks);
#endif
  for(n=0;n<num_blocks;n++)
    {
      if(fread(st,sizeof(double),2,infile)!=2)
	{
	  printf("Error reading start and end. Exiting.\n");
	  exit(0);
	}
#ifdef DEBUG
      printf("File %s: st[0]=%20.18e st[1]=%20.18e\n",ifname,st[0],st[1]);
#endif
      if(st[0]==0) // there is a segment missing so drill down
	{
	  char fname[100];
	  long unsigned int num_blocks,n1;
	  double st[2];
	  sprintf(fname,"zeros_%ld.dat",(long int) next_start);
	  FILE *newfile=fopen(fname,"rb");
	  if(!newfile)
	    {
	      printf("Failed to open file %s\n",fname);
	      exit(0);
	    }
	  printf("Checking file %s\n",fname);
	  if(fread(&num_blocks,sizeof(long int),1,newfile)!=1)
	    {
	      printf("Error reading num blocks. Exiting.\n");
	      exit(0);
	    }
#ifdef DEBUG
	  printf("file %s: num blocks=%ld\n",fname,num_blocks);
#endif
	  for(n1=0;n1<num_blocks;n1++)
	    {
	      if(fread(st,sizeof(double),2,newfile)!=2)
		{
		  printf("Error reading start and end. Exiting.\n");
		  exit(0);
		}
#ifdef DEBUG
	      printf("File %s: st[0]=%20.18e st[1]=%20.18e\n",fname,st[0],st[1]);
#endif

	      if(n!=0)
		{
		  if(st[0]!=next_start)
		    printf("Mismatch at start, expected %20.18e got %20.18e\n",next_start,st[0]);
		}
	      if(fread(zeros,sizeof(long unsigned int),2,newfile)!=2)
		{
		  printf("Error reading zero start and end. Exiting.\n");
		  exit(0);
		}
	      if(n==0)
		first_zero=zeros[0];
	      else
		{
		  if(last_zero!=zeros[0])
		    printf("Mismatch in zero count, expected %lu got %lu\n",last_zero,zeros[0]);
		}
	      last_zero=zeros[1];
	      read_zeros(fname,newfile,zeros[1]-zeros[0],st[0],st[1]);
	      next_start=st[1];
	    }
	  if(fread(zeros,sizeof(unsigned long int),1,infile)!=1)
	    {
	      printf("Error reading dead zero count after empty segment. Exiting.\n");
	      exit(0);
	    }
	  fclose(newfile);
	}
      else
	{
	  if(n!=0)
	    {
	      if(next_start!=st[0])
		printf("Mismatch in starting points.Expected %12.10e got %12.10e\n",next_start,st[0]);
	    }
	  next_start=st[1];
	  if(fread(zeros,sizeof(long int),2,infile)!=2)
	    {
	      printf("Error reading zero counts. Exiting.\n");
	      exit(0);
	    }
	  if(n==0)
	    first_zero=zeros[0];
	  else
	    {
	      if(last_zero!=zeros[0])
		printf("expecting to start at zero %ld, starting at %ld\n",last_zero,zeros[0]);
	    }
	  last_zero=zeros[1];
	  read_zeros(ifname,infile,zeros[1]-zeros[0],st[0],st[1]);
	}
    }
  //printf("last zero=%lu\n",last_zero);
  return(next_start);
}

int main(int argc, char **argv)
{
  double ans;
  char fname[100],*iname;
  FILE *infile;
  long int i;
  init_r_s();
  if(argc!=2)
    {
      printf("Usage:- read_zeros <filename>\n");
      exit(0);
    }
  for(i=0;argv[1][i]!=0;i++)
    fname[i]=argv[1][i];
  fname[i]=0;

  infile=fopen(fname,"rb");
  if(!infile)
    {
      printf("failed to open first file %s. exiting.\n",fname);
      exit(0);
    }
  ans=check_file(fname,infile);
  printf("check_file returned %12.0f\nfirst zero was %lu\nlast zero was %lu\nzero count was %lu\n",ans,first_zero,last_zero,zero_count);
  if(zero_count!=(last_zero-first_zero))
    printf("Error in zero count.\n");
  printf("Worst t was %20.18e where |Z(t)|=%20.18e\n",worst_t,worst_Z);
  fclose(infile);
  return(0);
}
