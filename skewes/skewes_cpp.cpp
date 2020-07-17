#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
using namespace std;



#include "../includes/int_double12.0.h"
#define OP_ACC (101)

double one_101,one_37;

void init_next_rho()
{
  one_101=1.0;
  for(long unsigned int i=0;i<101;i++)
    one_101/=2.0;
  one_37=1.0;
  for(long unsigned int i=0;i<37;i++)
    one_37/=2.0;
}

int_double next_rho(FILE *infile)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
  int_double res,res1;

  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  unsigned long int top_40=((unsigned long int)c<<32)+b;

  res=top_40;
  res*=one_37;
  res1=a;
  res1*=one_101;
  return(res+res1);
}


int main(int argc, char **argv)
{
  printf("Command Line:-");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>\n",argv[0]);
      exit(0);
    }

  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Failed to open file %s for binary input.\n",argv[1]);
      exit(0);
    }

  _fpu_rndd();
  init_next_rho();

  int_double omega=727951335795;
  omega/=1000000000L;
  int_double A=30610046000.0;
  int_double T=3001346000.0;
  int_double alpha=A*A/omega/4;
  int_double one_over_2_alpha=0.5/alpha;
  int_double eta=8*omega/A;

  long int num_its,it,z;
  double st[2];
  long int zs[2];
  int_double gam,t0,del_t,res=0,gam2,s,c;
  fread(&num_its,sizeof(long int),1,infile);
  //printf("Doing %ld iterations.\n",num_its);
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,infile);
      fread(&zs[0],sizeof(long int),1,infile);
      if(st[0]==0.0)
	continue;
      fread(&zs[1],sizeof(long int),1,infile);
      //printf("Processing zero %ld to %ld=%ld in total.\nt0=%f\n",zs[0]+1,zs[1],zs[1]-zs[0],st[0]);
      t0=st[0];
      gam=t0;
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  del_t=next_rho(infile);
          if(del_t.left==0.0)
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  gam+=del_t;
	  gam.left=(gam.left-gam.right)/2;gam.right=-gam.left;
	  gam2=gam*gam;
	  //printf("Zero %lu in block was at ",z);print_int_double(gam);printf("\n");

	  sin_cos(gam*omega,&s,&c);
	  res+=exp(-one_over_2_alpha*gam2)*(c+2*gam*s)/(0.25+gam2);

	}
    }
  printf("Result=");print_int_double(res);printf("\n");
  return(0);
}
