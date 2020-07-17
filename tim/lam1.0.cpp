//
// File: cos_sum.c
//
// Code to do Tim Trudgian's cosine sum
//
// 
//
// DJ Platt 
//
// version original
// version 1.0 iterates towards minimum point
//
#include "stdlib.h"
#include "stdio.h"
#include "malloc.h"
#include "inttypes.h"
#include "../includes/int_double12.0.h"

//#define h ((uint64_t) 23)


int main(int argc, char** argv)
{
  if(argc!=6)
    {
      printf("Usage %s <etan> <etad> <infile stub> <log2 nfiles> <h>\n",argv[0]);
      exit(0);
    }

  _fpu_rndd();
  char buff[1024];
  int_double res[3]={0.0,0.0,0.0};
  int_double log2=log(int_double(2));
  double eta[3];
  uint64_t num=atol(argv[1]);
  uint64_t den=atol(argv[2]);
  int64_t log2nfiles=atol(argv[4]);
  int64_t h=atol(argv[5]);
  if((log2nfiles<=0)||(log2nfiles>10))
    {
      printf("log2 nfiles must be in [1,10]\n");
      exit(0);
    }
  uint64_t nfiles=1<<log2nfiles;

  if((h<=0)||(h>32))
    {
      printf("h must be in [1,32].\n");
      exit(0);
    }

  eta[1]=num;
  eta[0]=eta[1]-1;
  eta[2]=eta[1]+1;
  for(uint64_t i=0;i<3;i++)
    eta[i]/=den;
  for(uint64_t j=0;j<nfiles;j++)
    {
      sprintf(buff,"%s%lu.dat",argv[3],j);
      FILE *infile=fopen(buff,"rb");
      if(!infile)
	{
	  printf("Failed to open %s for binary input. Exiting.\n",buff);
	  exit(0);
	}
      for(uint64_t r=0;r<(1<<(h-log2nfiles));r++)
	{
	  int_double term,term1;
	  fread(&term,sizeof(double),2,infile);
	  for(uint64_t i=0;i<3;i++)
	    res[i]+=exp(term*eta[i]);
	}
      fclose(infile);
    }
  for(uint64_t i=0;i<3;i++)
    {
      res[i]=log(res[i]/((uint64_t) 1<<h));
      res[i]/=(log2*h);
      res[i]+=int_double(109)/154; // w/o GRH. otherwise use +=0.5
      res[i]*=log2;
      res[i]/=eta[i];
      printf("eta=%lu/%lu gives ",num+i-1,den);print_int_double_str("",res[i]);
    }
  fflush(stdout);


  while(den<((uint64_t) 1<<32))
    {
      den<<=1;
      if(res[0].left<res[1].left)
	{
	  res[1]=res[0];
	  num=(num-1)<<1;
	}
      else
	{
	  if(res[2].left<res[1].left)
	    {
	      res[1]=res[2];
	      num=(num+1)<<1;
	    }
	  else
	    num=num<<1;
	}
      for(uint64_t i=0;i<3;i++)
	{
	  eta[i]=num+i-1;
	  eta[i]/=den;
	}
      for(uint64_t j=0;j<nfiles;j++)
	{
	  sprintf(buff,"%s%lu.dat",argv[3],j);
	  FILE *infile=fopen(buff,"rb");
	  if(!infile)
	    {
	      printf("Failed to open %s for binary input. Exiting.\n",buff);
	      exit(0);
	    }
	  for(uint64_t r=0;r<(1<<(h-log2nfiles));r++)
	    {
	      int_double term,term1;
	      fread(&term,sizeof(double),2,infile);
	      res[0]+=exp(term*eta[0]);
	      res[2]+=exp(term*eta[2]);
	    }
	  fclose(infile);
	}
      for(uint64_t i=0;i<3;i+=2)
	{
	  res[i]=log(res[i]/((uint64_t) 1<<h));
	  printf("log(F(xi,h)=");print_int_double_str("",res[i]);

	  res[i]/=(log2*h);
	  res[i]+=0.5;//int_double(109)/154; // w/o GRH. otherwise use +=0.5
	  res[i]*=log2;
	  res[i]/=eta[i];
	}
      for(uint64_t i=0;i<3;i++)
	{printf("eta=%lu/%lu gives ",num+i-1,den);print_int_double_str("",res[i]);}
      fflush(stdout);
    }

  return(0);
}
