//
// File: cos_sum.c
//
// Code to do Tim Trudgian's cosine sum
//
// 
//
// DJ Platt 
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
  int_double res=0.0;
  int_double log2=log(int_double(2));
  double eta=atol(argv[1]);
  eta/=(double)atol(argv[2]);
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

  for(uint64_t i=0;i<nfiles;i++)
    {
      sprintf(buff,"%s%lu.dat",argv[3],i);
      FILE *infile=fopen(buff,"rb");
      if(!infile)
	{
	  printf("Failed to open %s for binary input. Exiting.\n",buff);
	  exit(0);
	}
      for(uint64_t r=0;r<(1<<(h-log2nfiles));r++)
	{
	  int_double term;
	  fread(&term,sizeof(double),2,infile);
	  term*=eta;
	  res+=exp(term);
	}
      fclose(infile);
    }
  res/=((uint64_t) 1<<h);
  res=log(res);
  res/=(log2*h);
  log2/=eta;
  print_int_double_str("eta=",eta);
  printf("h=%lu\n",h);
  print_int_double_str("Without GRH lam> ",log2*(res+int_double(109)/154));
  print_int_double_str("With GRH    lam> ",log2*(res+0.5));

  return(0);
}
