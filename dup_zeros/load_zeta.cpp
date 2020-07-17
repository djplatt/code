#include "table.h"
#define ZERO_FILE "/projects/Zeros_of_Lfunctions/zzeros/zeros_14.dat"
#define DUP_FILE "zeros_dup.dat"
#define DIR_STUB "/projects/Zeros_of_Lfunctions/lzeros/zeros_%lu.dat"

int main()
{
  hash_t h;
  h.z=(uint64_t *)calloc(sizeof(uint64_t),HASH_LEN);
  h.N=(uint64_t *)calloc(sizeof(uint64_t),HASH_LEN);
  h.c=(uint64_t *)calloc(sizeof(uint64_t),HASH_LEN);
  h.t=(uint8_t *)calloc(sizeof(uint8_t),HASH_LEN);

  if(!h.z || !h.N || !h.c || !h.z)
    {
      printf("Error allocating memory. Exiting.\n");
      exit(0);
    }
  
  FILE *infile=fopen(ZERO_FILE,"r");
  if(!infile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",ZERO_FILE);
      exit(0);
    }

  uint64_t nits,zs[2];
  double st[2];
  if(fread(&nits,sizeof(long int),1,infile)!=1)
    {
      printf("Failed to read number of iterations. Exiting.\n");
      exit(0);
    }
  if(fread(st,sizeof(double),2,infile)!=2)
    {
      printf("Failed to read start and end. Exiting.\n");
      exit(0);
    }
  if(st[0]!=14.0)
    {
      printf("Expected start =14.0. Got %f. Exiting.\n",st[0]);
      exit(0);
    }
  if(fread(zs,sizeof(long int),2,infile)!=2)
    {
      printf("Failed to read zeros counts. Exiting.\n");
      exit(0);
    }
  uint64_t z=read_zero(infile);
  fclose(infile);
  z|=14ll<<60;
  if(!insert_zero(z,1,0,ZETA,&h))
    {
      printf("Failed to insert first zero for zeta. Exiting.\n");
      exit(0);
    }

  char buff[1024];

  for(uint64_t q=3;q<=10000;q++)
    {
      printf("q=%lu\n",q);
      if((q%4)==2) continue;
      sprintf(buff,DIR_STUB,q);
      FILE *infile=fopen(buff,"r");
      if(!infile)
	{
	  printf("Failed to open Dirichlet file %s. Exiting.\n",buff);
	  exit(0);
	}
      uint64_t qq;
      if(fread(&qq,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Failed to read q. Exiting.\n");
	  exit(0);
	}
      if(qq!=q)
	{
	  printf("q read %lu does not match q expected %lu. Exiting.\n",qq,q);
	  exit(0);
	}
      uint64_t c;
      while(fread(&c,sizeof(uint64_t),1,infile)==1)
	{
	  uint64_t nz;
	  if(fread(&nz,sizeof(uint64_t),1,infile)!=1)
	    {
	      printf("Failed to read number of zeros. Exiting.\n");
	      exit(0);
	    }
	  z=read_zero(infile);
	  if(q==3)
	    z|=(8LL<<60);
	  if(!insert_zero(z,q,c,DIRICHLET,&h))
	    {
	      printf("Failed to insert first Dirichlet %lu %lu. Exiting.\n",q,c);
	      exit(0);
	    }
	  fseek(infile,13*(nz-1),SEEK_CUR);
	}
      fclose(infile);
    }

  FILE *outfile=fopen(DUP_FILE,"w");
  if(!outfile)
    {
      printf("Failed to open %s for output. Exiting.\n",DUP_FILE);
      exit(0);
    }
  if(fwrite(h.z,sizeof(uint64_t),HASH_LEN,outfile)!=HASH_LEN)
    {
      printf("Failed to write zeros. Exiting.\n");
      exit(0);
    }
  if(fwrite(h.N,sizeof(uint64_t),HASH_LEN,outfile)!=HASH_LEN)
    {
      printf("Failed to write conductors. Exiting.\n");
      exit(0);
    }
  if(fwrite(h.c,sizeof(uint64_t),HASH_LEN,outfile)!=HASH_LEN)
    {
      printf("Failed to write characters. Exiting.\n");
      exit(0);
    }
  if(fwrite(h.t,sizeof(uint8_t),HASH_LEN,outfile)!=HASH_LEN)
    {
      printf("Failed to write types. Exiting.\n");
      exit(0);
    }

  return 0;
}


