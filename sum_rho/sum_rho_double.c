#include "stdio.h"
#include "inttypes.h"
#include "stdlib.h"

// read our 13 byte structure representing a zero gap
// into a double
double in_bytes(FILE *infile)
{
  uint64_t a;uint32_t b;uint8_t c;
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

  double r1,r2;
  r1=(double) c/32.0 + (double) b/(32.0*256.0*256.0*256.0*256.0)+(double) (a&0xfff8000000000000)/(32.0*65536.0*65536.0*65536.0*65536.0*65536.0*65536.0); // /2^5, /2^37, /2^101
  r2=(double) (a&0x0007ffffffffffff)/(32.0*65536.0*65536.0*65536.0*65536.0*65536.0*65536.0); // /2^101

  return r1+r2;
}


int main(int argc, char ** argv)
{
  printf("Command line: ");
  int i;
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");

  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>\n",argv[0]);
      exit(0);
    }
  FILE*zfile=fopen(argv[1],"rb");
  if(!zfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }


  uint64_t num_its,it;
  fread(&num_its,sizeof(uint64_t),1,zfile);
  double st[2];
  uint64_t zs[2];
  unsigned long int ptr=0;
  double del_t,t,res=0.0;
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,zfile);
      fread(&zs[0],sizeof(uint64_t),1,zfile);
      if(st[0]==0.0)
	continue;
      t=st[0];
      fread(&zs[1],sizeof(long int),1,zfile);
      uint64_t z;
      for(z=zs[0]+1;z<=zs[1];z++) // zero number
	{
	  del_t=in_bytes(zfile);
	  if(del_t==0.0)
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  t+=del_t;
	  res+=1.0/t;
	}
    }
  printf("sum 1/rho = %20.18e\n",res);
  printf("Last zero was near %20.18e\n",t);
  return 0;
}
