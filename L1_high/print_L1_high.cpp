//
// print out data in L1s_<q>.dat file
// produced by L1_high.cpp
//
// prints out as doubles although data held to +/-2^{-102}
//
#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "gmp.h"

double two_64,two_32;

// datum consists of:
//       double got from mpfi_get_d(thingy) doesn't handle 0 very well
//       14 byte representation of thingy*2^{101}
//          64 bits least sig
//          32 bits next least
//          16 bits most sig
//       if thingy is -ve then 16 bit quatity is negative
//

double get_double(char* name, FILE* infile, char* fname)
{
  double rd;
  if(fread(&rd,sizeof(double),1,infile)!=1)
    {
      printf("Error reading %s from file %s. Exiting.\n",name,fname);
      exit(0);
    }
  uint64_t ra;uint32_t rb;int16_t rc;
  if(fread(&ra,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading uint64 part of %s from %s. Exiting.\n",name,fname);
      exit(0);
    }
  if(fread(&rb,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Error reading uint32 part of %s from %s. Exiting.\n",name,fname);
      exit(0);
    }
  if(fread(&rc,sizeof(int16_t),1,infile)!=1)
    {
      printf("Error reading int16 part of %s from %s. Exiting.\n",name,fname);
      exit(0);
    }
  
  rd=ra;
  rd*=two_64;
  rd+=rb;
  rd*=two_32;
  if(rc<0)
    {
      rd+=-rc;
      rd/=32.0;
      rd=-rd;
    }
  else
    {
      rd+=rc;
      rd/=32.0;
    }
  return(rd);
}

void print_entry(char* name, FILE* infile, char* fname)
{
  double rd,id;
  rd=get_double(name,infile,fname);
  id=get_double(name,infile,fname);
  printf("%s: %18.16e %18.16e\n",name,rd,id);
}


int main(int argc, char **argv)
{
  if(argc!=2)
    {
      printf("Usage:- %s <L1_file>.\n",argv[0]);
      exit(0);
    }

  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  uint64_t q,ch;

  two_32=1.0;
  for(uint64_t i=0;i<32;i++)
    two_32/=2.0;
  two_64=two_32*two_32;

  if(fread(&q,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading modulus from %s. Exiting.\n",argv[1]);
      exit(0);
    }
  printf("Modulus %lu\n",q);
  while(fread(&ch,sizeof(uint64_t),1,infile)==1)
    {
      printf((char *)"Character %lu\n",ch);
      print_entry((char *)"Root Number",infile,argv[1]);
      print_entry((char *)"L(1/2)",infile,argv[1]);
      print_entry((char *)"L(1)",infile,argv[1]);
    }
  fclose(infile);
  return(0);
}
