/*
  File f_hat_merge.cpp
*/

#include "../includes/int_double12.0.h"
#include "../includes/fft_defs.h"
#include "../includes/im_s.h"
#include "f_defs.h"
#include "stdio.h"
#include "assert.h"

//#define PRINT
//#define PRINT1

int NUM_FILES;

FILE *infile1,*infile2,*outfile;

unsigned int check_uint ()
{
  unsigned int i,v; 
  assert(fread(&v,sizeof(unsigned int),1,infile1));
  assert(fread(&i,sizeof(unsigned int),1,infile2));
  assert(v==i);
  assert(fwrite(&i,sizeof(unsigned int),1,outfile));
  return(i);
}

long unsigned int check_luint ()
{
  long unsigned int i,v; 
  assert(fread(&v,sizeof(long unsigned int),1,infile1));
  assert(fread(&i,sizeof(long unsigned int),1,infile2));
  assert(v==i);
  assert(fwrite(&i,sizeof(long unsigned int),1,outfile));
  return(i);
}

bool check_bool ()
{
  bool i,v; 
  assert(fread(&v,sizeof(bool),1,infile1));
  assert(fread(&i,sizeof(bool),1,infile2));
  assert(v==i);
  assert(fwrite(&i,sizeof(bool),1,outfile));
  return(i);
}

int_complex check_int_complex () 
{
  int_complex i,v; 
  assert(fread(&v,sizeof(int_complex),1,infile1));
  assert(fread(&i,sizeof(int_complex),1,infile2));
  assert(contains_zero(v-i));
  assert(fwrite(&i,sizeof(int_complex),1,outfile));
  return(i);
}

double check_double() 
{
  double i,v; 
  assert(fread(&v,sizeof(double),1,infile1));
  assert(fread(&i,sizeof(double),1,infile2));
  assert(v==i);
  assert(fwrite(&i,sizeof(double),1,outfile));
  return(i);
}

void print_usage()
{
  printf("Usage:- f_hat_merge in_file1 in_file2 out_file even_p \n");
  printf("        even_p = 0 => odd characters\n");
  printf("         otherwise => even characters\n"); 
  printf("Exiting.\n");
  exit(1);
}

int main(int argc, char **argv)
{
  unsigned int i,num_s1,num_s2,num_prims;
  bool real_p,neg_one,even_p;
  int_complex z;
  _fpu_rndd();
  //printf("argc=%d\n",argc);
  if(argc!=5)
    print_usage();
  if(!(infile1=fopen(argv[1],"rb")))
    {
      printf("Failed to open %s for input. Exiting.\n",argv[1]);
      exit(1);
    }
  if(!(infile2=fopen(argv[2],"rb")))
    {
      printf("Failed to open %s for input. Exiting.\n",argv[2]);
      exit(1);
    }
  if(!(outfile=fopen(argv[3],"wb")))
    {
      printf("Failed to open %s for output. Exiting.\n",argv[3]);
      exit(1);
    }
  i=atoi(argv[4]);
  even_p=(i!=0);
  if(even_p)
    printf("Processing even characters.\n");
  else
    printf("Processing odd characters.\n");
  i=check_uint(); // q
  printf("Processing q=%d\n",i);
  num_prims=check_uint(); // num_prims
  check_luint(); // N
  check_double(); // eta
  assert(fread(&num_s1,sizeof(unsigned int),1,infile1));
  assert(fread(&num_s2,sizeof(unsigned int),1,infile2));
  i=num_s1+num_s2;
  assert(fwrite(&i,sizeof(unsigned int),1,outfile));
  for(i=0;i<num_prims;i++)
    {
      neg_one=check_bool();
      real_p=check_bool();
      if(!real_p)
	i++;
      if(neg_one==even_p)
	continue;
      check_uint(); // index
      check_int_complex(); // omega
      for(int n=0;n<num_s1;n++)
	{
	  assert(fread(&z,sizeof(int_complex),1,infile1));
	  assert(fwrite(&z,sizeof(int_complex),1,outfile));
	}
      for(int n=0;n<num_s2;n++)
	{
	  assert(fread(&z,sizeof(int_complex),1,infile2));
	  assert(fwrite(&z,sizeof(int_complex),1,outfile));
	}
      if(real_p)
	continue;
      check_uint(); // index
      check_int_complex(); // omega
      for(int n=0;n<num_s1;n++)
	{
	  assert(fread(&z,sizeof(int_complex),1,infile1));
	  assert(fwrite(&z,sizeof(int_complex),1,outfile));
	}
      for(int n=0;n<num_s2;n++)
	{
	  assert(fread(&z,sizeof(int_complex),1,infile2));
	  assert(fwrite(&z,sizeof(int_complex),1,outfile));
	}
    }

  return(0);
}

