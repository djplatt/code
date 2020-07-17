// **************************************************************
// check_proth.cpp
//
// DJ Platt 2013
//
// Take file written by proth1.4 and check that the
// Proth numbers stored are indeed prime and that they
// form a ladder with rungs close enough together
//
// NB This uses the CLN library built without GMP
//    and thus acts as a double check on GMP.
//
// 
#include <stdio.h>
#include <iostream>
#include "inttypes.h"
#include <stdlib.h>
#include <cln/integer.h>
#include <cln/integer_io.h>


using namespace cln;
using namespace std;


#define STEP_SIZE (4000000000000000000LL)

inline cl_I pow_mod(uint64_t a, cl_I& b, cl_I& m)
{
  cl_I ex=b;
  cl_I res=1,pow=a;
  cl_I_div_t div;
  while(ex>0)
    {
      if(oddp(ex))
	res=mod(res*pow,m);
      ex>>=1;
      pow=mod(pow*pow,m);
    }
  return(res);
}
  

int proth_p (uint64_t h, uint64_t n, uint64_t aa)
{
  cl_I pp=h;
  pp<<=n;
  cl_I pe=pp>>1;
  pp++;
  cl_I res=pow_mod(aa,pe,pp);
  return(res+1==pp);
}

int main(int argc, char **argv)
{

  if(argc!=2)
    {printf("Usage:- %s <infile>.\n",argv[0]);exit(0);}

  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Failed to open %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }


  uint64_t n,h0,h1,a,h;
  if(!fread(&n,sizeof(uint64_t),1,infile)||
     !fread(&h0,sizeof(uint64_t),1,infile)||
     !fread(&h1,sizeof(uint64_t),1,infile)||
     !fread(&a,sizeof(uint64_t),1,infile)||
     !fread(&h,sizeof(uint64_t),1,infile))
    {
      printf("Error reading data file. Exiting.\n");
      exit(0);
    }

  uint64_t max_delta=STEP_SIZE>>n;
  uint64_t first_delta=max_delta>>1;
  uint64_t last_h=h;
  if(h-h0>first_delta)
    {
      printf("First h chosen=%lu was too far from starting h=%lu,\n",h,h0);
      exit(0);
    }
  if(a!=0) // not a PARI prime
    if(!proth_p(h,n,a))
      {
	printf("Proth number at %lu*2^%lu+1 with a=%lu failed.\n",h,n,a);
	exit(0);
      }

       while(fread(&a,sizeof(uint64_t),1,infile)==1)
	 {
	   if(fread(&h,sizeof(uint64_t),1,infile)!=1)
	     {
	       printf("Error reading infile.\n");
	       exit(0);
	     }
	   // check the rung width
	   if(h-last_h>max_delta)
	     {
	       printf("Gap from h=%lu to h=%lu too large for n=%lu.\n",last_h,h,n);
	       exit(0);
	     }
	   last_h=h;
	   if(a!=0) // not a PARI prime, so check Proth prime
	     if(!proth_p(h,n,a))
	       {
		 printf("Proth number at %lu*2^%lu+1 with a=%lu failed.\n",h,n,a);
		 exit(0);
	       }
	 }

       if(h1-first_delta>last_h)
	 {
	   printf("Final prime at or after h=%lu was too far from end at h=%lu for n=%lu.\n",last_h,h1,n);
	   exit(0);
	 }
       fclose(infile);
       return(0);
}
