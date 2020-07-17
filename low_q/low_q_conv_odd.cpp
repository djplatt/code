// low_q_conv_odd.cpp
// convert file from output of fft to input for zeros.cpp
// MS Stuff, does no harm elsewhere
// created 16/02/2010
//
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <assert.h>
#include "../includes/int_double11.0.h"

#define VERSION ""
// 16/02/10
// 1.0 Original



void print_usage()
/* called when wrong arguments passed via command line */
{
	printf("Usage: low_q_conv_odd%s (in) (out)\n",VERSION);
	printf("  (in)  - file from low_q_odd\n");
	printf("  (out) - output file\n");
	exit(1);
}

void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
	std::cout << error_string << " Exiting." << endl;
	exit(1);
}


unsigned int conj_j (unsigned int j, unsigned int num_chi, unsigned int q)
{
	if(q&7) // q not divisible by 8
	  return(num_chi-j-1);
	if(j<(num_chi>>1))
	  return((num_chi>>1)-j-1);
	return(num_chi+(num_chi>>1)-j-1);
}


int main(int argc, char **argv)
{
  unsigned int num_s,q,n0,N,last_ind,*prims,prim_ptr;
  unsigned int s,num_prims,prim,*indices,prim1,num_odds;
  FILE *hur_file,*out_file;
  
  int_complex *omegas;
  int_complex *zs;
  bool *neg_ones,real_p;
  double eta;
  
  _fpu_rndd();
  
  clock_t no_clicks;
  
  no_clicks=clock(); // start timing
  
  if(argc!=3)
    print_usage();
  
  hur_file=fopen(argv[1],"rb");
  if(!hur_file)
    fatal_error("Couldn't open data file. Exiting.\n");
  
  out_file=fopen(argv[2],"wb");
  if(!out_file)
    fatal_error("Couldn't open output file. Exiting.\n");

  assert(fread(&q,sizeof(unsigned int),1,hur_file));
  fwrite(&q,sizeof(unsigned int),1,out_file);
  assert(fread(&num_s,sizeof(unsigned int),1,hur_file));    // number of s values in file
  //printf("num_s=%d\n",num_s);
  assert(fread(&num_prims,sizeof(unsigned int),1,hur_file)); // number of primitives we will encounter
  fwrite(&num_prims,sizeof(unsigned int),1,out_file);
  //printf("num_prims=%d\n",num_prims);
  assert(fread(&n0,sizeof(unsigned int),1,hur_file));
  printf("in low_q_conv_odd with q=%d and n0=%d\n",q,n0);
  assert(fread(&N,sizeof(unsigned int),1,hur_file));
  fwrite(&N,sizeof(unsigned int),1,out_file);
  //printf("N=%d\n",N);
  assert(fread(&eta,sizeof(double),1,hur_file));
  //printf("eta=%f\n",eta);
  fwrite(&eta,sizeof(double),1,out_file);
  fwrite(&num_s,sizeof(unsigned int),1,out_file);
  fwrite(&n0,sizeof(unsigned int),1,out_file);

 
  if(!(prims=(unsigned int *) malloc(num_prims*sizeof(unsigned int))))
    fatal_error("Couldn't allocate memory for prims.\n");
  if(!(indices=(unsigned int *) malloc(num_prims*sizeof(unsigned int))))
    fatal_error("Couldn't allocate memory for indices.\n");
  if(!(omegas=(int_complex *) _aligned_malloc(num_prims*sizeof(int_complex),16)))
    fatal_error("Couldn't allocate memory for omegas.\n");
  if(!(neg_ones=(bool *) malloc(num_prims*sizeof(bool))))
    fatal_error("Couldn't allocate memory for neg_ones.\n");
  if(!(zs=(int_complex *) _aligned_malloc(((num_prims+1)/2)*num_s*sizeof(int_complex),16)))
    fatal_error("Couldn't allocate memory for zs.\n");
  prim_ptr=0;
  for(prim=0;prim<num_prims;prim++)
    {
      assert(fread(&indices[prim],sizeof(unsigned int),1,hur_file));
      //printf("indices[%d]=%d\n",prim,indices[prim]);
      assert(fread(&neg_ones[prim],sizeof(bool),1,hur_file));
      if(!neg_ones[prim])
	continue;
      prims[prim]=prim_ptr;
      assert(fread(&omegas[prim],sizeof(int_complex),1,hur_file));
      assert(fread(&zs[prims[prim]*num_s],sizeof(int_complex),1,hur_file));
      prim_ptr++;
    }

  for(s=1;s<num_s;s++)
    for(prim=0;prim<num_prims;prim++)
      if(neg_ones[prim])
	assert(fread(&zs[prims[prim]*num_s+s],sizeof(int_complex),1,hur_file));
  /*
  for(prim=0;prim<num_prims;prim++)
    {
      if(neg_ones[prim])
	{
	  printf("character - %d is odd\n",prim);
	  printf("zs[%d][0]=",prim);
	  print_int_complex_str("",zs[prims[prim]*num_s]);
	}
      else
	printf("character - %d is even\n",prim);
    }
  */

  for(prim=0;prim<num_prims;prim++)
    {
      prim1=conj_j(prim,num_prims,q);
      real_p=(prim1==prim);
      //printf("prim=%d pr=%d prim1=%d pr1=%d\n",prim,pr,prim1,pr1);
      if(prim1<prim)
	continue;
      //printf("Saving q:%d ind:%d\n",q,prim);
      fwrite(&neg_ones[prim],sizeof(bool),1,out_file);
      fwrite(&real_p,sizeof(bool),1,out_file);
      if(!neg_ones[prim])
	continue;
      fwrite(&indices[prim],sizeof(unsigned int),1,out_file);
      fwrite(&omegas[prim],sizeof(int_complex),1,out_file);
      //print_int_complex_str("saving",zs[prims[prim]*num_s+1]);
      fwrite(&zs[prims[prim]*num_s],sizeof(int_complex),num_s,out_file);
      if(!real_p) // it has a conjugate so save it
	{
	  //printf("Saving composite ind:%d\n",prim1);
	  fwrite(&indices[prim1],sizeof(unsigned int),1,out_file);
	  fwrite(&omegas[prim1],sizeof(int_complex),1,out_file);
	  //print_int_complex_str("saving",zs[prims[prim1]*num_s+1]);
	  fwrite(&zs[prims[prim1]*num_s],sizeof(int_complex),num_s,out_file);
	}
      //else
      //printf("\n");
    }
  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  return(0);
}
