/*
  File: merge_terms.cpp

  Created: 14th December 2010

  Version: <v> = 1.0

  1.0 Initial implementation

  Dialect: C++


  Implementation notes: 

  Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

  merges two (consecutive) f_hat_odd or even files

  By: DJ Platt
  Bristol University

  Copyright 2010.

  This work is funded by the UK ESPRC. */

#define VERSION ""
//#define PRINT
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft.h"
#include "f_defs.h"

void print_usage()
/* called when wrong arguments passed via command line */
{
  printf("Usage: merge_terms%s (ifname1) (ifname2) (ofname)\n",VERSION);
  printf("  (ifname1)     - output from F_hat_odd/even_terms.c\n");
  printf("  (ifname2)     - ditto\n");
  printf("  (ofname)      - output\n");
  exit(1);
}


void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};

int main(int argc, char **argv)
{

  FILE *infile1,*infile2,*out_file;
  unsigned int q;

  if(argc!=4)
    print_usage();

  
  infile1=fopen(argv[1],"rb");
  if(!infile1)
    fatal_error("Couldn't open infile1. Exiting.\n");
  infile2=fopen(argv[2],"rb");
  if(!infile2)
    fatal_error("Couldn't open infile2. Exiting.\n");
  out_file=fopen(argv[3],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");

  fread(&q,sizeof(unsigned int),1,infile);
  check_uint

  if(!(factors=(factor *) calloc(q+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  if(!read_factors(factors,q,facs_file))
    fatal_error("Error reading factor file. Exiting.\n");

  facs_file.close();

  fread(&eta_odd,sizeof(double),1,infile);  
  fread(&n0,sizeof(unsigned int),1,infile);  
  fread(&n1,sizeof(unsigned int),1,infile);  
  fread(&N,sizeof(unsigned int),1,infile);  

  B=N*one_over_A;
  //printf("B=%f\n",B);
  phi_q=factors[q].phi;

  num_s=n1-n0;
  //printf("num_s=%d\n",num_s);
  n_prims=num_prims(q,factors); // no of prim characters mod q
  //n_prims=8;
  if(!(prims=(unsigned int *) malloc(n_prims*sizeof(unsigned int))))
    fatal_error("Couldn't allocate memory for prims.\n");
  if(!(indices=(unsigned int *) malloc(n_prims*sizeof(unsigned int))))
    fatal_error("Couldn't allocate memory for indices.\n");
  if(!(f_omegas=(int_complex *) _aligned_malloc(n_prims*sizeof(int_complex),16)))
    fatal_error("Couldn't allocate memory for omegas.\n");
  if(!(neg_ones=(bool *) malloc(n_prims*sizeof(bool))))
    fatal_error("Couldn't allocate memory for neg_ones.\n");
  if(!(zs=(int_complex *) _aligned_malloc(((n_prims+1)/2)*num_s*sizeof(int_complex),16)))
    fatal_error("Couldn't allocate memory for zs.\n");

  //printf("num_prims=%d\n",n_prims);

  fwrite(&q,sizeof(unsigned int),1,out_file);
  fwrite(&n_prims,sizeof(unsigned int),1,out_file);
  fwrite(&N,sizeof(int),1,out_file);
  fwrite(&eta_odd,sizeof(double),1,out_file);
  fwrite(&num_s,sizeof(int),1,out_file);
  //fwrite(&n0,sizeof(int),1,out_file);

  fwrite(&q,sizeof(unsigned int),1,out_file2);
  fwrite(&n_prims,sizeof(unsigned int),1,out_file2);
  fwrite(&N,sizeof(int),1,out_file2);
  fwrite(&eta_odd,sizeof(double),1,out_file2);

  if(fnum==MIDDLE)
    num_s2=num_s;
  else
    {
      //printf("This is the first or last file.\n");
      num_s2=num_s-1;
    }
  fwrite(&num_s2,sizeof(int),1,out_file2);
  //fwrite(&n0,sizeof(int),1,out_file2);


  if(!(omegas=(int_complex *) _aligned_malloc(phi_q*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  //printf("Initialising ws...\n");
  init_ws(_ws);

  //printf("Starting...\n");
  make_l_odd(q,num_s,factors,q_n,offset_vec,a_n,omegas,B,n0,eta_odd,infile);

  for(prim=0;prim<n_prims;prim++)
    {
      prim1=conj_j(prim,n_prims,q);
      real_p=(prim1==prim);
      if(prim1<prim)
	continue;
      //printf("Saving q:%d ind1:%d ind2:%d\n",q,prim,prim1);
      fwrite(&neg_ones[prim],sizeof(bool),1,out_file);
      fwrite(&real_p,sizeof(bool),1,out_file);
      fwrite(&neg_ones[prim],sizeof(bool),1,out_file2);
      fwrite(&real_p,sizeof(bool),1,out_file2);
      if(!neg_ones[prim])
	continue;
      fwrite(&indices[prim],sizeof(unsigned int),1,out_file);
      fwrite(&f_omegas[prim],sizeof(int_complex),1,out_file);
      fwrite(&indices[prim],sizeof(unsigned int),1,out_file2);
      fwrite(&f_omegas[prim],sizeof(int_complex),1,out_file2);
      //print_int_complex_str("saving",zs[prims[prim]*num_s+1]);
      fwrite(&zs[prims[prim]*num_s],sizeof(int_complex),num_s,out_file);
      if(fnum==FIRST)
	{
	  conj_reverse(&zs[prims[prim]*num_s+1],num_s2);
	  fwrite(&zs[prims[prim]*num_s+1],sizeof(int_complex),num_s2,out_file2);
	}
      else
	{
	  conj_reverse(&zs[prims[prim]*num_s],num_s2);
	  fwrite(&zs[prims[prim]*num_s],sizeof(int_complex),num_s2,out_file2);
	}
      if(!real_p) // it has a conjugate so save it
	{
	  //printf("Saving composite ind:%d\n",prim1);
	  //printf("prims[%d]=%d\n",prim1,prims[prim1]);

	  fwrite(&indices[prim1],sizeof(unsigned int),1,out_file);
	  fwrite(&f_omegas[prim1],sizeof(int_complex),1,out_file);
	  fwrite(&indices[prim1],sizeof(unsigned int),1,out_file2);
	  fwrite(&f_omegas[prim1],sizeof(int_complex),1,out_file2);
	  //print_int_complex_str("saving",zs[prims[prim1]*num_s+1]);
	  fwrite(&zs[prims[prim1]*num_s],sizeof(int_complex),num_s,out_file);
	  if(fnum==FIRST)
	    {
	      // the zero'th element is not replicated
	      conj_reverse(&zs[prims[prim1]*num_s+1],num_s2);
	      fwrite(&zs[prims[prim1]*num_s+1],sizeof(int_complex),num_s2,out_file2);
	    }
	  else
	    {
	      // num_s2 takes care of the last element where applicable
	      conj_reverse(&zs[prims[prim1]*num_s],num_s2);
	      fwrite(&zs[prims[prim1]*num_s],sizeof(int_complex),num_s2,out_file2);
	    }
	}
      //else
      //printf("\n");
    }


  fclose(out_file);
  fclose(out_file2);
  fclose(infile);

  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  printf("f_hat_odd1.2 successful completion on q=%d n=[%d,%d)\n",q,n0,n1);
  return(0);
}
