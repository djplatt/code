//
// low_upsam.cpp
//
// version 1.0
//  22 February 2010
//
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <complex>
#include <limits.h>
#include "fftw3.h"
#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"

#define UP_SAMPLE_RATE (32)

#include "../includes/upsample.h"

#define one_over_A ((double) 5.0/64.0)

void print_usage()
{
  printf("Usage:- low_upsam <infile> <outfile>\n");
  exit(1);
}





int main(int argc, char **argv)
{
  FILE *in_file,*out_file;
  im_s *im_s_vec;
  unsigned int index,index2;
  unsigned int q,N0,t0,t1;
  int_double arg_omega,arg_omega2;
  int_double *re_zs1,*re_zs2;
  bool neg_one,real_p,first_q=true;
  int_complex omega,omega2;
  int_double gamma;
  unsigned int diff,t_sam_start,t_sam_end,t_sam_start_orig;
  unsigned int diff1,cz1,cz2,no_zeros=0;

  _fpu_rndd();

  //printf("Running low_upsam with error=%e\n",d_inter_err);

  if(argc!=3)
    print_usage();

  in_file=fopen(argv[1],"rb");
  if(!in_file)
    {
      printf("Failed to open %s for input. Exiting.\n",argv[1]);
      exit(1);
    }

  out_file=fopen(argv[2],"wb");
  if(!out_file)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[2]);
      exit(1);
    }

  while(fread(&q,sizeof(unsigned int),1,in_file))
    {
      //printf("q=%d\n",q);
      fread(&index,sizeof(unsigned int),1,in_file);
      //printf("index=%d\n",index);
      fread(&omega,sizeof(int_complex),1,in_file);
      //print_int_complex_str("omega=%d",omega);
      fread(&neg_one,sizeof(bool),1,in_file);
      fread(&real_p,sizeof(bool),1,in_file);
      fread(&N0,sizeof(unsigned int),1,in_file);
      //printf("N0=%d\n",N0);

      t_sam_start_orig=N0-30.0/one_over_A;
      t_sam_start=t_sam_start_orig;
      t_sam_end=t_sam_start+10.0/one_over_A;
      if(first_q)
	{
	  re_zs1=(int_double *) _aligned_malloc(sizeof(int_double)*N0,16);
	  if(!re_zs1) 
	    fatal_error("Fatal error allocating memory for re_zs1. Exting.\n");

	  re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*N0,16);
	  if(!re_zs2)
	    fatal_error("Fatal error allocating memory for re_zs2. Exting.\n");
	  //im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*N0,16);
	  //if(!im_s_vec)
	  //fatal_error("Fatal error allocating memory for im_s_vec. Exting.\n");
	  //im_s_vec[0].im_s=0.0;
	  //for(int n=1;n<N0;n++)
	  //im_s_vec[n].im_s=im_s_vec[n-1].im_s+one_over_A;

	  first_q=false;
	}
      fread(re_zs1,sizeof(int_double),N0,in_file);
      rig_setup(im_s_vec);
      arg_omega=argument(omega);
      if(real_p) // it's a real character
	{
	  while(sign(re_zs1[t_sam_start])==UNK)
	    {
	      printf("Moving Turing Start\n");
	      t_sam_start++;
	      t_sam_end++;
	    }

	  t0=t_sam_start*one_over_A;
	  t1=t_sam_end*one_over_A;
	  set_d_inter_err(q,t1);
	  if(neg_one)
	    gamma=int_odd(t0,t1);
	  else
	    gamma=int_even(t0,t1);

	  diff=rig_calc_zeros_re(t_sam_start,t_sam_end,q,gamma,arg_omega,re_zs1,im_s_vec,index);
	  if(diff==BAD_DIFF) // Error in Turing Zone
	    {
	      printf("q:%d index:%d check_zeros (real) returned:%d\n",q,index,diff);
	      save_state1(q,index,true,neg_one,omega,omega2,N0,re_zs1,re_zs1,t_sam_start,t_sam_end,gamma,out_file);
	      continue;
	    }
	  //printf("Looking for %d zeros.\n",diff);
	  diff1=rig_num_zeros_first(0,t_sam_start,re_zs1,re_zs1);

	  if(diff1>diff) // oh no!
	    {
	      save_state1(q,index,true,neg_one,omega,omega2,N0,re_zs1,re_zs1,t_sam_start,t_sam_end,gamma,out_file);
	      printf("q:%d index: %d too many zeros (%d) located.\n",q,index,diff);
	      continue;  
	    } 
	  if(diff1<diff) // missed some zeros
	    {
	      save_state1(q,index,true,neg_one,omega,omega2,N0,re_zs1,re_zs1,t_sam_start,t_sam_end,gamma,out_file);
	      printf("Do some more checking on q:%d, index:%d, (real) difference:%d\n",q,index,diff-diff1);
	      continue;
	    } 
	  printf("All zeros located for q:%d index:%d\n",q,index);
	}
      else
	{
	  fread(&index2,sizeof(int),1,in_file);
	  fread(&omega2,sizeof(int_complex),1,in_file);
	  fread(re_zs2,sizeof(int_double),N0,in_file);
	  arg_omega2=argument(omega2);
	  while((sign(re_zs1[t_sam_start])==UNK)||(sign(re_zs2[t_sam_start])==UNK))
	    {
	      t_sam_start++;
	      t_sam_end++;
	    }
	  t0=t_sam_start*one_over_A;
	  t1=t_sam_end*one_over_A;
	  set_d_inter_err(q,t1);
	  if(neg_one)
	    gamma=int_odd(t0,t1);
	  else
	    gamma=int_even(t0,t1);
	  diff=rig_calc_zeros_cmplx(t_sam_start,t_sam_end,q,gamma,re_zs1,re_zs2,im_s_vec,index);
	  if(diff==BAD_DIFF) // problem in Turing method
	    {
	      printf("q:%d index:%d check_zeros (cmplx) returned:%d\n",q,index,diff);
	      save_state1(q,index,false,neg_one,omega,omega2,N0,re_zs1,re_zs2,t_sam_start,t_sam_end,gamma,out_file);
	      continue;
	    }
	  //printf("Looking for %d zeros.\n",diff);
	  diff1=rig_num_zeros_first(0,t_sam_start,re_zs1,re_zs2)+rig_num_zeros_first(0,t_sam_start,re_zs2,re_zs1);
	   if(diff1>diff) // what th...
	     {
	       printf("q:%d index: %d %d too many zeros (%d) located.\n",q,index,index2,diff1);
	       save_state1(q,index,false,neg_one,omega,omega2,N0,re_zs1,re_zs2,t_sam_start,t_sam_end,gamma,out_file);
	       continue;
	     }
	   if(diff1<diff) // missed some zeros
	     {
	       save_state1(q,index,false,neg_one,omega,omega2,N0,re_zs1,re_zs2,t_sam_start,t_sam_end,gamma,out_file);
	       printf("Do some more checking on q:%d, index:%d, (cmpx) difference:%d\n",q,index,diff-diff1);
	       continue;
	     }
	   printf("All zeros located for q:%d index:%d %d\n",q,index,index2);
	}
    }
  printf("upsampling executed successfully.\n");
}
