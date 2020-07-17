//
// low_upsammoreD.cpp
//
// use more detailed upsampling to find zeros missed by low_upsammoreC
// saves data suitable for processing by upsamdouble.cpp and upsamhigh.c
//
// version 1.0
// Created: 25th February 2010
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
#include "../includes/int_double11.0.h"
#include "../includes/im_s.h"

#define UP_SAMPLE_RATE (8192)

#include "../includes/upsample.h"



int main(int argc, char **argv)
{
  FILE *infile,*out_file;
  im_s *im_s_vec;
  int i,j,num_s,index,index2;
  int q,turing_starta,turing_enda;
  int_double gamma,arg_omega,arg_omega2;
  int_double *re_zs1,*re_zs2;
  bool real_p,neg_one,first=true;
  int_complex omega,omega2;
  int diff;
  int diff1,cz1,cz2,no_zeros=0;
  q_state qs;
  qs.rate=UP_SAMPLE_RATE;
  qs.gap=one_over_two_B;

  _fpu_rndd();

  if(argc!=3)
    {
      printf("usage low_upsammoreD <infile> <outfile>\n");
      fatal_error("Incorrect command line. Exiting.\n");
    }

  infile=fopen(argv[1],"rb");
  if(!infile)
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

   printf("Running low_upsammoreD, upsampling at x %d Upsample Error=%e.\n",UP_SAMPLE_RATE,d_inter_err);




  while(fread(&q,sizeof(int),1,infile))
    {
      if(q>MAX_Q)
	{
	  printf("q=%d exceeds MAX_Q. Exiting.\n",q);
	  exit(1);
	}
      fread(&index,sizeof(int),1,infile);
      printf("processing q=%d index %d\n",q,index);
      fread(&real_p,sizeof(bool),1,infile);
      //fread(&n_zeros,sizeof(int),1,infile);  // not used
      fread(&neg_one,sizeof(bool),1,infile);
      fread(&omega,sizeof(int_complex),1,infile);
      qs.q=q;qs.index=index;qs.neg_one=neg_one;
      qs.omega[0]=omega.real.left;
      qs.omega[1]=-omega.real.right;
      qs.omega[2]=omega.imag.left;
      qs.omega[3]=-omega.imag.right;
      //print_int_complex_str("Omega=",omega);
      fread(&num_s,sizeof(unsigned int),1,infile);
      //printf("num_s=%d\n",num_s);
      if(first)
	{
	  im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*num_s,16);
	  if(!im_s_vec)
	    fatal_error("Fatal error allocating memory for im_s_vec. Exting.\n");

	  re_zs1=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
	  if(!re_zs1) 
	    fatal_error("Fatal error allocating memory for re_zs1. Exting.\n");

	  re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
	  if(!re_zs2)
	    fatal_error("Fatal error allocating memory for re_zs2. Exting.\n");

	  im_s_vec[0].im_s=0.0;
	  for(int n=1;n<num_s;n++)
	    im_s_vec[n].im_s=im_s_vec[n-1].im_s+one_over_two_B;
	  rig_setup(im_s_vec);
	  first=false;
	}
      fread(re_zs1,sizeof(int_double),num_s,infile);
      fread(&turing_starta,sizeof(unsigned int),1,infile);
      fread(&turing_enda,sizeof(unsigned int),1,infile);
      printf("TZ=[%d,%d]\n",turing_starta,turing_enda);
      fread(&gamma,sizeof(int_double),1,infile);
      //print_int_double_str("gamma=",gamma);


      arg_omega=argument(omega);
      if(!real_p)
	{
	  fread(&omega2,sizeof(int_complex),1,infile);
	  fread(re_zs2,sizeof(int_double),num_s,infile);
	  arg_omega2=argument(omega2);
	  //print_int_complex_str("Arg(omega1)=",omega);
	  //print_int_complex_str("Arg(omega2)=",omega2);
	}
      if(real_p) // it's a real character
	{
	  diff1=rig_calc_zeros_re1(turing_starta,turing_enda,q,gamma,
			     arg_omega,re_zs1,im_s_vec,index,out_file,qs);
	  if(diff1==BAD_DIFF) // Error in Turing Zone
	    {
	      continue;
	    }
	  // how many can we find?
	  cz1=rig_num_zeros1(0,turing_starta,re_zs1,re_zs1,out_file,qs);
	  no_zeros+=diff1;
	  if(diff1!=cz1)
	    printf("Do some more checking on q:%d, index:%d, (real) difference:%d\n",q,index,diff1-cz1);
	  else
	    printf("All zeros accounted for on q:%d, index:%d, (real)\n",q,index);
	}
      else
	{
	  diff1=rig_calc_zeros_cmplx1(turing_starta,turing_enda,q,gamma,
				      re_zs1,re_zs2,im_s_vec,omega2,index,out_file,qs);
	  // omegas and index now relate to second character
	  if(diff1==BAD_DIFF) // problem in Turing method
	    {
	      continue;
	    }
	  cz2=rig_num_zeros1(0,turing_starta,re_zs2,re_zs1,out_file,qs);//,im_s_vec,re_zs1);
	  // swap omega and index back to first character
	  qs_omega(qs,omega);
	  cz1=rig_num_zeros1(0,turing_starta,re_zs1,re_zs2,out_file,qs);//,im_s_vec,re_zs2);
	  if(diff1!=(cz1+cz2))
	    printf("Do some more checking on q:%d, index:%d, (cmpx) difference:%d\n",q,index,diff1-cz1-cz2);
	  else
	    printf("All zeros accounted for on q:%d, index:%d, (cmpx)\n",q,index);
	}
    }
}
