//
// upsammoreD.cpp
//
// use more detailed upsampling to find zeros missed by upsammoreC
// saves data suitable for processing by upsamdouble.cpp and upsamhigh.c
//
// version 1.0
// Created: 30th December 2009
//
// Last Modified: 30th December 2009

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

#define UP_SAMPLE_RATE (8192)

#include "../includes/upsample.h"



int main(int argc, char **argv)
{
  FILE *infile,*out_file;
  im_s *im_s_vec;
  int i,j,num_s,index,index2,n_zeros;
  int q,turing_start,turing_end,turing_starta,turing_enda;
  int_double gamma,gamma_a,arg_omega,arg_omega2;
  int_double *re_zs1,*re_zs2,gammaa,gamma_aa;
  int_double gamma_del,gamma_a_del;
  bool real_p,neg_one,neg_one2;
  int_complex omega,omega2;
  int diff,t_sam_start;
  int diff1,cz1,cz2,no_zeros=0;
  q_state qs;
  qs.rate=UP_SAMPLE_RATE;
  qs.gap=one_over_two_B;

  _fpu_rndd();

  if(argc!=3)
    {
      printf("usage upsammoreB <infile> <outfile>\n");
      fatal_error("Incorrect command line. Exiting.\n");
    }

  infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Failed to open %s for input. Exiting.\n",argv[1]);
      exit(1);
    }
  if(fread(&num_s,sizeof(int),1,infile)!=1)
    {
      printf("Empty data file.\n");
      exit(1);
    }
  out_file=fopen(argv[2],"wb");
  if(!out_file)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[2]);
      exit(1);
    }
  fread(&gamma,sizeof(int_double),1,infile);
  fread(&gamma_a,sizeof(int_double),1,infile);
  fread(&gamma_del,sizeof(int_double),1,infile);
  fread(&gamma_a_del,sizeof(int_double),1,infile);
  fread(&turing_start,sizeof(int),1,infile);
  fread(&turing_end,sizeof(int),1,infile);

  //im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*num_s,16);
  //if(!im_s_vec)
  //fatal_error("Fatal error allocating memory for im_s_vec. Exting.\n");

  re_zs1=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
  if(!re_zs1) 
    fatal_error("Fatal error allocating memory for re_zs1. Exting.\n");

  re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
  if(!re_zs2)
    fatal_error("Fatal error allocating memory for re_zs2. Exting.\n");

  //fread(im_s_vec,sizeof(im_s),num_s,infile);
  rig_setup(im_s_vec);

  while(fread(&q,sizeof(int),1,infile))
    {
      if(q>MAX_Q)
	{
	  printf("q=%d exceeds MAX_Q. Exiting.\n",q);
	  exit(1);
	}
      set_d_inter_err(q,(turing_end+1)*one_over_two_B);//im_s_vec[turing_end+1].im_s);
      fread(&index,sizeof(int),1,infile);
      printf("processing q=%d index %d\n",q,index);
      fread(&real_p,sizeof(bool),1,infile);
      fread(&n_zeros,sizeof(int),1,infile);  // not used
      fread(&neg_one,sizeof(bool),1,infile);
      fread(&omega,sizeof(int_complex),1,infile);
      qs.q=q;qs.index=index;qs.neg_one=neg_one;
      qs.omega[0]=omega.real.left;
      qs.omega[1]=-omega.real.right;
      qs.omega[2]=omega.imag.left;
      qs.omega[3]=-omega.imag.right;
      //print_int_complex_str("Omega=",omega);
      fread(re_zs1,sizeof(int_double),num_s,infile);
      arg_omega=argument(omega);
      if(!real_p)
	{
	  fread(&neg_one2,sizeof(bool),1,infile);
	  fread(&omega2,sizeof(int_complex),1,infile);
	  fread(re_zs2,sizeof(int_double),num_s,infile);
	  arg_omega2=argument(omega2);
	  //print_int_complex_str("Arg(omega1)=",omega);
	  //print_int_complex_str("Arg(omega2)=",omega2);
	}
      if(index!=831)
	continue;
      if(real_p) // it's a real character
	{
	  if(sign(re_zs1[turing_start])==UNK)
	    {
	      printf("Moving Turing Start at q=%d index = %d.\n",
		     q,index);
	      turing_starta=turing_start+1;
	      turing_enda=turing_end+1;
	      if(sign(re_zs1[turing_starta])==UNK)
		{
		printf("Two unknowns in succession at turing_start (real). Exiting.\n");
		exit(1);
		}
	      gamma_aa=gamma_a+gamma_a_del; // diff between int T1..T2
	      gammaa=gamma+gamma_del;     // and T1+5/64..T2+5/64
	    }
	  else
	    {
	      turing_starta=turing_start;
	      turing_enda=turing_end;
	      gamma_aa=gamma_a;
	      gammaa=gamma;
	    }
	  // how many should there be?
	  diff1=rig_calc_zeros_re1(turing_starta,turing_enda,q,(neg_one ? gamma_aa : gammaa),
			     arg_omega,re_zs1,im_s_vec,index,out_file,qs);
	  if(diff1==BAD_DIFF) // Error in Turing Zone
	    {
	      continue;
	    }
	  // how many can we find?
	  cz1=rig_num_zeros_D(0,turing_starta,re_zs1,re_zs1,out_file,qs);
	  no_zeros+=diff1;
	  if(diff1!=cz1)
	    printf("Do some more checking on q:%d, index:%d, (real) difference:%d\n",q,index,diff1-cz1);
	  else
	    printf("All zeros accounted for on q:%d, index:%d, (real)\n",q,index);
	}
      else
	{
	  if((sign(re_zs1[turing_start])==UNK)||(sign(re_zs2[turing_start])==UNK))
	    {
	      printf("Moving Turing Start at q=%d index = %d.\n",
		     q,index);
	      turing_starta=turing_start+1;
	      turing_enda=turing_end+1;
	      if((sign(re_zs1[turing_starta])==UNK)||(sign(re_zs2[turing_starta])==UNK))
		{
		  printf("Two unknowns in succession at turing_start (cmplx). Exiting.\n");
		  exit(1);
		}
	      gamma_aa=gamma_a+gamma_a_del; // diff between int T1..T2
	      gammaa=gamma+gamma_del;     // and T1+5/64..T2+5/64
	    }
	  else
	    {
	      turing_starta=turing_start;
	      turing_enda=turing_end;
	      gamma_aa=gamma_a;
	      gammaa=gamma;
	    }

	  diff1=rig_calc_zeros_cmplx1(turing_starta,turing_enda,q,(neg_one ? gamma_aa : gammaa),
				      re_zs1,re_zs2,im_s_vec,omega2,index,out_file,qs);
	  // omegas and index now relate to second character
	  if(diff1==BAD_DIFF) // problem in Turing method
	    {
	      continue;
	    }
	  cz2=rig_num_zeros_D(0,turing_starta,re_zs2,re_zs1,out_file,qs);//,im_s_vec,re_zs1);
	  // swap omega and index back to first character
	  qs_omega(qs,omega);
	  cz1=rig_num_zeros_D(0,turing_starta,re_zs1,re_zs2,out_file,qs);//,im_s_vec,re_zs2);
	  if(diff1!=(cz1+cz2))
	    printf("Do some more checking on q:%d, index:%d, (cmplx) difference:%d\n",q,index,diff1-cz1-cz2);
	  else
	    printf("All zeros accounted for on q:%d, index:%d, (cmplx)\n",q,index);
	}
    }
}
