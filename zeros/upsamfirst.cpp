//
// upsampler1.4.cpp
//
// version 1.0
// 3rd August 2009
// version 1.1
// 14th August 2009
// changed positioning of in up_n_twiddle
// version 1.2
// 9th September 2009
// changed gauss*sinc to be once and for all calculation
// version 1.3
// 24th September 2009
// using stationary points to locate zeros
//
// version 1.4
// 14th October 1009
// share code with upsammoreA and B via upsample.h
//
// 4th January 2009
// fixed bug in up_num_zeros
// now uses rigorous upsampling if non-rigorous suggests
// more than one sign change in an interval
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
//#include "fftw3.h"
#include "../includes/int_double13.0.h"
#include "../includes/im_s.h"

#define UP_SAMPLE_RATE (8)

#include "../includes/upsample.h"

#define MAX_F_NAME (256)



int main(int argc, char **argv)
{
  FILE *spec_file,**infiles,*out_file;
  im_s *im_s_vec;
  int num_files,num_top_files,i,j,num_s,index,index2,*num_ss;
  int q,num_chi,re_z_ptr,turing_start,turing_starta,turing_end,turing_enda,prim,prim1;
  int_double gamma,gammaa,gamma_aa,gamma_a,arg_omega,arg_omega2;
  int_double *re_zs1,*re_zs2,gamma_del,gamma_a_del;
  char fname[MAX_F_NAME];
  bool neg_one,neg_one2;
  int_complex omega,omega2;
  int diff;//,t_sam_start;
  int diff1,cz1,cz2,no_zeros=0;
  int num_procs,this_proc;

  _fpu_rndd();
  if(argc!=5)
    {
      fatal_error("Usage upsample1.4 <infile> <outfile> <num procs> <this proc>. Exiting.\n");
    }

  spec_file=fopen(argv[1],"r");
  if(!spec_file)
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
  num_procs=atoi(argv[3]);
  if((num_procs<1)||(num_procs>8)) // 8 is max number of cores available
    fatal_error("Bad num_procs specified, must be in [1,8]. Exiting.\n");
  this_proc=atoi(argv[4]);
  if((this_proc<0)||(this_proc>=num_procs))
    fatal_error("Bad this_procs specified, must be in [0,num_procs-1]. Exiting.\n");

  // first line is number of zeros files, index of start of turing zone and index of end of turing zone
  fscanf(spec_file,"%d %d %d\n",&num_files,&turing_start,&turing_end);
  //turing_end=turing_start+256;
  infiles=(FILE **) malloc(sizeof(spec_file)*num_files);
  if(!infiles)
    fatal_error("Fatal error allocating memory for infiles. Exiting.\n");
  num_ss=(int *)malloc(sizeof(int)*num_files);
  if(!num_ss)
    fatal_error("Fatal error allocating memory for num_ss. Exiting.\n");
  for(i=0;i<num_files;i++)
    {
      fscanf(spec_file,"%s\n",fname);
      infiles[i]=fopen(fname,"rb");
      if(!infiles[i])
	{
	  printf("Failed to open file |%s| for binary input. Exiting.\n",fname);
	  perror("");
	  exit(1);
	}
    }
  fclose(spec_file);
  num_s=0;
  for(i=0;i<num_files;i++)
    {
      fread(&num_ss[i],sizeof(int),1,infiles[i]);
      //printf("File No. %d num_s=%d\n",i,num_ss[i]);
      num_s+=num_ss[i];
    }
  printf("Running upsamfirst with num_s=%d, t_start=%d, t_end=%d.\n",
	 num_s,turing_start,turing_end);

  re_zs1=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
  if(!re_zs1)
    fatal_error("Fatal error allocating memory for re_zs1. Exting.\n");

  re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
  if(!re_zs2)
    fatal_error("Fatal error allocating memory for re_zs2. Exting.\n");

  index=0;
  //for(i=0;i<num_files;i++)
  //{
  //fread(&im_s_vec[index],sizeof(im_s),num_ss[i],infiles[i]);
  //index+=num_ss[i];
  //}

  gamma=int_even(turing_start*one_over_two_B,turing_end*one_over_two_B);//int_even(im_s_vec[turing_start].im_s,im_s_vec[turing_end].im_s);
  gamma_a=int_odd((turing_start+1)*one_over_two_B,(turing_end+1)*one_over_two_B);//gamma_a=int_odd(im_s_vec[turing_start].im_s,im_s_vec[turing_end].im_s);
  gamma_del=int_even(turing_start*one_over_two_B,turing_end*one_over_two_B)-gamma;//gamma_del=int_even(im_s_vec[turing_start+1].im_s,im_s_vec[turing_end+1].im_s)-gamma;
  gamma_a_del=int_odd((turing_start+1)*one_over_two_B,(turing_end+1)*one_over_two_B)-gamma_a;//gamma_a_del=int_odd(im_s_vec[turing_start+1].im_s,im_s_vec[turing_end+1].im_s)-gamma_a;
  //print_int_double_str("Int Gamma Del=",gamma_del);
  //print_int_double_str("Int Gamma Del odd=",gamma_a_del);
  fwrite(&num_s,sizeof(int),1,out_file);
  fwrite(&gamma,sizeof(int_double),1,out_file);
  fwrite(&gamma_a,sizeof(int_double),1,out_file);
  fwrite(&gamma_del,sizeof(int_double),1,out_file);
  fwrite(&gamma_a_del,sizeof(int_double),1,out_file);
  fwrite(&turing_start,sizeof(int),1,out_file);
  fwrite(&turing_end,sizeof(int),1,out_file);
  //fwrite(im_s_vec,sizeof(im_s),num_s,out_file);
  rig_setup(im_s_vec);
  while(fread(&q,sizeof(int),1,infiles[0]))
    {
      if(q>MAX_Q)
	{
	  printf("q=%d exceeds MAX_Q. Exiting.\n",q);
	  exit(1);
	}
      set_d_inter_err(q,(turing_end+1)*one_over_two_B/*im_s_vec[turing_end+1].im_s*/);
      fread(&num_chi,sizeof(int),1,infiles[0]);
      printf("processing q=%d\n",q);
      for(i=1;i<num_files;i++)
	{
	  read_check(q,i,infiles);
	  read_check(num_chi,i,infiles);
	}
      for(prim=0;prim<num_chi;prim++)
	{
	  prim1=conj_j(prim,num_chi,q);
	  if(prim<=prim1)
	    {
	      re_z_ptr=0;
	      for(i=0;i<num_files;i++)
		{
		  fread(&index,sizeof(int),1,infiles[i]);
		  fread(&omega,sizeof(int_complex),1,infiles[i]);
		  fread(&neg_one,sizeof(bool),1,infiles[i]);
		  fread(&re_zs1[re_z_ptr],sizeof(int_double),num_ss[i],infiles[i]);
		  re_z_ptr+=num_ss[i];
		}
	      arg_omega=argument(omega);
	      if(prim==prim1) // it's a real character
		{
		  if((prim%num_procs)!=this_proc)
		    continue;
		  if(sign(re_zs1[turing_start])==UNK)
		    {
		      printf("Moving Turing Start at q=%d index = %d.\n",
			     q,index);
		      turing_starta=turing_start+1;
		      turing_enda=turing_end+1;
		      if(sign(re_zs1[turing_starta])==UNK)
			{
			  printf("Two unknowns in succession at turing_start.\n");
			  save_state(q,index,true,neg_one,omega,neg_one2,omega2,
				     0,num_s,re_zs1,re_zs1,out_file);
			  continue;
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
		  diff=rig_calc_zeros_re(turing_starta,turing_enda,q,(neg_one ? gamma_aa: gammaa),
					 arg_omega,re_zs1,im_s_vec,index);
		  if(diff==BAD_DIFF) // Error in Turing Zone
		    {
		      printf("q:%d index:%d (real) error in Turing Zone\n",q,index);
		      save_state(q,index,true,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs1,out_file);
		      continue;
		    }
		  printf("Looking for %d zeros.\n",diff);
		  cz1=rig_num_zeros_first(0,turing_starta,re_zs1,re_zs1);
		  printf("Found %d zeros.\n",cz1);
		  if(diff==cz1)
		    {
		      printf("All zeros found on q:%d, index %d, (real).\n",q,index);
		      continue;
		    }
		  if(diff<cz1) // missed some zeros
		    printf("Found too many zeros on q:%d, index:%d, (real) excess: %d\n",q,index,cz1-diff);
		  else
		    printf("Do some more checking on q:%d, index:%d, (real) difference:%d\n",q,index,diff-cz1);
		  save_state(q,index,true,neg_one,omega,neg_one2,omega2,diff-cz1,num_s,re_zs1,re_zs1,out_file);
		}
	      else
		{
		  // it has a conjugate so read it
		  re_z_ptr=0;
		  for(i=0;i<num_files;i++)
		    {
		      fread(&index2,sizeof(int),1,infiles[i]);
		      fread(&omega2,sizeof(int_complex),1,infiles[i]);
		      fread(&neg_one2,sizeof(bool),1,infiles[i]);
		      fread(&re_zs2[re_z_ptr],sizeof(int_double),num_ss[i],infiles[i]);
		      re_z_ptr+=num_ss[i];
		    }
		  if((prim%num_procs)!=this_proc)
		    continue;
		  arg_omega2=argument(omega2);
		  if((sign(re_zs1[turing_start])==UNK)||(sign(re_zs2[turing_start])==UNK))
		    {
		      printf("Moving Turing Start at q=%d index =%d.\n",
			     q,index);
		      turing_starta=turing_start+1;
		      turing_enda=turing_end+1;
		      if((sign(re_zs1[turing_starta])==UNK)||(sign(re_zs2[turing_starta])==UNK))
			{
			  printf("Two unknowns in succession at turing_start.\n");
			  save_state(q,index,false,neg_one,omega,neg_one2,omega2,
				     0,num_s,re_zs1,re_zs2,out_file);
			  continue;
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
		  diff=rig_calc_zeros_cmplx(turing_starta,turing_enda,q,(neg_one ? gamma_aa: gammaa),
					    re_zs1,re_zs2,im_s_vec,index);
		  printf("Looking for %d zeros.\n",diff);
		  if(diff==BAD_DIFF) // Error in Turing Zone
		    {
		      printf("q:%d index:%d (cmplx) error in Turing Zone\n",q,index);
		      save_state(q,index,false,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs2,out_file);
		      continue;
		    }
		  cz1=rig_num_zeros_first(0,turing_starta,re_zs1,re_zs2);
		  cz2=rig_num_zeros_first(0,turing_starta,re_zs2,re_zs1);
		  printf("Found %d+%d=%d zeros.\n",cz1,cz2,cz1+cz2);
		  no_zeros+=diff;
		  cz1+=cz2;
		  if(diff==cz1)
		    {
		      printf("All zeros found on q:%d, index %d, (cmplx).\n",q,index);
		      continue;
		    }
		  if(diff<cz1) // missed some zeros
		    printf("Found too many zeros on q:%d, index:%d, (cmplx) excess=%d\n",q,index,cz1-diff);
		  else
		    printf("Do some more checking on q:%d, index:%d, (cmplx) difference:%d\n",q,index,diff-cz1);
		  save_state(q,index,false,neg_one,omega,neg_one2,omega2,diff-cz1,num_s,re_zs1,re_zs2,out_file);
		}
	    }
	}
    }
  printf("upsampling executed successfully.\n");
}
