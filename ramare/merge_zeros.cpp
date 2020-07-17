//
// merge_zeros.cpp
//
// version original
// 4th Jan 2012
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
#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"
#include "../includes/upsamdefs.h"

#define MAX_F_NAME (256)

void fatal_error(const char* str)
{
  printf("%s. Exiting.\n",str);
  exit(0);
}

void read_check( int expect,  int i, FILE **infiles)
{
	int tmp;
	if(fread(&tmp,sizeof(int),1,infiles[i]))
		if(tmp==expect)
			return;
	fatal_error("Data mismatch between input files");
}


inline int conj_j ( int j,  int num_chi,  int q)
{
	if(q&7) // q not divisible by 8
		return(num_chi-j-1);
	if(j<(num_chi>>1))
		return((num_chi>>1)-j-1);
	return(num_chi+(num_chi>>1)-j-1);
}



void save_state(int q, int index, bool real_p, bool neg_one, int_complex &omega, 
				bool neg_one2, int_complex & omega2, int num_zeros, int num_s, 
				int_double *re_zs1, int_double *re_zs2, FILE *out_file)
{
  //save_p=true;
	fwrite(&q,sizeof(int),1,out_file);
	fwrite(&index,sizeof(int),1,out_file);
	fwrite(&real_p,sizeof(bool),1,out_file);
	fwrite(&num_zeros,sizeof(int),1,out_file);
	fwrite(&neg_one,sizeof(bool),1,out_file);
	fwrite(&omega,sizeof(int_complex),1,out_file);
	//fwrite(&num_s,sizeof(int),1,out_file);
	fwrite(re_zs1,sizeof(int_double),num_s,out_file);
	if(!real_p)
	{
		fwrite(&neg_one2,sizeof(bool),1,out_file);
		fwrite(&omega2,sizeof(int_complex),1,out_file);
		fwrite(re_zs2,sizeof(int_double),num_s,out_file);
	}
}

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
      fatal_error("Usage merge_zeros <spec file> <outfile> <num procs> <this proc>. Exiting.\n");
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
  //rig_setup(im_s_vec);
  while(fread(&q,sizeof(int),1,infiles[0]))
    {
      if(q>MAX_Q)
	{
	  printf("q=%d exceeds MAX_Q. Exiting.\n",q);
	  exit(1);
	}
      //set_d_inter_err(q,(turing_end+1)*one_over_two_B/*im_s_vec[turing_end+1].im_s*/);
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
		  if((prim%num_procs)==this_proc)
		    save_state(q,index,true,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs1,out_file);
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
		  if((prim%num_procs)==this_proc)
		    save_state(q,index,false,neg_one,omega,neg_one2,omega2,
			       0,num_s,re_zs1,re_zs2,out_file);
		}
	    }
	}
    }
  printf("merging executed successfully.\n");
}
