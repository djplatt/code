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
#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"


#include "../includes/upsamdefs.h"

void fatal_error(const char *str)
{
  printf(str);
  exit(1);
}


int main(int argc, char **argv)
{
  FILE *infile,*out_file;
  int i,j,num_s,index,index2,n0,Q,INDEX;
  int q,turing_starta,turing_enda;
  int_double gamma,arg_omega,arg_omega2;
  int_double *re_zs1,*re_zs2;
  bool real_p,neg_one,first=true;
  int_complex omega,omega2;
  int diff;
  int diff1,cz1,cz2,no_zeros=0;
  q_state qs;
  qs.gap=one_over_two_B;

  _fpu_rndd();

  if(argc!=7)
    {
      printf("usage low_writeD <infile> <outfile> <rate> <n0> <Q> <INDEX>\n");
      printf("Incorrect command line. Exiting.\n");
      exit(1);
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
  qs.rate=atoi(argv[3]);
  n0=atoi(argv[4]);
  //n02=atoi(argv[5]);
  //off1=atoi(argv[6]);
  //off2=atoi(argv[7]);
  Q=atoi(argv[5]);
  INDEX=atoi(argv[6]);

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
      //print_int_complex_str("Omega=",omega);
      fread(&num_s,sizeof(unsigned int),1,infile);
      //printf("num_s=%d\n",num_s);
      if(first)
	{
	  re_zs1=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
	  if(!re_zs1) 
	    fatal_error("Fatal error allocating memory for re_zs1. Exting.\n");

	  re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
	  if(!re_zs2)
	    fatal_error("Fatal error allocating memory for re_zs2. Exting.\n");
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

      if((q!=Q)||(index!=INDEX))
	continue;
      qs.q=q;qs.index=index;qs.neg_one=neg_one;
      qs.omega[0]=omega.real.left;
      qs.omega[1]=-omega.real.right;
      qs.omega[2]=omega.imag.left;
      qs.omega[3]=-omega.imag.right;
      qs.n0=n0;
      //qs.n02=n02;
      //qs.offset1=off1;
      //qs.offset2=off2;
      fwrite(&qs,sizeof(q_state),1,out_file);
      /*
      qs.index=-index;
      qs.omega[2]=omega.imag.right;
      qs.omega[3]=-omega.imag.left;
      fwrite(&qs,sizeof(q_state),1,out_file);
      */
      fclose(out_file);
      fclose(infile);
      exit(0);
    }
}
