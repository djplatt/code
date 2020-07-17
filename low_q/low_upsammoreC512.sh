//
// low_upsammoreC.cpp
//
// use more detailed upsampling to find zeros missed first time
// reads output from low_upsam.cpp
// version 1.0
// 23rd Feb 2010
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
#include <unistd.h>
#include "fftw3.h"
using namespace std;

#include "../includes/int_double11.0.h"
#include "../includes/im_s.h"

#define UP_SAMPLE_RATE (512)

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
   bool real_p,neg_one,neg_one2,first_time=true;
   int_complex omega,omega2;
   int diff,t_sam_start;
   int diff1,cz1,cz2,no_zeros=0;

   _fpu_rndd();

   if(argc!=3)
     {
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

   printf("Running upsammoreC, upsampling at x %d Upsample Error=%e.\n",UP_SAMPLE_RATE,d_inter_err);

   while(fread(&q,sizeof(unsigned int),1,infile))
     {
       if(q>MAX_Q)
	 {
	   printf("q=%d exceeds MAX_Q. Exiting.\n",q);
	   exit(1);
	 }

       fread(&index,sizeof(unsigned int),1,infile);
       fread(&real_p,sizeof(bool),1,infile);
       fread(&neg_one,sizeof(bool),1,infile);
       fread(&omega,sizeof(int_complex),1,infile);
       arg_omega=argument(omega);
       fread(&num_s,sizeof(unsigned int),1,infile);
       if(first_time)
	 {
	   first_time=false;

	   im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*num_s,16);
	   if(!im_s_vec)
	     fatal_error("Fatal error allocating memory for im_s_vec. Exting.\n");
	   im_s_vec[0].im_s=0.0;
	   for(int n=1;n<num_s;n++)
	     im_s_vec[n].im_s=im_s_vec[n-1].im_s+one_over_two_B;
   
	   re_zs1=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
	   if(!re_zs1) 
	     fatal_error("Fatal error allocating memory for re_zs1. Exting.\n");
  
	   re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
	   if(!re_zs2)
	     fatal_error("Fatal error allocating memory for re_zs2. Exting.\n");

	   rig_setup(im_s_vec);
	 }
       fread(re_zs1,sizeof(int_double),num_s,infile);
       //for(int n=0;n<num_s;n+=100)
       // {
       //   print_int_double_str("",re_zs1[n]);
       // }
       //exit(0);
       fread(&turing_starta,sizeof(unsigned int),1,infile);
       fread(&turing_enda,sizeof(unsigned int),1,infile);
       fread(&gamma,sizeof(int_double),1,infile);
       if(real_p)
	 {
	   diff=rig_calc_zeros_re(turing_starta,turing_enda,q,gamma,arg_omega,re_zs1,im_s_vec,index);
	   if(diff==BAD_DIFF) // Error in Turing Zone
	     {
	       printf("q:%d index:%d (real) error in Turing Zone\n",q,index);
	       save_state1(q,index,true,neg_one,omega,omega2,num_s,re_zs1,re_zs1,turing_starta,turing_enda,gamma,out_file);
	       continue;
	     }
	   printf("Looking for %d zeros.\n",diff);
	   cz1=rig_num_zeros(0,turing_starta,re_zs1,re_zs1);
	   printf("Found %d zeros.\n",cz1);
	   no_zeros+=diff;
	   if(diff!=cz1) // missed some zeros
	     {
	       save_state1(q,index,true,neg_one,omega,omega2,num_s,re_zs1,re_zs1,turing_starta,turing_enda,gamma,out_file);
	       printf("Do some more checking on q:%d, index:%d, (real) difference:%d\n",q,index,diff-cz1);
	     }
	   else
	     printf("All zeros found on q:%d, index %d, (real).\n",q,index);
	 }
       else
	 {
	   fread(&omega2,sizeof(int_complex),1,infile);
	   fread(re_zs2,sizeof(int_double),num_s,infile);
	   arg_omega2=argument(omega2);
	   if(!(contains_zero(arg_omega+arg_omega2)))
	     {
	       printf("Omegas are not complex conjugates.\n");
	       print_int_double_str("Arg(omega1)=",arg_omega);
	       print_int_double_str("Arg(omega2)=",arg_omega2);
	       omega2=conj(omega);
	       arg_omega2=argument(omega2);
	     }

	   diff=rig_calc_zeros_cmplx(turing_starta,turing_enda,q,gamma,re_zs1,re_zs2,im_s_vec,index);
	   printf("Looking for %d zeros.\n",diff);
	   if(diff==BAD_DIFF) // Error in Turing Zone
	     {
	       printf("q:%d index:%d (cmplx) error in Turing Zone\n",q,index);
	       save_state1(q,index,false,neg_one,omega,omega2,num_s,re_zs1,re_zs2,turing_starta,turing_enda,gamma,out_file);
	       continue;
	     }
	   cz1=rig_num_zeros(0,turing_starta,re_zs1,re_zs2);
	   cz2=rig_num_zeros(0,turing_starta,re_zs2,re_zs1);
	   printf("Found %d+%d=%d zeros.\n",cz1,cz2,cz1+cz2);
	   no_zeros+=diff;
	   if((diff)!=(cz1+cz2)) // missed some zeros
	     {
	       save_state1(q,index,false,neg_one,omega,omega2,num_s,re_zs1,re_zs2,turing_starta,turing_enda,gamma,out_file);
	       printf("Do some more checking on q:%d, index:%d, (cmplx) difference:%d\n",q,index,diff-cz1-cz2);
	     }
	   else
	     printf("All zeros found on q:%d, index %d, (cmplx).\n",q,index);
	 }
     }
 }
