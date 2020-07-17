//
// upsammoreC.cpp
//
// use more detailed upsampling to find zeros missed first time
// reads output from upsampler1.4.cpp
// version 1.0
// 28 December 2009
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

#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"

#define UP_SAMPLE_RATE (2048)

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

   _fpu_rndd();

   if(argc!=3)
     {
       fatal_error("Incorrect command line. Exiting.\n");
     }

   infile=fopen(argv[1],"r");
   if(!infile)
     {
       printf("Failed to open %s for input. Exiting.\n",argv[1]);
       exit(1);
     }
   if(fread(&num_s,sizeof(int),1,infile)!=1)
     {
       printf("Empty data file. Exiting.\n");
       exit(0);
     }
   
   out_file=fopen(argv[2],"wb");
   if(!out_file)
     {
       printf("Failed to open %s for binary output. Exiting.\n",argv[2]);
       exit(1);
     }

   //printf("Running upsammoreC, upsampling at x %d Upsample Error=%e.\n",UP_SAMPLE_RATE,d_inter_err);


   fread(&gamma,sizeof(int_double),1,infile);
   fread(&gamma_a,sizeof(int_double),1,infile);
   fread(&gamma_del,sizeof(int_double),1,infile);
   fread(&gamma_a_del,sizeof(int_double),1,infile);
   //fread(&foo,sizeof(int_double),1,infile);
   fread(&turing_start,sizeof(int),1,infile);
   fread(&turing_end,sizeof(int),1,infile);   
   fwrite(&num_s,sizeof(int),1,out_file);
   fwrite(&gamma,sizeof(int_double),1,out_file);
   fwrite(&gamma_a,sizeof(int_double),1,out_file);
   //print_int_double_str("gamma=",gamma);
   //print_int_double_str("gamma_a=",gamma_a);
   fwrite(&gamma_del,sizeof(int_double),1,out_file);
   fwrite(&gamma_a_del,sizeof(int_double),1,out_file);
   fwrite(&turing_start,sizeof(int),1,out_file);
   fwrite(&turing_end,sizeof(int),1,out_file);

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
   //fwrite(im_s_vec,sizeof(im_s),num_s,out_file);
   rig_setup(im_s_vec);

   while(fread(&q,sizeof(int),1,infile))
     {
       if(q>MAX_Q)
	 {
	   printf("q=%d exceeds MAX_Q. Exiting.\n",q);
	   exit(1);
	 }
       set_d_inter_err(q,(turing_end+1)*one_over_two_B/*im_s_vec[turing_end+1].im_s*/);
       fread(&index,sizeof(int),1,infile);
       printf("processing q=%d index %d\n",q,index);
       fread(&real_p,sizeof(bool),1,infile);
       fread(&n_zeros,sizeof(int),1,infile);  // not used
       fread(&neg_one,sizeof(bool),1,infile);
       fread(&omega,sizeof(int_complex),1,infile);
       fread(re_zs1,sizeof(int_double),num_s,infile);
       arg_omega=argument(omega);
       if(!real_p)
	 {
	   fread(&neg_one2,sizeof(bool),1,infile);
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
	 }
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
	   cz1=rig_num_zeros_C(0,turing_starta,re_zs1,re_zs1);
	   printf("Found %d zeros.\n",cz1);
	   if(diff==cz1)
	     {
	       printf("All zeros found on q:%d, index %d, (real).\n",q,index);
	       continue;
	     }
	   if(diff<cz1) // missed some zeros
	     printf("Found too many zeros on q:%d, index:%d, (real) excess:%d\n",q,index,cz1-diff);
	   else
	       printf("Do some more checking on q:%d, index:%d, (real) difference:%d\n",q,index,diff-cz1);
	   save_state(q,index,true,neg_one,omega,neg_one2,omega2,diff-cz1,num_s,re_zs1,re_zs1,out_file);
	 }
       else // it's a complex character
	 {
	   if((sign(re_zs1[turing_start])==UNK)||(sign(re_zs2[turing_start])==UNK))
	     {
	       printf("Moving Turing Start at q=%d index = %d.\n",
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
	   printf("Checking index %d\n",index);
	   cz1=rig_num_zeros_C(0,turing_starta,re_zs1,re_zs2);
	   printf("Checking its conjugate.\n");
	   cz2=rig_num_zeros_C(0,turing_starta,re_zs2,re_zs1);
	   printf("Found %d+%d=%d zeros.\n",cz1,cz2,cz1+cz2);
	   no_zeros+=diff;
	   cz1+=cz2;
	   if(diff==cz1)
	     {
	       printf("All zeros found on q:%d, index %d, (cmplx).\n",q,index);
	       continue;
	     }
	   if(diff<cz1) // missed some zeros
	     printf("Found too many zeros on q:%d, index:%d, (cmplx) excess:%d\n",q,index,cz1-diff);
	   else
	       printf("Do some more checking on q:%d, index:%d, (cmplx) difference:%d\n",q,index,diff-cz1);
	   save_state(q,index,false,neg_one,omega,neg_one2,omega2,diff-cz1,num_s,re_zs1,re_zs2,out_file);
	 }
     }
   if(!save_p) // we never saved anything so "erase" out_file
     {
       fclose(out_file);
       out_file=fopen(argv[2],"w");
       fclose(out_file);
     }
 }
