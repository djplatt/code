//
// upsammoreA.cpp
//
// use more detailed upsampling to fing zeros missed first time
// version 1.0
// 21st September 2009
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

#define UP_SAMPLE_RATE (32)

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

   if(argc!=2)
     {
       fatal_error("Incorrect command line. Exiting.\n");
     }

   infile=fopen(argv[1],"r");
   if(!infile)
     {
       printf("Failed to open %s for input. Exiting.\n",argv[1]);
       exit(1);
     }
   


   fread(&num_s,sizeof(int),1,infile);
   fread(&gamma,sizeof(int_double),1,infile);
   fread(&gamma_a,sizeof(int_double),1,infile);
   fread(&gamma_del,sizeof(int_double),1,infile);
   fread(&gamma_a_del,sizeof(int_double),1,infile);
   //fread(&foo,sizeof(int_double),1,infile);
   fread(&turing_start,sizeof(int),1,infile);
   fread(&turing_end,sizeof(int),1,infile);   


   t_sam_start=num_s-UP_SAMPLE_IN_WIDTH;

   if((num_s-turing_end)<UP_SAMPLE_SACRIFICE)
     fatal_error("End of Turing region to close to end of data. Exiting.\n");
   if((turing_start-t_sam_start)<UP_SAMPLE_SACRIFICE)
     fatal_error("Start of Turing region to close to start of up_sample. Exiting.\n");
   
   
   im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*num_s,16);
   if(!im_s_vec)
     fatal_error("Fatal error allocating memory for im_s_vec. Exting.\n");
   
   re_zs1=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
   if(!re_zs1) 
     fatal_error("Fatal error allocating memory for re_zs1. Exting.\n");
  
   re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
   if(!re_zs2)
     fatal_error("Fatal error allocating memory for re_zs2. Exting.\n");

   fread(im_s_vec,sizeof(im_s),num_s,infile);
   printf("setting up \n");
   setup(im_s_vec);
   printf("all setup\n");

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
       fread(&n_zeros,sizeof(int),1,infile);  // not used
       fread(&neg_one,sizeof(bool),1,infile);
       fread(&omega,sizeof(int_complex),1,infile);
       fread(re_zs1,sizeof(int_double),num_s,infile);
       if(re_zs1[0].right>0.0)
	 {
	   omega=-omega;
	   for(i=0;i<num_s;i++)
	     re_zs1[i]=-re_zs1[i];
	 }
       arg_omega=argument(omega);
       //for(i=0;i<=turing_start;i++)
       // printf("%10.8e\n",re_zs1[i].left);
       //exit(0);
       /*
       if(arg_omega<d_pi/2)
	 arg_omega+=d_pi;
       if(arg_omega>d_pi/2)
	 arg_omega-=d_pi;
       */
       if(!real_p)
	 {
	   fread(&neg_one2,sizeof(bool),1,infile);
	   fread(&omega2,sizeof(int_complex),1,infile);
	   fread(re_zs2,sizeof(int_double),num_s,infile);
	   if(re_zs2[0].right>0.0)
	     {
	       omega2=-omega2;
	       for(i=0;i<num_s;i++)
		 re_zs2[i]=-re_zs2[i];
	     }
	   arg_omega2=argument(omega2);
	   /*
	   if((arg_omega>d_pi/2)||(arg_omega<-d_pi/2))
	     {
	       print_int_double_str("re_zs1[0]=",re_zs1[0]);
	       print_int_double_str("re_zs2[0]=",re_zs2[0]);
	       print_int_double_str("re_zs1[turing_start]=",re_zs1[turing_start]);
	       print_int_double_str("re_zs2[turing_start]=",re_zs2[turing_start]);
	       print_int_double_str("arg(omega1)=",arg_omega);
	       print_int_double_str("arg(omega2)=",arg_omega2);
	     }
	   */
	 }
       if(real_p) // it's a real character
	 {
	   if(sign(re_zs1[turing_start])==UNK)
	     {
	       printf("Moving Turing Start at q=%d index = %d. Only works for 1250..1260\n",
		      q,index);
	       turing_starta=turing_start+1;
	       turing_enda=turing_end+1;
	       if(sign(re_zs1[turing_starta])==UNK)
		 printf("Two unknowns in succession at turing_start.\n");
	       gamma_aa=gamma_a+gamma_a_del; // diff between int 1250..1260
	       gammaa=gamma+gamma_del;     // and 1250+5/64..1260+5/64
	     }
	   else
	     {
	       turing_starta=turing_start;
	       turing_enda=turing_end;
	       gamma_aa=gamma_a;
	       gammaa=gamma;
	     }
	 }
       else
	 {
	   if((sign(re_zs1[turing_start])==UNK)||(sign(re_zs2[turing_start])==UNK))
	     {
	       printf("Moving Turing Start at q=%d index = %d. Only works for 1250..1260\n",
		      q,index);
	       turing_starta=turing_start+1;
	       turing_enda=turing_end+1;
	       if((sign(re_zs1[turing_starta])==UNK)||(sign(re_zs2[turing_starta])==UNK))
		 printf("Two unknowns in succession at turing_start.\n");
	       gamma_aa=gamma_a+int_double(2516,2517)/1000; // diff between int 1250..1260
	       gammaa=gamma+int_double(2516,2517)/1000;     // and 1250+5/64..1260+5/64
	     }
	   else
	     {
	       turing_starta=turing_start;
	       turing_enda=turing_end;
	       gamma_aa=gamma_a;
	       gammaa=gamma;
	     }
	 }
       for(i=0;i<UP_SAMPLE_IN_WIDTH-2;i++)
	 {
	   in_array[i]=re_zs1[i].left-re_zs1[i].right;
	   in_array[i]/=2.0;
	   printf("%10.8e\n",real(in_array[i]));
	 }
       exit(0);
       for(i=0;i<UP_SAMPLE_IN_WIDTH/2-1;i++)
	 {
	   in_array[UP_SAMPLE_IN_WIDTH/2-i-1]*=gaussians[i];
	   in_array[UP_SAMPLE_IN_WIDTH/2+i]*=gaussians[i];
	 }
       in_array[0]*=gaussians[i];
       in_array[UP_SAMPLE_IN_WIDTH-1]=(in_array[0]+in_array[UP_SAMPLE_IN_WIDTH-2])/2.0;
       
       fftw_execute(p_fwd);
       for(i=0;i<UP_SAMPLE_IN_WIDTH;i++)
	 printf("%10.8e\n",abs(in_array[i]));
       exit(0);
     }
   
 }
