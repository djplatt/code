//
// low_upsammoreC_print.cpp
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

#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"

#define UP_SAMPLE_RATE (1)

#include "../includes/upsample.h"

#define INDEX (1573)
#define N0_START (0)
#define N0_END (130560)

void print_change(int counter, sign_t this_sign, sign_t last_sign, int n0, int offset)
{
  printf("%d changing from ",counter);print_sign(last_sign);printf(" to ");print_sign(this_sign);
  printf(" at %d %d %f\n",n0,offset,((double)n0+(double)offset/(double)UP_SAMPLE_RATE)*one_over_two_B);
}

 int main(int argc, char **argv)
 {
   FILE *infile,*out_file;
   im_s *im_s_vec;
   int counter,i,j,num_s,index,index2,n_zeros;
   int q,turing_start,turing_end,turing_starta,turing_enda;
   int_double gamma,gamma_a,arg_omega,arg_omega2;
   int_double *re_zs1,*re_zs2,gammaa,gamma_aa;
   int_double gamma_del,gamma_a_del;
   bool real_p,neg_one,neg_one2,first_time=true;
   int_complex omega,omega2;
   int diff,t_sam_start,n0,offset;
   int diff1,cz1,cz2,no_zeros=0;
   sign_t this_sign,last_sign;

   _fpu_rndd();

   if(argc!=4)
     {
       printf("usage infile outfile error\n");
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

   //d_inter_err=atof(argv[3]);

   //   printf("Running upsammoreC, upsampling at x %d Upsample Error=%e.\n",UP_SAMPLE_RATE,d_inter_err);

   while(fread(&q,sizeof(unsigned int),1,infile))
     {
       if(q>MAX_Q)
	 {
	   printf("q=%d exceeds MAX_Q. Exiting.\n",q);
	   exit(1);
	 }

       fread(&index,sizeof(unsigned int),1,infile);
       printf("Index=%d\n",index);
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
	   if(INDEX!=index)
	     continue;
	   /*	   last_sign=check_sign(N0_START,0,re_zs1,re_zs1);
	   for(n0=N0_START;n0<N0_END;n0++)
	     {
	       for(offset=1;offset<=UP_SAMPLE_RATE;offset++)
		 {
		   this_sign=check_sign(n0,offset,re_zs1,re_zs1);
		   if((this_sign!=last_sign)&&(this_sign!=UNK))
		     {
		       printf("sign change at %10.8e\n",((double)n0+(double)offset/(double)UP_SAMPLE_RATE)*one_over_two_B);
		       last_sign=this_sign;
		     }
		 }
		 }*/
	   counter=0;
	   last_sign=sign(re_zs1[N0_START]);
	   for(n0=N0_START;n0<N0_END;n0++)
	     {
	       this_sign=sign(re_zs1[n0]);
	       if(last_sign!=this_sign)
		 {
		   print_change(++counter,this_sign,last_sign,n0,0);
		   last_sign=this_sign;
		 }
	       for(offset=1;offset<UP_SAMPLE_RATE;offset++)
		 {
		   this_sign=check_sign(n0,offset,re_zs1,re_zs2);
		   if(last_sign!=this_sign)
		     {
		       print_change(++counter,this_sign,last_sign,n0,offset);
		       last_sign=this_sign;
		     }
		 }
	     }

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
	   if(INDEX!=index)
	     continue;
	   /*
	   int counter=0;
	   bool unk=false;
	   last_sign=sign(re_zs1[N0_START]);
	   for(n0=N0_START;n0<N0_END;n0++)
	     {
	       this_sign=sign(re_zs1[n0]);
	       if(this_sign==UNK)
		 {
		   unk=true;
		 }
	       else
		 {
		   if(unk)
		     {
		       unk=false;
		       if(this_sign==last_sign)
			 printf("unknown portion at %d %d %10.8e\n",n0,offset,((double)n0+(double)offset/(double)UP_SAMPLE_RATE)*one_over_two_B);
		     }
		   else
		     {
		       if(this_sign!=last_sign)
			 {
			   counter++;
			   printf("sign change %d at %d %d %10.8e\n",counter,n0,offset,((double)n0+(double)offset/(double)UP_SAMPLE_RATE)*one_over_two_B);
			   last_sign=this_sign;
			 }
		     }
		 }
	       for(offset=1;offset<UP_SAMPLE_RATE;offset++)
		 {
		   this_sign=check_sign(n0,offset,re_zs1,re_zs2);
		   if(this_sign==UNK)
		     {
		       unk=true;
		     }
		   else
		     {
		       if(unk)
			 {
			   unk=false;
			   if(this_sign==last_sign)
			     printf("unknown portion at %d %d %10.8e\n",n0,offset,((double)n0+(double)offset/(double)UP_SAMPLE_RATE)*one_over_two_B);
			 }
		       else
			 {
			   if(this_sign!=last_sign)
			     {
			       counter++;
			       printf("sign change %d at %d %d %10.8e\n",counter,n0,offset,((double)n0+(double)offset/(double)UP_SAMPLE_RATE)*one_over_two_B);
			       last_sign=this_sign;
			     }
			 }
		     } 
		 }
	     }
	 }
     }
 }
	   */
	   counter=0;
	   printf("re_zs1[%d]=",N0_START);print_int_double_str("",re_zs1[N0_START]);
	   last_sign=sign(re_zs1[N0_START]);
	   for(n0=N0_START;n0<N0_END;n0++)
	     {
	       this_sign=sign(re_zs1[n0]);
	       if(last_sign!=this_sign)
		 {
		   print_change(++counter,this_sign,last_sign,n0,0);
		   last_sign=this_sign;
		 }
	       for(offset=1;offset<UP_SAMPLE_RATE;offset++)
		 {
		   this_sign=check_sign(n0,offset,re_zs1,re_zs2);
		   if(last_sign!=this_sign)
		     {
		       print_change(++counter,this_sign,last_sign,n0,offset);
		       last_sign=this_sign;
		     }
		 }
	     }
	   counter=0;
	   printf("re_zs2[%d]=",N0_START);print_int_double_str("",re_zs2[N0_START]);

	   last_sign=sign(re_zs2[N0_START]);
	   for(n0=N0_START;n0<N0_END;n0++)
	     {
	       this_sign=sign(re_zs2[n0]);
	       if(last_sign!=this_sign)
		 {
		   print_change(++counter,this_sign,last_sign,n0,0);
		   last_sign=this_sign;
		 }
	       for(offset=1;offset<UP_SAMPLE_RATE;offset++)
		 {
		   this_sign=check_sign(n0,offset,re_zs2,re_zs1);
		   if(last_sign!=this_sign)
		     {
		       print_change(++counter,this_sign,last_sign,n0,offset);
		       last_sign=this_sign;
		     }
		 }
	     }
	 }
     }
 }
     
