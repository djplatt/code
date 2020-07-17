//
// low_upsammoreC_find_unks.cpp
//
// upsample at x512 and list unknown regions
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

#define UP_SAMPLE_RATE (512)

#include "../includes/upsample.h"


sign_t get_sign(int_double *re_zs1, int_double *re_zs2, int i, int j)
{
  if(j==0)
    return(sign(re_zs1[i]));
  else
    return(check_sign(i,j,re_zs1,re_zs2));
}

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
   int diff,t_sam_start,last_i,last_j;
   int diff1,cz1,cz2,no_zeros=0,len;
   sign_t this_sign,last_sign;
   int in_index;
   _fpu_rndd();

   if(argc!=4)
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

   in_index=atol(argv[3]);

   //printf("Running upsammoreC, upsampling at x %d Upsample Error=%e.\n",UP_SAMPLE_RATE,d_inter_err);

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
       printf("q=%u index=%u num_s=%u.\n",q,index,num_s);
       if(first_time)
	 {
	   first_time=false;

	   //im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*num_s,16);
	   //if(!im_s_vec)
	   //fatal_error("Fatal error allocating memory for im_s_vec. Exting.\n");
	   //im_s_vec[0].im_s=0.0;
	   //for(int n=1;n<num_s;n++)
	   //im_s_vec[n].im_s=im_s_vec[n-1].im_s+one_over_two_B;
   
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
       set_d_inter_err(q,turing_enda*one_over_two_B);

       fread(&gamma,sizeof(int_double),1,infile);


       if(!real_p)
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
	 }
       if(index!=in_index)
	 continue;
       if(real_p)
	 {
	   last_sign=sign(re_zs1[0]);
	   last_i=0;
	   last_j=0;
	   for(i=0;i<turing_starta;i++)
	     for(j=0;j<UP_SAMPLE_RATE;j++)
	       {
		 this_sign=get_sign(re_zs1,re_zs1,i,j);
		 if(this_sign!=last_sign)
		   {
		     if(last_sign==UNK) // end of unknown bit
		       printf("Real:- unknown section found from %d %d to %d %d\n",last_i,last_j,i,j);
		     if(this_sign==UNK) // new unkown bit
		       {
			 last_i=i;
			 last_j=j;
		       }
		     last_sign=this_sign;
		   }
	       }
	 }
       else
	 {
	   last_sign=sign(re_zs1[0]);
	   last_i=0;
	   last_j=0;
	   for(i=0;i<turing_starta;i++)
	     for(j=0;j<UP_SAMPLE_RATE;j++)
	       {
		 this_sign=get_sign(re_zs1,re_zs2,i,j);
		 if(this_sign!=last_sign)
		   {
		     if(last_sign==UNK) // end of unknown bit
		       {
			 len=(i-last_i)*UP_SAMPLE_RATE+j-last_j;
			 if(len>10)
			   printf("Cmplx:- unknown section found from %d %d length %d\n",last_i,last_j,len);
		       }
		     if(this_sign==UNK) // new unkown bit
		       {
			 last_i=i;
			 last_j=j;
		       }
		     last_sign=this_sign;
		   }
	       }
	   last_sign=sign(re_zs2[0]);
	   last_i=0;
	   last_j=0;
	   for(i=0;i<turing_starta;i++)
	     for(j=0;j<UP_SAMPLE_RATE;j++)
	       {
		 this_sign=get_sign(re_zs2,re_zs1,i,j);
		 if(this_sign!=last_sign)
		   {
		     if(last_sign==UNK) // end of unknown bit
		       {
			 len=(i-last_i)*UP_SAMPLE_RATE+j-last_j;
			 if(len>10)
			   {
			     printf("Cmplx_bar:- unknown section found from %d %d length %d\n",last_i,last_j,len);
			     //print_int_double_str("F(t)=",re_zs2[last_i]);
			   }
		       }
		     if(this_sign==UNK) // new unkown bit
		       {
			 last_i=i;
			 last_j=j;
		       }
		     last_sign=this_sign;
		   }
	       }
	 }
     }   
 }
 
