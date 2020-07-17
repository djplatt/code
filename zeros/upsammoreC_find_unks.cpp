//
// upsammoreC_find_unks.cpp
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

#define UP_SAMPLE_RATE (4096)

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

   if(fread(&num_s,sizeof(int),1,infile)!=1)
     {
       printf("Empty data file. Exiting.\n");
       exit(0);
     }
 
   printf("File contains %d data points.\n",num_s);
   re_zs1=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
   if(!re_zs1) 
     fatal_error("Fatal error allocating memory for re_zs1. Exting.\n");
  
   re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
   if(!re_zs2)
     fatal_error("Fatal error allocating memory for re_zs2. Exting.\n");

   
   out_file=fopen(argv[2],"wb");
   if(!out_file)
     {
       printf("Failed to open %s for binary output. Exiting.\n",argv[2]);
       exit(1);
     }

   in_index=atol(argv[3]);

   //printf("Running upsammoreC, upsampling at x %d Upsample Error=%e.\n",UP_SAMPLE_RATE,d_inter_err);

   fread(&gamma,sizeof(int_double),1,infile);
   fread(&gamma_a,sizeof(int_double),1,infile);
   fread(&gamma_del,sizeof(int_double),1,infile);
   fread(&gamma_a_del,sizeof(int_double),1,infile);
   //fread(&foo,sizeof(int_double),1,infile);
   fread(&turing_starta,sizeof(int),1,infile);
   fread(&turing_enda,sizeof(int),1,infile);   

   rig_setup(im_s_vec);

   while(fread(&q,sizeof(unsigned int),1,infile))
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
       if(real_p) printf("Character is real.\n"); else printf("Character is complex.\n");
       fread(&n_zeros,sizeof(int),1,infile);  // not used
       printf("n_zeros read as %d\n",n_zeros);
       fread(&neg_one,sizeof(bool),1,infile);
       if(neg_one) printf("Character is odd.\n"); else printf("Character is even.\n");
       fread(&omega,sizeof(int_complex),1,infile);
       print_int_complex_str("Omega=",omega);
       fread(re_zs1,sizeof(int_double),num_s,infile);
       arg_omega=argument(omega);
       if(!real_p)
	 {
	   fread(&neg_one2,sizeof(bool),1,infile);
	   fread(&omega2,sizeof(int_complex),1,infile);
	   print_int_complex_str("Omega_bar=",omega2);
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
       /*
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
	 }
       */
       if(index!=in_index)
	 {
	   printf("Skipping index %d\n",index);
	   continue;
	 }
       if(real_p)
	 {
	   last_sign=sign(re_zs1[0]);
	   last_i=0;
	   last_j=0;
	   for(i=0;i<turing_enda;i++)
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
	   for(i=0;i<turing_enda;i++)
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
	   // now do the conjugate
	   last_sign=sign(re_zs2[0]);
	   last_i=0;
	   last_j=0;
	   for(i=0;i<turing_enda;i++)
	     for(j=0;j<UP_SAMPLE_RATE;j++)
	       {
		 this_sign=get_sign(re_zs2,re_zs1,i,j);
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
	 }
     }
 }
