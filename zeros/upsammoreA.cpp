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
#include <unistd.h>
#include "fftw3.h"
using namespace std;

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
   
   out_file=fopen(argv[2],"wb");
   if(!out_file)
     {
       printf("Failed to open %s for binary output. Exiting.\n",argv[2]);
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
   fwrite(&num_s,sizeof(int),1,out_file);
   fwrite(&gamma,sizeof(int_double),1,out_file);
   fwrite(&gamma_a,sizeof(int_double),1,out_file);
   fwrite(&gamma_del,sizeof(int_double),1,out_file);
   fwrite(&gamma_a_del,sizeof(int_double),1,out_file);
   fwrite(&turing_start,sizeof(int),1,out_file);
   fwrite(&turing_end,sizeof(int),1,out_file);

   //printf("%d %d %d\n",num_s,turing_start,turing_end);
   //print_int_double_str("int(Im(Log(Gamma(s/2))))=",gamma);
   //print_int_double_str("int(Im(Log(Gamma((s+1)/2))))=",gamma_a);


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
   fwrite(im_s_vec,sizeof(im_s),num_s,out_file);
   setup(im_s_vec);

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
	   cz1=num_zeros1(0,turing_starta,re_zs1);
	   no_zeros+=cz1;
	   // how many should there be?
	   diff1=check_zeros(cz1,turing_starta,turing_enda,q,(neg_one ? gamma_aa : gammaa),
			     arg_omega,re_zs1,im_s_vec,index,t_sam_start);
	   if(diff1==BAD_DIFF) // Error in Turing Zone
	     {
	       //printf("q:%d index:%d check_zeros (real) returned:%d\n",q,index,diff1);
	       save_state(q,index,true,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs1,out_file);
	       continue;
	     }
	   if(diff1!=0)
	     diff1=find_more_zeros(re_zs1,re_zs1,turing_starta,diff1,im_s_vec);
	   if(diff1!=0) // missed some zeros
	     {
	       no_zeros+=diff1;
	       save_state(q,index,true,neg_one,omega,neg_one2,omega2,cz1+diff1,num_s,re_zs1,re_zs1,out_file);
	       printf("Do some more checking on q:%d, index:%d, (real) difference:%d\n",q,index,diff1);
	     }
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
		   printf("Two unknowns in succession at turing_start.\n");
		   save_state(q,index,false,neg_one,omega,neg_one2,omega2,
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
	       
	   cz1=num_zeros1(0,turing_starta,re_zs1);//,im_s_vec,re_zs2);
	   cz2=num_zeros1(0,turing_starta,re_zs2);//,im_s_vec,re_zs1);
	   
	   printf("cz1 %d cz2 %d\n",cz1,cz2);
	   
	   no_zeros+=cz1+cz2;
	   diff1=check_zeros(cz1,turing_starta,turing_enda,q,(neg_one ? gamma_aa : gammaa),
			     arg_omega,re_zs1,im_s_vec,index,t_sam_start);
	   if(diff1==BAD_DIFF) // problem in Turing method
	     {
	       printf("q:%d index:%d check_zeros (cmplx) returned:%d\n",q,index,diff1);
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs2,out_file);
	       continue;
	     }
	   diff=check_zeros(cz2,turing_starta,turing_enda,q,(neg_one2 ? gamma_aa : gammaa),
			    arg_omega2,re_zs2,im_s_vec,index2,t_sam_start);
	   if(diff==BAD_DIFF) // problem in Turing method with conjugate
	     {
	       printf("q:%d index:%d check_zeros (cmplx) returned:%d\n",q,index,diff1);
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs2,out_file);
	       continue;
	     }
	   printf("diff1=%d diff=%d\n",diff1,diff);
	   diff1=find_more_zeros(re_zs1,re_zs2,turing_starta,diff1,im_s_vec);
	   printf("diff1=%d\n",diff1);
	   diff=find_more_zeros(re_zs2,re_zs1,turing_starta,diff,im_s_vec);
	   printf("diff=%d\n",diff);
	   diff1+=diff;
	   
	   if(diff1<0) // what th...
	     {
	       printf("q:%d index: %d %d too many zeros located.\n",q,index,index2);
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs2,out_file);
	     }
	   if(diff1>0) // missed some zeros
	     {
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,cz1+cz2+diff1,num_s,re_zs1,re_zs2,out_file);
	       printf("Do some more checking on q:%d, index:%d, (cmpx) difference:%d\n",q,index,diff1);
	       no_zeros+=diff1;
	     }
	 }
       //printf("There were %d zeros for q=%d and index=%d T<=%f.\n",no_zeros,q,index,im_s_vec[turing_start].im_s);
       //printf("Stat_count was %d.\n",stat_count); 
       //stat_count=0;
       no_zeros=0;
       //printf("Check_sign was called %d times.\n",check_count);
       check_count=0;
     }
   if(!save_p)
     {
       fclose(out_file);
       out_file=fopen(argv[2],"w");
       fclose(out_file);
     }

 }
