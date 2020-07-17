//
// summarise_errorC.cpp
//
// read in a C format file
// and summarise the width
// of the intervals contained

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include "../includes/int_double12.0.h"

void fatal_error(const char *str)
{
  printf("%s. Exiting.\n",str);
  exit(0);
}

int main(int argc, char **argv)
 {
   FILE *infile,*out_file;
   int i,j,num_s,index,index2,n_zeros;
   int q,turing_start,turing_end,turing_starta,turing_enda;
   int_double gamma,gamma_a,arg_omega,arg_omega2;
   int_double *re_zs1,*re_zs2,gammaa,gamma_aa;
   int_double gamma_del,gamma_a_del;
   bool real_p,neg_one,neg_one2,first_time=true;
   int_complex omega,omega2;
   int diff,t_sam_start,last_i,last_j;
   int diff1,cz1,cz2,no_zeros=0,len;
   int in_index;
   _fpu_rndd();

   if(argc!=2)
     {
       printf("Usage:- summarise_error <infile>. Exiting.\n");
       exit(0);
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
     fatal_error("Fatal error allocating memory for re_zs1");
  
   re_zs2=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
   if(!re_zs2)
     fatal_error("Fatal error allocating memory for re_zs2");

   fread(&gamma,sizeof(int_double),1,infile);
   fread(&gamma_a,sizeof(int_double),1,infile);
   fread(&gamma_del,sizeof(int_double),1,infile);
   fread(&gamma_a_del,sizeof(int_double),1,infile);
   //fread(&foo,sizeof(int_double),1,infile);
   fread(&turing_starta,sizeof(int),1,infile);
   fread(&turing_enda,sizeof(int),1,infile);   

   while(fread(&q,sizeof(unsigned int),1,infile))
     {

       fread(&index,sizeof(int),1,infile);
       printf("processing q=%d index %d\n",q,index);
       fread(&real_p,sizeof(bool),1,infile);
       if(real_p) printf("Character is real.\n"); else printf("Character is complex.\n");
       fread(&n_zeros,sizeof(int),1,infile);  // not used
       //printf("n_zeros read as %d\n",n_zeros);
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
       int r_err=-1000,a_err=-1000;
       for(long unsigned int j=0;j<num_s;j++)
	 {
	   int this_r_err=rel_error(re_zs1[j]);
	   int this_a_err=abs_error(re_zs1[j]);
	   if(this_r_err>r_err) r_err=this_r_err;
	   if(this_a_err>a_err) a_err=this_a_err;
	 }
       if(!real_p)
	 for(long unsigned int j=0;j<num_s;j++)
	   {
	     int this_r_err=rel_error(re_zs2[j]);
	     int this_a_err=abs_error(re_zs2[j]);
	     if(this_r_err>r_err) r_err=this_r_err;
	     if(this_a_err>a_err) a_err=this_a_err;
	   }

       printf("relative error 10^{%d} absolute error 10^{%d}\n",r_err,a_err);
     }
 }
