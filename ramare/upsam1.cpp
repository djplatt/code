//
// upsam1.cpp
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
//#include "fftw3.h"
using namespace std;

#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"

#define UP_SAMPLE_RATE (2048)
#define FFT_IN_LEN (512) // = 40 in t
#define FFT_OUT_LEN (FFT_IN_LEN*UP_SAMPLE_RATE)
#define FFT_SACRIFICE (3*FFT_IN_LEN/8)

#include "../includes/upsample.h"

int_complex Fft_in[FFT_IN_LEN],Fft_out[FFT_OUT_LEN],ws_in[FFT_IN_LEN/2],ws_out[FFT_OUT_LEN/2],Gaussians[FFT_IN_LEN/2];


long unsigned int simple_count(int_double *re_zs,int len)
{
  long unsigned int k=0,count=0;
  while(sign(re_zs[k])==UNK) k++;
  sign_t last_sign=sign(re_zs[k++]);
  while(k<=len)
    {
      sign_t this_sign=sign(re_zs[k++]);
      if((this_sign==UNK)||(this_sign==last_sign))
	continue;
      count++;
      last_sign=this_sign;
    }
  return(count);
}

void init()
{
  int i;
  int_double two_pi_n;
  ws_in[0]=c_one;
  two_pi_n=int_double(2.0)/FFT_IN_LEN;
  for(i=1;i<FFT_IN_LEN/2;i++)
    sin_cospi(two_pi_n*i,&ws_in[i].imag,&ws_in[i].real);
  ws_out[0]=c_one;
  two_pi_n=int_double(2.0)/FFT_OUT_LEN;
  for(i=1;i<FFT_OUT_LEN/2;i++)
    sin_cospi(two_pi_n*i,&ws_out[i].imag,&ws_out[i].real);
}



// perform an in place FFT ,n a power of 2
// F[k]=sum_{j=0..N-1} F[j]e(-kj/N)
void fft(int_complex *x,unsigned int n,int_complex *w)
{
	unsigned int i,j,k,l;
	int_complex *p,*xend=x+n,ctemp;

	/* swap each element with one with bit-reversed index */
	for (i=0,l=n>>1;i<l;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k)
				j |= 1;
		}
		if (i < j)
		{
			ctemp=x[i];
			x[i]=x[j];
			x[j]=ctemp;
		}
		else if (i > j)
		{
			ctemp=x[n-1-i];
			x[n-1-i]=x[n-1-j];
			x[n-1-j]=ctemp;
		}
		++i, j |= l;
		ctemp=x[i];
		x[i]=x[j];
		x[j]=ctemp;
	}

	for (k=1,l=n/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<n/2;j+=l,p++) {
				ctemp=p[k]*w[j];
				p[k]=p[0]-ctemp;
				p[0]+=ctemp;
			}
}

/* perform an in place inverse FFT */
void ifft(int_complex *x,unsigned int n,int_complex *w) {
	unsigned int i,l=n>>1;
	int_complex ctemp;

	fft(x,n,w);
	x[0]=x[0]/n;
	x[l]=x[l]/n;

	for (i=1;i<l;i++) {
		x[i]=x[i]/n;
		x[n-i]=x[n-i]/n;
		ctemp=x[i];
		x[i]=x[n-i];
		x[n-i]=ctemp;
	}
}

void up_sample(unsigned long int start, int_double *re_zs1, int_double *re_zs2)
{
  int in_ptr=0,i;
  if(start==0) // use re_zs2
    {
      for(i=FFT_SACRIFICE-1;i>=0;i--,in_ptr++)
	Fft_in[in_ptr]=re_zs2[i];
    }
  for(i=0;in_ptr<FFT_IN_LEN;i++,in_ptr++)
    Fft_in[in_ptr]=re_zs1[i];
  for(int k=0;k<FFT_IN_LEN;k++)
    {printf("%f ",(k-FFT_SACRIFICE)*5.0/64.0);print_int_complex_str("",Fft_in[k]);}
  printf("\n");
  double t=5.0/128.0;
  for(int k=0;k<FFT_IN_LEN/2;k++,t+=5.0/64.0)
    {
      int_double g=exp(-sqr(int_double(t))/10);
      Fft_in[FFT_IN_LEN/2+k]*=g;
      Fft_in[FFT_IN_LEN/2-1-k]*=g;
    }
  for(int k=0;k<FFT_IN_LEN;k++)
    print_int_complex_str("",Fft_in[k]);
  printf("\n");

  // multiply by Gaussians here
  fft(Fft_in,FFT_IN_LEN,ws_in);
  for(int k=0;k<FFT_IN_LEN;k++)
    print_int_complex_str("",Fft_in[k]);
  printf("\n");
  for(i=0;i<FFT_IN_LEN/2;i++)
    {
      Fft_out[i]=Fft_in[i];
      Fft_out[FFT_OUT_LEN-1-i]=Fft_in[FFT_IN_LEN-1-i];
    }
  for(;i<FFT_OUT_LEN-FFT_IN_LEN/2;i++)
    Fft_out[i]=c_zero; // should be a little error term
  for(int k=0;k<FFT_OUT_LEN;k++)
    print_int_complex_str("",Fft_out[k]);
  printf("\n");
  ifft(Fft_out,FFT_OUT_LEN,ws_out);
  for(int k=0;k<FFT_OUT_LEN;k++)
    {printf("%f ",(k-FFT_SACRIFICE*UP_SAMPLE_RATE)*5.0/64.0/UP_SAMPLE_RATE);print_int_complex_str("",Fft_out[k]);}
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
   init();

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
       //for(int k=turing_start;k<turing_end;k++) print_int_double_str("re_zs1[]=",re_zs1[k]); exit(0);
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
	   printf("Looking for %d zeros.\n",diff);continue;
	   cz1=rig_num_zeros_first(0,turing_starta,re_zs1,re_zs1);
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
	   if(diff==BAD_DIFF) // Error in Turing Zone
	     {
	       printf("q:%d index:%d (cmplx) error in Turing Zone\n",q,index);
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs2,out_file);
	       continue;
	     }
	   printf("Looking for %d zeros.\n",diff);

	   cz1=rig_num_zeros_first(0,turing_starta,re_zs1,re_zs2);
	   cz2=rig_num_zeros_first(0,turing_starta,re_zs2,re_zs1);
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
