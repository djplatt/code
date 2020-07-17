//
// upsammore.cpp
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

typedef complex<double> mcomplex;

#define MAX_Q (100000)
#define POS 0
#define NEG 1
#define UNK 2
#define UP 0
#define DOWN 1
#define UP_SAMPLE_RATE (8192) // factor by which we upsample
#define UP_SAMPLE_IN_WIDTH (512)
#define UP_SAMPLE_OUT_WIDTH (UP_SAMPLE_IN_WIDTH*UP_SAMPLE_RATE)
#define UP_SAMPLE_H (90.0) // controls decay of Gaussian, small h, fast decay
#define UP_SAMPLE_USABLE_WIDTH (UP_SAMPLE_IN_WIDTH/2)
#define UP_SAMPLE_SACRIFICE ((UP_SAMPLE_IN_WIDTH-UP_SAMPLE_USABLE_WIDTH)/2)
#define FFTW_PLAN_MODE FFTW_MEASURE // we only plan once, so get the best plan poss

#define INTER_H (7.0/32.0)
#define TWO_H_SQR (2.0*INTER_H*INTER_H)
#define one_over_two_B (5.0/64.0)
#define INTER_N (20)
#define d_inter_err ((double) 1.06e-7)
#define BAD_DIFF (INT_MAX)

#define sign_t int

bool testing=false;

void fatal_error (const char *str)
{
  printf("%s",str);
  exit(1);
}

mcomplex *in_array;
mcomplex *out_array;
fftw_plan p_fwd,p_bkwd;

// this contains the gaussian factors to use when upsampling
double gaussians[UP_SAMPLE_IN_WIDTH/2];


int_double expsincs[(UP_SAMPLE_RATE-1)*INTER_N*2];
int_double twoBpi;
// consists of UP_SAMPLE_RATE-1 "rows" of gauss*sinc
// first row centred at t0-n/2B = 1/UP_SAMPLE_RATE
// last row at                  = (UP_SAMPLE_RATE-1)/UP_SAMPLE_RATE
void setexpsincs()
{
	int offset,i,ptr=0;
	double del;
	int_double ex,si;

	for(offset=1;offset<UP_SAMPLE_RATE;offset++)
	{
	  del=(((double) offset)/((double) UP_SAMPLE_RATE)+(double) INTER_N - 1.0)*one_over_two_B;
		for(i=0;i<INTER_N*2;i++,del-=one_over_two_B)
		  {
		    
		    ex=exp(-int_double(del*del)/TWO_H_SQR);
		    
		    si=sinc(twoBpi*del);
		    
		    expsincs[ptr++]=ex*si;
		  }
	}
}

// please note double exp(double) does not work if you mess with rounding
// ergo use double exp1(double)
void set_gaussians()
{
  double i_over_h;
	for(int i=0;i<UP_SAMPLE_IN_WIDTH/2;i++)
	{
		i_over_h=(double) i/UP_SAMPLE_H;
		gaussians[i]=exp1(-(i_over_h*i_over_h));
	}
}

double exppittwos[UP_SAMPLE_SACRIFICE];
int_double exppittwos_d[UP_SAMPLE_SACRIFICE];

void setexppittwos(im_s *im_s_vec)
{
  for(int i=0;i<UP_SAMPLE_SACRIFICE;i++)
    {
      exppittwos[i]=exp1(-im_s_vec[i+1].im_s*d_pi.left/2);
      exppittwos_d[i]=exp(d_pi*im_s_vec[i+1].im_s/(-2.0));
    }
}

void setup(im_s* im_s_vec)
{
  twoBpi=d_pi/one_over_two_B;
	set_gaussians();
	setexpsincs();
	setexppittwos(im_s_vec);
	in_array=(mcomplex*) fftw_malloc(sizeof(fftw_complex) * UP_SAMPLE_IN_WIDTH);
	out_array=(mcomplex*) fftw_malloc(sizeof(fftw_complex) * UP_SAMPLE_OUT_WIDTH);

	p_fwd=fftw_plan_dft_1d(UP_SAMPLE_IN_WIDTH, reinterpret_cast<fftw_complex*>(in_array),
		reinterpret_cast<fftw_complex*>(in_array), FFTW_FORWARD, FFTW_PLAN_MODE);

	p_bkwd=fftw_plan_dft_1d(UP_SAMPLE_OUT_WIDTH, reinterpret_cast<fftw_complex*>(out_array),
		reinterpret_cast<fftw_complex*>(out_array), FFTW_BACKWARD, FFTW_PLAN_MODE);
}

// in_array[0..UP_SAMPLE_IN_WIDTH-2] contains the sampled data
// this is multiplied by the Gaussian
// the in_array[UP_SAMPLE_WIDTH-1] is set to avg of [0], [U_S_W-2]
// then FFT'd in place
// then zero padded into out_array
// then inverse FFT'd inplace
void up_sample ()
{
  int i,j;
	for(i=0;i<UP_SAMPLE_IN_WIDTH/2-1;i++)
	{
		in_array[UP_SAMPLE_IN_WIDTH/2-i-1]*=gaussians[i];
		in_array[UP_SAMPLE_IN_WIDTH/2+i]*=gaussians[i];
	}
	in_array[0]*=gaussians[i];
	in_array[UP_SAMPLE_IN_WIDTH-1]=(in_array[0]+in_array[UP_SAMPLE_IN_WIDTH-2])/2.0;

	fftw_execute(p_fwd);

	// put in first half of frequencies
	for(i=0;i<(UP_SAMPLE_IN_WIDTH>>1);i++)
	  out_array[i]=in_array[i];
	// remember where we are in in_array
	j=i;
	// zero out central part of frequencies
	for(;i<(UP_SAMPLE_OUT_WIDTH-(UP_SAMPLE_IN_WIDTH>>1));i++)
	  out_array[i]=mcomplex(0.0,0.0);
	// put in second half of frequencies
	for(;i<UP_SAMPLE_OUT_WIDTH;i++)
		out_array[i]=in_array[j++];

	fftw_execute(p_bkwd);
}

inline sign_t sign(int_double &x)
{
  if(x.left>0.0)
    return(POS);
  if(x.right>0.0)
    return(NEG);
  return(UNK);
}

inline sign_t d_sign(double t)
{
  if(t<0.0)
    return(NEG);
  if(t>0.0)
    return(POS);
  return(UNK);
}

void print_sign(sign_t sign)
{
  if(sign==POS)
    printf("+");
  else
    {
      if(sign==NEG)
	printf("-");
      else
	printf("?\n");
    }
}

void d_print_signs( int n, mcomplex *re_zs,  int step)
{
  int i;
  sign_t last_sign,this_sign;
  last_sign=d_sign(real(re_zs[0]));
  print_sign(last_sign);
  for(i=step;i<n;i+=step)
    {
      this_sign=(d_sign(real(re_zs[i])));
      if(this_sign!=last_sign)
	{
	  //printf("\n");
	  last_sign=this_sign;
	}
      print_sign(this_sign);
    }
  printf("\n");
}

void print_signs( int n, int_double *re_zs,  int step)
{
  int i;
  print_sign(sign(re_zs[0]));
  for(i=step;i<n;i+=step)
    print_sign(sign(re_zs[i]));
  printf("\n");
}

void print_zeros(int n, int_double *re_zs,im_s *im_s_vec)
{
  int i=1;
  sign_t this_sign,last_sign=sign(re_zs[0]);
  
  while(i<n)
    {
      this_sign=sign(re_zs[i]);
      if(this_sign==UNK)
	{i++;continue;}
      if(this_sign!=last_sign)
	{
	  printf("%7.4f\n",im_s_vec[i].im_s);
	  last_sign=this_sign;
	}
      i++;
    }
}

int check_count=0;

// check that the sign of Z(t) is as suspected = rough_sign
inline bool check_sign(int n0, int offset, sign_t rough_sign, int_double *re_zs1, int_double *re_zs2)
{
  int n,n1;
  int_double *expsinc,res=int_double(-d_inter_err,d_inter_err);
  bool bres;
  check_count++;
  //printf("In check sign with n0=%d, offset =%d\n",n0,offset);
  while(offset>UP_SAMPLE_RATE)
    {
      n0++;
      offset-=UP_SAMPLE_RATE;
    }
  // t is now between [n0] and [n0+1]
  expsinc=&expsincs[(offset-1)*INTER_N*2];
  if(n0-INTER_N+1<0)
    {
      //if(sign(re_zs1[0])==sign(re_zs2[0]))
	for(n=n0-INTER_N+1,n1=0;n<0;n++,n1++)
	  res+=re_zs2[-n]*expsinc[n1]*exppittwos_d[-n-1];
	//else
	//for(n=n0-INTER_N+1,n1=0;n<0;n++,n1++)
	//res+=(-re_zs2[-n])*expsinc[n1]*exppittwos_d[-n-1];
      for(n=0;n<=n0+INTER_N;n++,n1++)
	res+=re_zs1[n]*expsinc[n1];
      return(sign(res)==rough_sign);
    }
  for(n=n0-INTER_N+1,n1=0;n<=n0+INTER_N;n++,n1++)
    {
      res+=re_zs1[n]*expsinc[n1];
    }
  bres=(sign(res)==rough_sign);
  if(!bres)
    printf("Check_sign failed.\n");
  return(bres);
}


// return the number of zeros in out_array between 0 and UP_SAMPLE_RATE*width inclusive
// the signs of the end points are known to be st_sign and en_sign, neither of which is UNK
// re_z_ptr points at the Z value at the start of the sample
int up_num_zeros (mcomplex *out_array, im_s *im_s_vec, sign_t st_sign, sign_t en_sign, int width,
		  int re_z_ptr, int_double *re_zs1, int_double *re_zs2)
{
  int ptr=1,count=0,last_change=0,mid_pt;
  sign_t sign_bar,this_sign,last_sign;
  
  //printf("In up_num_zeros with ptr=%d sign= ",re_z_ptr);
  //print_sign(st_sign);printf("\n");
  //d_print_signs(width*UP_SAMPLE_RATE,out_array,1);
  if(en_sign==POS)
    sign_bar=NEG;
  else
    sign_bar=POS;
  last_sign=st_sign;
  for(ptr=1;ptr<=width*UP_SAMPLE_RATE;ptr++)
    {
      this_sign=d_sign(real(out_array[ptr]));
      if(this_sign==last_sign)
	continue;
      if(this_sign==UNK)
	continue;
      if(last_change==0) // first change
	{
	  count++; // or =1?
	  last_change=ptr;
	  last_sign=this_sign;
	  continue;
	}
      //printf("sign change at %d\n",ptr);
      mid_pt=(ptr+last_change-1)>>1;
      if(check_sign(re_z_ptr,mid_pt,sign_bar,re_zs1,re_zs2))
	count++;
      last_change=ptr;
      last_sign=this_sign;
    }
  //printf("up_num_zeros at %d returning %d\n",re_z_ptr,count);
  return(count);
}

inline void up_sample_load (int_double *re_zs)
{
  for(int i=0;i<UP_SAMPLE_IN_WIDTH-1;i++)
    in_array[i]=avg_int(re_zs[i]);
  up_sample();
}

inline void up_sample_load_neg (int_double *re_zs1, int_double *re_zs2, int ptr)
{
  int i;
  /*
  if(sign(re_zs1[0])!=sign(re_zs2[0]))
    for(i=0;i<UP_SAMPLE_SACRIFICE;i++)
      in_array[i]=-avg_int(re_zs2[UP_SAMPLE_SACRIFICE-i])*exppittwos[UP_SAMPLE_SACRIFICE-i-1];
  else
  */
    for(i=0;i<UP_SAMPLE_SACRIFICE;i++)
      in_array[i]=avg_int(re_zs2[UP_SAMPLE_SACRIFICE-i])*exppittwos[UP_SAMPLE_SACRIFICE-i-1];
  
  for(;i<UP_SAMPLE_IN_WIDTH-1;i++)
    in_array[i]=avg_int(re_zs1[i-UP_SAMPLE_SACRIFICE]);
  up_sample();
}

// work out the integral of N(T) from re_zs[a] to [b]
// assumes the section a to b fits into the usable
// portion of a single up_sample if re_zs[num_s] goes
// in the end
int_double up_n_twiddle( int a,  int b, int_double *re_zs, im_s *im_s_vec, int t_sam_start)
{
  int i,j,n,out_ptr=t_sam_start,n1=0;
  int_double end_point=int_double(im_s_vec[b].im_s),res=int_double(0.0);
  sign_t st_sign=sign(re_zs[a]),en_sign;

  while(st_sign==UNK)
      st_sign=sign(re_zs[++a]);
  while(sign(re_zs[b])==UNK)
    b--;
  //printf("Loading from re_zs[%d]\n",out_ptr);
  up_sample_load(&re_zs[out_ptr]);
  //for(i=0;i<UP_SAMPLE_OUT_WIDTH;i++)
  //  printf("%10.8e\n",real(out_array[i]));

  out_ptr=(a-out_ptr)*UP_SAMPLE_RATE;
  //d_print_signs((b-a)*UP_SAMPLE_RATE,&out_array[out_ptr],1);
  //printf("out_ptr=%d\n",out_ptr);
  i=a;
  j=a+1;
  while(j<=b)
    {
      en_sign=sign(re_zs[j]);
      while(en_sign==UNK)
	en_sign=sign(re_zs[++j]);
      n=up_num_zeros(&out_array[out_ptr],
		     im_s_vec,st_sign,en_sign,(j-i),i,re_zs,re_zs); // last arg is dummy here
      if(n!=0)
	{
	  //printf("Found zeros here\n");
	  //d_print_signs((j-i)*UP_SAMPLE_RATE,&out_array[out_ptr],1);
	  res+=(end_point-int_double(im_s_vec[i].im_s,im_s_vec[j].im_s))*n;
	  n1+=n;
	}
      out_ptr+=UP_SAMPLE_RATE*(j-i);
      i=j;
      j++;
      st_sign=en_sign;
    }
  //exit(0);
  //printf("Found %d zeros in Turing Zone.\n",n1);
  return(res);
}

int_double s_of_t( int q, double t2) // Rumely Th. 2 p 429
{
  int_double tmp2=int_double (t2),tmp=int_double(1242), res=int_double(18397);
  tmp/=10000;res/=10000;
  tmp2*=q;
  tmp2/=2;
  tmp2/=d_pi;
  res=res+tmp*log(tmp2);
  res.left=res.right;
  return(res);
}

// returns (t1+t2)*ln(q/pi)/4
inline int_double ln_term (double t1, double t2,  int q)
{
  return(log(int_double(q)/d_pi)*(t1+t2)/4);
}


sign_t up_sam(int ptr, int_double delta, int_double *re_zs)
{
  int_double res=int_double(-d_inter_err,d_inter_err);
  int_double del;
  int n;
  print_int_double_str("up_sam called with delta =",delta);
  printf("ptr= %d.\n",ptr);
  del=delta-one_over_two_B;
  for(n=1;n<=INTER_N;n++)
    {
      res+=re_zs[ptr+n]*exp(-del*del/TWO_H_SQR)*sinc(twoBpi*del);
      print_int_double_str("re_zs= ",re_zs[ptr+n]);
      print_int_double_str("res  = ",res);
      del-=one_over_two_B;
    }
  del=delta;
  for(n=0;n<INTER_N;n++)
    {
      res+=re_zs[ptr-n]*exp(-del*del/TWO_H_SQR)*sinc(twoBpi*del);
      print_int_double_str("re_zs= ",re_zs[ptr+n]);
      print_int_double_str("res  = ",res);
      del+=one_over_two_B;
    }
  print_int_double_str("up_sam computed ",res);
  exit(0);
  return(sign(res));
}

inline bool check_sign1(int ptr, int_double x, sign_t exp_sign, int_double *re_zs)
{
  //printf("In check_sign1 with ptr=%d\n.",ptr);
  if(x.right>0.0) // negative
    {
      ptr--;
      x+=1;
    }
  // now our point is between re_zs[ptr] and re_zs[ptr+1]
  //return(up_sam(ptr,x*one_over_two_B,re_zs)==exp_sign);
  return(check_sign(ptr,x.left*UP_SAMPLE_RATE,exp_sign,re_zs,re_zs));
}

// count zeros between a and b
inline int num_zeros1(int a, int b, int_double *re_zs)
{
  int zero_count=0,ptr=a+1;
  sign_t this_sign,last_sign=sign(re_zs[a]);
  while(last_sign==UNK)
    {
      last_sign=sign(re_zs[ptr++]);
      if(ptr>b)
	fatal_error("Catastrophic failure in num_zeros. Exiting.\n");
    }

  while(ptr<=b)
    {
      this_sign=sign(re_zs[ptr++]);
      while((this_sign==UNK)&&(ptr<=b))
	this_sign=sign(re_zs[ptr++]);
      if(this_sign==UNK) // file finished with one or more unknowns
	return(zero_count);
      if(this_sign!=last_sign)
	{
	  zero_count++;
	  last_sign=this_sign;
	}
    }
  return(zero_count);
}	

int stat_count=0;

inline int num_zeros2 (int b, int_double *re_zs,im_s *im_s_vec, int_double *re_zs2) 
{
  int n,count=0,ptr=0;
  int_double y1,y2=re_zs[ptr];
  int_double y3=re_zs[ptr+1];
  sign_t this_sign;
  while(ptr<=b-2)
    {
      y1=y2;
      y2=y3;
      y3=re_zs[ptr+2];
      this_sign=sign(y2);
      if((this_sign==POS)&&(y1>y2)&&(y3>y2)) // min above zero
	{ 
	  //printf("Positive stat point near %d\n",ptr);
	  stat_count++;
	  if(ptr-UP_SAMPLE_SACRIFICE<0)
	    {
	      up_sample_load_neg(re_zs,re_zs2,0);
	      n=up_num_zeros(&out_array[(UP_SAMPLE_SACRIFICE+ptr)*UP_SAMPLE_RATE],im_s_vec,POS,POS,
			     2,ptr,re_zs,re_zs2);
	    }
	  else
	    {
	      up_sample_load(&re_zs[ptr-UP_SAMPLE_SACRIFICE]);
	      n=up_num_zeros(&out_array[UP_SAMPLE_SACRIFICE*UP_SAMPLE_RATE],im_s_vec,POS,POS,
			     2,ptr,re_zs,re_zs2);
	    }
	  if(n==0)
	    {
	      printf("Failed to find zero near ptr=%d.\n",ptr);
	      //d_print_signs(2*UP_SAMPLE_RATE,&out_array[UP_SAMPLE_SACRIFICE*UP_SAMPLE_RATE],1);
	    }
	  else
	    count+=n;
	  ptr++;
	  continue;
	}
      if((this_sign==NEG)&&(y2>y1)&&(y2>y3)) // max below zero
	{
	  //printf("Negative stat point near %d\n",ptr);
	  stat_count++;
	  if(ptr-UP_SAMPLE_SACRIFICE<0)
	    {
	      up_sample_load_neg(re_zs,re_zs2,0);
	      n=up_num_zeros(&out_array[(UP_SAMPLE_SACRIFICE+ptr)*UP_SAMPLE_RATE],im_s_vec,NEG,NEG,
			     2,ptr,re_zs,re_zs2);
	    }
	  else
	    {
	      up_sample_load(&re_zs[ptr-UP_SAMPLE_SACRIFICE]);
	      n=up_num_zeros(&out_array[UP_SAMPLE_SACRIFICE*UP_SAMPLE_RATE],im_s_vec,NEG,NEG,
			     2,ptr,re_zs,re_zs2);
	    }
	  if(n==0)
	    {
	      printf("Failed to find zero near ptr=%d.\n",ptr);
	      //d_print_signs(2*UP_SAMPLE_RATE,&out_array[UP_SAMPLE_SACRIFICE*UP_SAMPLE_RATE],1);
	    }
	  else
	    count+=n;
	  ptr++;
	  continue;
	}
      ptr++;
    }
  return(count);
}

inline int num_zeros(int a, int b, int_double *re_zs, im_s *im_s_vec, int_double *re_zs2)
{
  return(num_zeros1(a,b,re_zs)+num_zeros2(b,re_zs,im_s_vec,re_zs2));
}

// compare number of zeros observed with that calculate by Turing's method
int check_zeros(int count_zeros, int t1_ptr, int t2_ptr, int q,
		int_double &gamma_int, int_double &arg_omega,int_double *re_zs, im_s *im_s_vec,
		int index, int num_s)
{
  double t1=im_s_vec[t1_ptr].im_s;
  double t2=im_s_vec[t2_ptr].im_s;
  int_double ntwiddle=up_n_twiddle(t1_ptr,t2_ptr,re_zs,im_s_vec,num_s);

  int_double s_t=s_of_t(q,t2);
  int_double lnt=ln_term(t1,t2,q);
  int_double calc_zeros;
  int calc_zeros_int;

  calc_zeros=(arg_omega+lnt+gamma_int/(t2-t1))/d_pi-ntwiddle/(t2-t1)+s_t/(t2-t1); // Turing's estimate
	
  calc_zeros_int=(int)ceil(calc_zeros.left);

  if(calc_zeros_int!=floor(-calc_zeros.right)) // Turing's estimate did not bracket a unique integer
    {
      printf("q: %ld index: %ld Problem with Turing's estimate:- ",q,index);
      print_int_double(calc_zeros);
      printf("\n");
      return(BAD_DIFF);
    }
  return(calc_zeros_int-count_zeros);
}

int find_more_zeros(int_double *re_zs,  int_double *re_zs2, int turing_start,  int missing, im_s *im_s_vec)
{
  int n,up_n;
  int t2=turing_start;
  int t1=t2-UP_SAMPLE_USABLE_WIDTH+1;
  int t3=t2+UP_SAMPLE_SACRIFICE;
  int t0=t3-UP_SAMPLE_IN_WIDTH+1;
  int still_missing=missing;
  int ptr_hi,ptr_lo;
  sign_t st_sign,en_sign;
  while(t0>0)
    {
      up_sample_load(&re_zs[t0]);
      while(sign(re_zs[t1])==UNK) t1--;
      en_sign=sign(re_zs[t2]);
      while(en_sign==UNK) en_sign=sign(re_zs[++t2]);
      ptr_hi=t2;
      ptr_lo=ptr_hi-1;
      while(ptr_hi>t1)
	{
	  st_sign=sign(re_zs[ptr_lo]);
	  while(st_sign==UNK)
	    st_sign=sign(re_zs[--ptr_lo]);
	  up_n=up_num_zeros(&out_array[(ptr_lo-t0)*UP_SAMPLE_RATE],im_s_vec,st_sign,en_sign,
			    ptr_hi-ptr_lo,ptr_lo,re_zs,re_zs2);
	  still_missing-=(up_n&0xFFFFFFFE); // if odd, will already know about one of them
	  if(still_missing==0)
	    return(0);
	  en_sign=st_sign;
	  ptr_hi=ptr_lo;
	  ptr_lo--;
	}
      t2=t1;
      t1=t2-UP_SAMPLE_USABLE_WIDTH+1;
      t3=t2+UP_SAMPLE_SACRIFICE;
      t0=t3-UP_SAMPLE_IN_WIDTH+1;
    }
  // if here then t0<0 so need to load the negative portion
  // we have checked from t2 upwards so now check
  // from 0 to t2-1
  // we have en_sign=sign of t2
  up_sample_load_neg(re_zs,re_zs2,t2);
  ptr_hi=t2;
  ptr_lo=ptr_hi-1;
  while(ptr_hi>0)
    {
      st_sign=sign(re_zs[ptr_lo]);
      while(st_sign==UNK)
	{
	  ptr_lo--;
	  if(ptr_lo<0)
	    {
	      printf("re_zs[0] is of unknown sign. Cannot resolve it.\n");
	      return(still_missing);
	    }
	  st_sign=sign(re_zs[ptr_lo]);
	}
      up_n=up_num_zeros(&out_array[(ptr_lo+UP_SAMPLE_SACRIFICE)*UP_SAMPLE_RATE],im_s_vec,
			st_sign,en_sign,ptr_hi-ptr_lo,ptr_lo,re_zs,re_zs2);
      // if up_n is odd, we will have counted one already
      still_missing-=(up_n&0xFFFFFFFE);
      
      if(still_missing==0)
	return(0);
      
      en_sign=st_sign;
      ptr_hi=ptr_lo;
      ptr_lo--;
    }
  return(still_missing);
}
	

inline int find_more_zeros2(int_double *re_zs1, int_double *re_zs2,  int turing_start,  int missing,
							  im_s *im_s_vec)
{
  int still_missing=missing;
  still_missing=find_more_zeros(re_zs1,re_zs2,turing_start,missing,im_s_vec);
  if(still_missing!=0)
    still_missing=find_more_zeros(re_zs2,re_zs1,turing_start,still_missing,im_s_vec);
  return(still_missing);
}

void save_state(int q, int index, bool real_p, bool neg_one, int_complex &omega, 
		bool neg_one2, int_complex & omega2, int num_zeros, int num_s, 
		int_double *re_zs1, int_double *re_zs2, FILE *out_file)
{
  fwrite(&q,sizeof(int),1,out_file);
  fwrite(&index,sizeof(int),1,out_file);
  fwrite(&real_p,sizeof(bool),1,out_file);
  fwrite(&num_zeros,sizeof(int),1,out_file);
  fwrite(&neg_one,sizeof(bool),1,out_file);
  fwrite(&omega,sizeof(int_complex),1,out_file);
  //fwrite(&num_s,sizeof(int),1,out_file);
  fwrite(re_zs1,sizeof(int_double),num_s,out_file);
  if(!real_p)
    {
      fwrite(&neg_one2,sizeof(bool),1,out_file);
      fwrite(&omega2,sizeof(int_complex),1,out_file);
      fwrite(re_zs2,sizeof(int_double),num_s,out_file);
    }
}



#define MAX_F_NAME (256)

 int main(int argc, char **argv)
 {
   FILE *infile,*out_file;
   im_s *im_s_vec;
   int i,j,num_s,index,index2,n_zeros;
   int q,turing_start,turing_end;
   int_double gamma,gamma_a,arg_omega,arg_omega2;
   int_double *re_zs1,*re_zs2;
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
   //fread(&gamma_a,sizeof(int_double),1,infile);
   //int_double foo;
   //fread(&foo,sizeof(int_double),1,infile);
   fread(&turing_start,sizeof(int),1,infile);
   fread(&turing_end,sizeof(int),1,infile);   
   fwrite(&num_s,sizeof(int),1,out_file);
   fwrite(&gamma,sizeof(int_double),1,out_file);
   fwrite(&gamma_a,sizeof(int_double),1,out_file);
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
       //printf("processing q=%d index %d\n",q,index);
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
	   //print_int_double_str("re_zs1[0]=",re_zs1[0]);
	   //print_int_double_str("re_zs2[0]=",re_zs2[0]);
	   //print_int_double_str("re_zs1[turing_start]=",re_zs1[turing_start]);
	   //print_int_double_str("re_zs2[turing_start]=",re_zs2[turing_start]);
	   //print_int_double_str("arg(omega1)=",arg_omega);
	   //print_int_double_str("arg(omega2)=",arg_omega2);
	 }
       if(real_p) // it's a real character
	 {
	   cz1=num_zeros(0,turing_start,re_zs1,im_s_vec,re_zs1);
	   no_zeros+=cz1;
	   diff1=check_zeros(cz1,turing_start,turing_end,q,(neg_one ? gamma_a : gamma),
			     arg_omega,re_zs1,im_s_vec,index,t_sam_start);
	   if(diff1==BAD_DIFF) // Error in Turing Zone
	     {
	       //printf("q:%d index:%d check_zeros (real) returned:%d\n",q,index,diff1);
	       save_state(q,index,true,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs1,out_file);
	       continue;
	     }
	   if(diff1!=0) // missed some zeros
	     {
	       no_zeros+=diff1;
	       save_state(q,index,true,neg_one,omega,neg_one2,omega2,cz1+diff1,num_s,re_zs1,re_zs1,out_file);
	       printf("Do some more checking on q:%d, index:%d, difference:%d\n",q,index,diff);
	     }
	 }
       else
	 {
	   //if(index!=31553)
	     //  continue;
	   cz1=num_zeros(0,turing_start,re_zs1,im_s_vec,re_zs2);
	   cz2=num_zeros(0,turing_start,re_zs2,im_s_vec,re_zs1);
	   //find_more_zeros(re_zs1,re_zs2,turing_start,diff1,im_s_vec);
	   //find_more_zeros(re_zs2,re_zs1,turing_start,diff1,im_s_vec);
	   //printf("Zeros found %d %d.\n",cz1,cz2);
	   no_zeros+=cz1+cz2;
	   diff1=check_zeros(cz1,turing_start,turing_end,q,(neg_one ? gamma_a : gamma),
			     arg_omega,re_zs1,im_s_vec,index,t_sam_start);   
	   if(diff1==BAD_DIFF) // problem in Turing method
	     {
	       printf("q:%d index:%d check_zeros (cmplx) returned:%d\n",q,index,diff1);
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs2,out_file);
	       continue;
	     }
	   diff=check_zeros(cz2,turing_start,turing_end,q,(neg_one2 ? gamma_a : gamma),
			    arg_omega2,re_zs2,im_s_vec,index2,t_sam_start);
	   if(diff==BAD_DIFF) // problem in Turing method with conjugate
	     {
	       printf("q:%d index:%d check_zeros (cmplx) returned:%d\n",q,index,diff1);
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs2,out_file);
	       continue;
	     }
	   diff1+=diff;
	   
	   if(diff1<0) // what th...
	     {
	       printf("q:%d index: %d %d too many zeros located.\n",q,index,index2);
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,0,num_s,re_zs1,re_zs2,out_file);
	     }
	   if(diff1>0) // missed some zeros
	     {
	       save_state(q,index,false,neg_one,omega,neg_one2,omega2,cz1+cz2+diff1,num_s,re_zs1,re_zs2,out_file);
	       printf("Do some more checking on q:%d, index:%d, difference:%d\n",q,index,diff1);
	       /*
	       print_int_complex_str("omega1=",omega);
	       print_int_complex_str("omega2=",omega2);
	       print_int_double_str("Arg(omega1)=",argument(omega));
	       print_int_double_str("Arg(omega2)=",argument(omega2));
	       print_int_double_str("Z1(1/2)=",re_zs1[0]);
	       print_int_double_str("Z2(1/2)=",re_zs2[0]);
	       */
	       no_zeros+=diff1;
	     }
	 }
       //printf("There were %d zeros for q=%d and index=%d T<=%f.\n",no_zeros,q,index,im_s_vec[turing_start].im_s);
       //printf("Stat_count was %d.\n",stat_count); 
       stat_count=0;
       no_zeros=0;
       //printf("Check_sign was called %d times.\n",check_count);
       check_count=0;
     }
 }
