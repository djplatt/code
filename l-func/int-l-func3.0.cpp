/*
File: int-l-func3.0.cpp

Created: 6th Jan 2009

Version: <v> = 3.0

Last Modified: 26th Jan 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: Takes Output from l-func-mpfi-1.0
                      num_s (int)
					  N (int)
					  rn (int)
					  NO_GAPS (int)

					  im_s (double)
					  gam_s (dcomplex)
					  gam_s_2 (dcomplex)
					  N*(NO_GAPS+1) h_rn (dcomplex)

Build instructions: 

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "3.0"
#define ERRORS
//#define DEBUGGING



#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;



#include "../includes/int_double10.0.h"
typedef int_complex dcomplex;

void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: int-l-func%s (q-start) (q-end) (ifname) (ofname) (N)\n",VERSION);
	printf("  (q-start)   - integer >=3\n");
	printf("  (q-end)     - integer >=(q-start)\n");
	printf("  (ifname)    - file with lattice values.\n");
	printf("  (ofname)    - output file.\n");
	printf("  (N)         - number of Taylor terms to use.\n");
	exit(0);
}


void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(0);
};

inline int gcd (unsigned int a, unsigned int b)
  /* Euclid algorithm gcd */
{
  unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

// make a the smaller
inline int co_prime(unsigned int a, unsigned int b)
{
  return(gcd(a,b)==1);
};

// only a vector here. last element used for taylor error
void create_s_array(dcomplex &s,dcomplex *s_array, unsigned int n, double &gap)
{
  unsigned int i;
  s_array[0].real=int_double(gap/2); // 0.5*gap is representable
  s_array[0].imag=s.imag*gap;
  for(i=1;i<n;i++)
	  s_array[i]=s_array[i-1]*(s+i)*gap/(i+1); // ditto and (i+1)
}

void print_r_n_vec(dcomplex *r_n_vec,unsigned int n)
{
  unsigned int i;

  for(i=0;i<n;i++)
    {
      cout << "r_n_vec[" << i << "] is ";
      print_int_complex(r_n_vec[i]);
      cout << endl;
    };
}

// sum (alpha+n)^(-s) from n=R to n=RN
dcomplex zeta (dcomplex &s, int_double &alpha, unsigned int start, unsigned int end)
{
	dcomplex tot=c_zero;
	int i,en=end,st=start;
	for(i=en;i>=st;i--)
		tot=tot+pow(alpha+(double)i,-s);
	return(tot);
}


dcomplex calc_r_n(unsigned int r_n_terms, dcomplex *s_array, int_double &delta,
		  dcomplex *r_n_vec,int_double &alpha)
{
  unsigned int i;
  dcomplex res,term;
  int_double delta_n=delta;

  res=r_n_vec[0];


  for(i=0;;)
    {
      term=s_array[i]*delta_n;
	  
	  i++;
      term=term*r_n_vec[i];
      res=res+term;
	  if(i==r_n_terms-1)
		  break;
	  else
		  delta_n=times_pos(delta_n,delta);
    }
  return(res);
}

dcomplex calc_r_n_neg(unsigned int r_n_terms, dcomplex *s_array, int_double &delta,
		  dcomplex *r_n_vec,int_double &alpha)
{
  unsigned int i;
//  int_double ln_delta=log(-delta);
  int_double delta_n=delta; // will be negative
//  int sign=-1;
  dcomplex res,term;

  res=r_n_vec[0];


  for(i=0;i<r_n_terms-1;)
    {
      term=s_array[i]*delta_n;
	  delta_n=delta_n*delta; // can't use quick multiply
      i++;
      term=term*r_n_vec[i];
	  res=res+term;
    }
  return(res);
}

#define PRINT_ERROR_STEP 16
void print_errors(int_complex &s,int_complex *r_n_vals, unsigned int N, unsigned int no_gaps,unsigned int M)
{
	unsigned int i,j;
	int err;
	printf("Rel Error in r_n_vals at Im(s)=");
	print_int_double(s.imag);
	printf("\n");
	for(i=0;i<=PRINT_ERROR_STEP;i++)
	{
		printf("At%2d/%2d ",i,PRINT_ERROR_STEP);
		for(j=0;j<N;j++)
		{
			if(j==M)
				printf("|");
			err=rel_error(r_n_vals[N*i*no_gaps/PRINT_ERROR_STEP+j].real);
			printf("%3d",err);
		}
		printf("\n");
	}
	
}

int_double read_int_double(FILE *in_file)
{
	double x[2];
	int_double ans;
	fread(&x,sizeof(double),2,in_file);
	ans.left=x[0];
	ans.right=-x[1];
	return(ans);
}


dcomplex read_dcomplex(FILE *in_file)
{
	double x[4];
	dcomplex ans;
	fread(x,sizeof(double),4,in_file);
	ans.real.left=x[0];
	ans.real.right=-x[1];
	ans.imag.left=x[2];
	ans.imag.right=-x[3];
	return(ans);
}

void skip_dcomplex(FILE *in_file,unsigned int n)
{
	unsigned int i;
	double x[4];
	for(i=0;i<n;i++)
		fread(x,sizeof(double),4,in_file);
}

#ifdef ERRORS
dcomplex taylor_error(dcomplex &s, unsigned int N,unsigned int RN,double gap,dcomplex *s_array)
{
	int_double a,r,err;

	// set a to max mod of H_RN(s+N,alpha) with alpha in [0,1]
	if(RN>1)
	{
		printf("Need to revisit error estimation for RN>1 case. Exiting.\n");
		exit(0);
		a=int_double(pow(RN-1,-(double)N-0.5))/(N-0.5);
	}
	else // RN=1
	{
		a=int_double((double)N+0.5);
		a=a/((double)N-0.5);
		// now multiply by s(s+1)..(s+N)*(gap/2)^N/N! assuming a/q is exactly half way
		a=a*sqrt(norm(s_array[N-1]))/(1<<N);
		r=sqrt(norm(s+N))*gap;
		r=r/(N+N+2);
	}	
	err=a/(r-1.0); // this is -ve what we want but never mind
	if(err.right>=0.0)  // interval in (-infty,0.0]
		err.right=err.left;
	else
	{
		if(err.left<=0.0) // interval in [0.0,infty)
			err.left=err.right;
		else
		{
			if(err.right>=err.left)
				err.left=-err.right;
			else
				err.right=-err.left;
		}
	}
	return(dcomplex(err,err));
}
#endif


int calc_hurwitz1(dcomplex &s,unsigned int r_n_terms,unsigned int file_N,double gap,dcomplex *s_array,
				  dcomplex *r_n_vals,unsigned int no_gaps,unsigned int q_start,
				  unsigned int q_end, FILE *in_file, FILE *out_file,unsigned int RN//,unsigned int counter
				  )
{
  dcomplex out_val;
#ifdef ERRORS
  dcomplex taylor_err;
#endif
  unsigned int i,j,gap_ptr,q,num,steps,num_fracs;
  //int r_err,i_err,diff_err;
  int_double x,frac;
  dcomplex q_minus_s_2,minus_s;
//  int rel_err;
  double dsteps,csteps;

  minus_s=-s;

  create_s_array(s,s_array,r_n_terms,gap); // end column used for error estimation
#ifdef ERRORS
  taylor_err=taylor_error(s,r_n_terms,RN,gap,s_array);
#endif
  gap_ptr=0;
  
  for(i=0;i<=no_gaps;i++)
  {

	  for(j=0;j<r_n_terms;j++)
		  r_n_vals[gap_ptr++]=read_dcomplex(in_file);
	  skip_dcomplex(in_file,file_N-r_n_terms);
  }
//#define TEST_LATTICE
#ifdef TEST_LATTICE
  FILE *zeta_vals;
  double re_z,im_z;
  if(!(zeta_vals=fopen("zeta_495_2990_3010.dat","r"))) // this file contains 250 rows of 15 values of zeta
	  exit(0);
  for(i=2990;i<=3010;i++)
	  for(j=0;j<r_n_terms;j++)
	  {
		  fscanf(zeta_vals,"%le",&re_z);
		  fscanf(zeta_vals,"%le",&im_z);
		  if((re_z!=r_n_vals[i*r_n_terms+j].real.left)&&(re_z!=-r_n_vals[i*r_n_terms+j].real.right))
			  printf("Re Mismatch at Row %d Column %d.\n",i,j);
		  if((im_z!=r_n_vals[i*r_n_terms+j].imag.left)&&(im_z!=-r_n_vals[i*r_n_terms+j].imag.right))
			  printf("Im Mismatch at Row %d Column %d.\n",i,j);
	  }
	  fclose(zeta_vals);
	  exit(0);
#endif  

#ifdef DEBUGGING
  print_errors(s,r_n_vals,r_n_terms,no_gaps,r_n_terms);
//  print_int_complex(s);printf(" ");print_int_complex(r_n_vals[r_n_terms*no_gaps/2]);printf("\n");
#endif
  num_fracs=0;
  for(q=q_start;q<=q_end;q++)
    {
		if((q&3)==2)
			continue;      // if q = 2 mod 4 then no primitive characters
		q_minus_s_2=pow(q,-s/2);
		for(num=1;num<q;num++)
		  if(co_prime(num,q))
			  {
				num_fracs++;
				frac=int_double(num)/q;
				dsteps=double(num)/(gap*q);
				csteps=ceil(dsteps);
				steps=(int) csteps;
				if((csteps-dsteps)>=0.5) // use lower row of R_ns
				{
					gap_ptr=r_n_terms*(no_gaps-steps);
					x=int_double(steps)-frac/gap;
					out_val=calc_r_n(r_n_terms,s_array,x,&r_n_vals[gap_ptr],frac);
				}
				else // use upper row
				{
					gap_ptr=r_n_terms*(no_gaps-steps+1);
					x=int_double(steps)-(frac/gap+1.0);
					out_val=calc_r_n_neg(r_n_terms,s_array,x,&r_n_vals[gap_ptr],frac);
				}
				// out_val is now ~ sum_{n=RN}^\infty (n+num/q)^(-s)
#ifdef ERRORS				
				out_val=out_val+taylor_err;
#endif
				out_val=out_val+pow(frac,minus_s);
				for(i=1;i<RN;i++)
					out_val=out_val+pow(frac+i,minus_s);
				out_val=out_val*q_minus_s_2; // out_val now contains Zeta(s,num/q)*q^(-s/2)
				// uncomment to reduce im_s by factor of two
//				if((counter&1)==0) // every other s value
				fwrite(&out_val,sizeof(dcomplex),1,out_file);
		  }
  }
  return(num_fracs);
}

int main(int argc, char **argv)
{
  _fpu_rndd();

  dcomplex s,gam,*s_array,*r_n_vals;
  int i,no_gaps,q_start,q_end,num_s,num_fracs,N,rn,file_N;
  int_double mgap,im_s_2_mod;
  double gap,im_s;
  FILE *in_file,*out_file;
  fpos_t out_file_ptr;
  int file_pos;
  clock_t no_clicks;

  no_clicks=clock(); // start timing

  if(argc!=6)
    print_usage();
  q_start=atoi(argv[1]);
  if(q_start<3)
    print_usage();
  q_end=atoi(argv[2]);
  if(q_end<q_start)
    print_usage();
  in_file=fopen(argv[3],"rb");
  if(!in_file)
    fatal_error("Couldn't open in_file. Exiting.\n");
  N=atoi(argv[5]);
  fread(&num_s,sizeof(int),1,in_file);
  fread(&file_N,sizeof(int),1,in_file);
  if((N>file_N)||(N<=0))
	  fatal_error("N<=0 or exceeds N used to save lattice file. Exiting.\n");
  fread(&rn,sizeof(int),1,in_file);
  if(rn<1)
	  fatal_error("RN must be at least 1. Exiting.\n");
  out_file=fopen(argv[4],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");
  fread(&no_gaps,sizeof(int),1,in_file);
  gap=1.0/no_gaps;
  cout << "Lattice is Zeta(1/2+it+n,alpha+"<< rn <<") for n=0.." << N-1 << ", alpha in [0,1] in steps of 1/" << no_gaps << endl;
  cout << "Q is running from " << q_start << " to " << q_end << endl;


  cout << "Allocating memory for s_array." << endl;
  if(!(s_array=(dcomplex *) _aligned_malloc(sizeof(dcomplex)*N,16)))
    fatal_error("Can't allocate memory for s_array. Exiting.");
  
  
  cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=(dcomplex *) _aligned_malloc(sizeof(dcomplex)*(no_gaps+1)*N,16)))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");


  fwrite(&q_start,sizeof(int),1,out_file);
  fwrite(&q_end,sizeof(int),1,out_file);
  fgetpos(out_file,&out_file_ptr);
  fwrite(&num_fracs,sizeof(int),1,out_file); // garbage at this point, will be overwritten
  // uncomment to reduce im_s by factor of two
  //i=num_s/2;
  fwrite(&num_s,sizeof(int),1,out_file);
  printf("Going to try to process %d values of s.\n",num_s);

  s.real=int_double(0.5);

  file_pos=ftell(in_file);  // remember where in file we are

  for(i=0;i<num_s;i++)
  {
		fread(&im_s,sizeof(double),1,in_file);
		// uncomment to reduce im_s by factor of two
//		if((i&1)==0)
		fwrite(&im_s,sizeof(double),1,out_file);
		gam=read_dcomplex(in_file);
		// uncomment to reduce im_s by factor of two
//		if((i&1)==0)
		fwrite(&gam,sizeof(dcomplex),1,out_file);
		gam=read_dcomplex(in_file);
		// uncomment to reduce im_s by factor of two
//		if((i&1)==0)
		fwrite(&gam,sizeof(dcomplex),1,out_file);
		fseek(in_file,ftell(in_file)+(no_gaps+1)*N*sizeof(dcomplex),SEEK_SET);
  }

  fseek(in_file,file_pos,SEEK_SET);

  for(i=0;i<num_s;i++)
    {
		fread(&im_s,sizeof(double),1,in_file);
		s.imag=int_double(im_s);
		gam=read_dcomplex(in_file);
		gam=read_dcomplex(in_file);
//		printf("Processing s=0.5+i%10.8e.\n",im_s);
		num_fracs=calc_hurwitz1(s,N,file_N,gap,s_array,r_n_vals,no_gaps,q_start,q_end,in_file,out_file,rn/*,i*/);
		if(!num_fracs)
			fatal_error("Error running Hurwitz routine. Exiting.");
    }

  
  fsetpos(out_file,&out_file_ptr);
  fwrite(&num_fracs,sizeof(int),1,out_file);
  fclose(in_file);
  fclose(out_file);

  printf("Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  printf("Fractions processed=%d.\n",num_fracs);
  return(0);
}


