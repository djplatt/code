/*

File: int-l-func2.0.cpp

Created: 23rd August 2008

Version: <v> = 2.0

Last Modified: 5th November 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: 

Build instructions: 

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "2.0"
//#define ERRORS
#define DEBUGGING



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;


#define INT_DOUBLE

#ifdef SSE
#include "IntervalComplex.h"  // the complex interval class
typedef MachineEstimate int_double;

#else

#ifdef INT_DOUBLE
#include "int_double6.0.h"
typedef int_complex dcomplex;
#define nextafter(x,y) nextafter(x)

#else
// this option isn't working yet
#include <complex>
using namespace std;
typedef double int_double;
typedef complex<double> dcomplex;

#endif
#endif



#define debug printf("got to line %d\n",__LINE__);

unsigned int RN; // lattice will contain R_{RN} values ie sum n=RN..infty (n+alpha)^(-s)

void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: int-l-func%s (q-start) (q-end) (R terms) (N) (M)\n",VERSION);
	printf("                    (ifname) (ofname) (gap) (zeta_terms)\n");
	printf("  (q-start)   - integer >=3\n");
	printf("  (q-end)     - integer >=(q-start)\n");
	printf("  (R terms)   - integer > 0, # of terms to use for lattice points.\n");
	printf("  (N)         - integer >0, <= (R terms), # of terms to use\n");
	printf("                                          for intermediate values.\n");
	printf("  (M)         - integer >0, < (R terms), # of terms to calculate\n");
	printf("                                         by Taylor expansion.\n");
	printf("  (ifname)    - file with zeta values.\n");
	printf("  (ofname)    - output file.\n");
	printf("  (gap)       - 0.0 < gap < 1.0, spacing of lattice points.\n");
	printf("  (zeta_terms)- integer >0, number of terms to use in zeta approximation.\n");
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

inline int co_prime(unsigned int a, unsigned int b)
{
  return(gcd(a,b)==1);
};

void create_s_array(dcomplex &s,dcomplex *s_array, unsigned int n, double gap)
{
  unsigned int i,del;
  double fact=2.0;

  s_array[0]=s*gap;
  for(i=1;i<n;i++)
      s_array[i*n+i]=dcomplex(real(s)+i,imag(s))*gap;

  del=1;
  while(del<n)
    {
      for(i=0;i<n-del;i++)
		  s_array[i*n+i+del]=s_array[i*n+i]*s_array[i*n+n+i+del]/fact;
      del++;
      fact=fact+1.0;
    };
};

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
	dcomplex tot=c_zero,p;
	unsigned int i;
	for(i=start;i>=end;i--)
		tot+=pow(alpha+(double)i,-s);
	return(tot);
}

// estimate the error truncating R1(s,alpha) to M terms
// power = Re(s)-1
// M_alpha = M+alpha
dcomplex zeta_error(double power, int_double &M_alpha)
{
	int_double temp;
	dcomplex c_temp;

	temp=pow(M_alpha,-power)/power;
	temp.left=temp.right;
	c_temp.imag=temp;
	c_temp.real=temp;
	return(c_temp);
};

// 
int_double R1_est(double re_s, int_double &alpha)
{
	return((alpha+re_s)/(pow(alpha+1,re_s)*(re_s-1)));
}

dcomplex rn_error(double re_s,unsigned int N, dcomplex &s_array, int_double &alpha, int_double &delta_ratio)
{
	dcomplex res;
	int_double temp;

	temp=sqrt(norm(s_array));
	temp*=pow(delta_ratio,N);
	temp*=R1_est(re_s+N,alpha);
	temp.left=temp.right;
	res.real=temp;
	res.imag=temp;
	return(res);
}

void make_r_ns (dcomplex *s_array,
		unsigned int n, // number of terms in lattice
		unsigned int M, // number to generate by Taylor
		dcomplex *this_r_n_vec,
		dcomplex *last_r_n_vec,
		dcomplex &s,
		int_double &alpha,
		dcomplex *zeta_terms,
		unsigned int R) // upper limit of hurwitz zeta sum
{
  unsigned int i;
  int j;
  int_double err;
  dcomplex s_m;

  for(i=0;i<M;i++)
    {
		this_r_n_vec[i]=c_zero;
		for(j=n-i-2;j>=0;j--)
		{
			this_r_n_vec[i]+=s_array[i*(n+1)+j]*last_r_n_vec[j+i+1];
			
		}
		this_r_n_vec[i]+=last_r_n_vec[i];

#ifdef ERRORS	  
	  err=sqrt(norm(s_array[(i+1)*n-1]))*R1_est(0.5+n,alpha);
	  err.left=err.right;
	  this_r_n_vec[i].real+=err;
	  this_r_n_vec[i].imag+=err;
#endif
    }
/* // this was how I did it before RN came in
  s_m=s+M;
  for(i=0;i<R;i++)
	  zeta_terms[i]=pow(alpha+(i+1),-s_m);
  this_r_n_vec[M]=zeta_terms[0];
	for(i=1;i<R;i++)
		this_r_n_vec[M]+=zeta_terms[i];

	for(i=M+1;i<n;i++)
	{
		for(j=0;j<(int)R;j++)
			zeta_terms[j]=zeta_terms[j]/(alpha+(j+1.0));
		this_r_n_vec[i]=zeta_terms[0];
		for(j=1;j<(int)R;j++)
			this_r_n_vec[i]+=zeta_terms[j];
	}
	*/
  for(i=M;i<n;i++)
	  this_r_n_vec[i]=zeta(s+i,alpha,RN,R);
#ifdef ERRORS
	for(i=M;i<n;i++) // from s=0.5+M-1 ...
		this_r_n_vec[i]+=zeta_error(i-0.5,alpha+R);
#endif
}

void make_r_ns_neg (dcomplex *s_array,
		unsigned int n,
		unsigned int M,
		dcomplex *this_r_n_vec,
		dcomplex *last_r_n_vec,
		dcomplex &s,
		int_double &alpha,
		dcomplex *zeta_terms,
		unsigned int R)
{
  unsigned int i;
  int j;
  int sign;
  dcomplex s_m;
  int_double err;

 
  for(i=0;i<M;i++)
    {
		this_r_n_vec[i]=c_zero;
		if((n-i)&1)
			sign=1;
		else
			sign=-1;
		for(j=n-i-2;j>=0;j--)
		{
		  if(sign==1)
			  this_r_n_vec[i]+=s_array[i*(n+1)+j]*last_r_n_vec[j+i+1];
		  else
			  this_r_n_vec[i]-=s_array[i*(n+1)+j]*last_r_n_vec[j+i+1];
		  sign=-sign;
		}
		this_r_n_vec[i]+=last_r_n_vec[i];

#ifdef ERRORS
	  
	  err=sqrt(norm(s_array[(i+1)*n-1]))*R1_est(0.5+n,alpha);
	  err.left=err.right;
	  this_r_n_vec[i].real+=err;
	  this_r_n_vec[i].imag+=err;
#endif
  }
/*
  s_m=s+M;
  for(i=0;i<R;i++)
	  zeta_terms[i]=pow(alpha+(i+1),-s_m);
  this_r_n_vec[M]=zeta_terms[0];
	for(i=1;i<R;i++)
		this_r_n_vec[M]+=zeta_terms[i];

	for(i=M+1;i<n;i++)
	{
		for(j=0;j<(int)R;j++)
			zeta_terms[j]=zeta_terms[j]/(alpha+(j+1.0));
		this_r_n_vec[i]=zeta_terms[0];
		for(j=1;j<(int)R;j++)
			this_r_n_vec[i]+=zeta_terms[j];
	}*/

  for(i=M;i<n;i++)
	  this_r_n_vec[i]=zeta(s+i,alpha,RN,R);
#ifdef ERRORS
		for(i=M;i<n;i++)
	    	this_r_n_vec[i]+=zeta_error(i-0.5,alpha+R);
#endif

}



dcomplex calc_r_n(unsigned int r_n_terms, dcomplex *s_array, int_double &delta,
		  dcomplex *r_n_vec,int_double &alpha)
{
  unsigned int i;
  dcomplex res,term;
  int_double ln_delta=log(delta);

  res=r_n_vec[0];


  for(i=0;i<r_n_terms-1;)
    {
      term=s_array[i]*exp(ln_delta*(i+1));
      i++;
      term*=r_n_vec[i];
      res+=term;
    };
#ifdef ERRORS
  res+=rn_error(0.5,r_n_terms,s_array[r_n_terms],alpha,delta);
#endif
  return(res);
};

dcomplex calc_r_n_neg(unsigned int r_n_terms, dcomplex *s_array, int_double &delta,
		  dcomplex *r_n_vec,int_double &alpha)
{
  unsigned int i;
  int_double ln_delta=log(-delta);
  int sign=-1;
  dcomplex res,term;

  res=r_n_vec[0];


  for(i=0;i<r_n_terms-1;)
    {
      term=s_array[i]*exp(ln_delta*(i+1));
      i++;
      term*=r_n_vec[i];
	  if(sign==1)
		  res+=term;
	  else
		  res-=term;
	  sign=-sign;
    }
#ifdef ERRORS
  res+=rn_error(0.5+r_n_terms,r_n_terms,s_array[r_n_terms],alpha,-delta);
#endif

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

void print_s_array_errors(int_complex &s,int_complex *s_array,int r_n_terms)
{
	int i,j,ptr=0;
	int err;
	printf("Im(s)=");print_int_double(s.imag);printf("\n");
	for(i=0;i<r_n_terms;i++)
		for(j=0;j<r_n_terms;j++)
		{
			if((ptr%r_n_terms)==0)
				printf("\n");
			if(j<i)
				printf("   ");
			else
			{
				err=rel_error(s_array[ptr].real);
				printf("%3d",err);
			}
			ptr++;
		}
		printf("\n");
}

int calc_hurwitz1(dcomplex &s, unsigned int r_n_terms,
				  unsigned int N, unsigned int M,
				  double gap,dcomplex *s_array,
				  dcomplex *r_n_vals,unsigned int no_gaps,
				  unsigned int q_start,
				  unsigned int q_end, FILE *out_file, dcomplex *zeta_terms, unsigned int R)
{
  dcomplex out_val;
  unsigned int i,gap_ptr,q,num,steps,num_fracs;
  int_double x,frac;
  dcomplex q_minus_s_2;
  double dsteps,csteps;

  create_s_array(s,s_array,r_n_terms,gap); // end column used for error estimation

#ifdef DEBUGGING
  // there is no problem with s_array accuracy, at worst 10^{-13}
  // print_s_array_errors(s,s_array,r_n_terms);
#endif

  num=no_gaps/2;
  for(i=1;i<=num;i++)
    make_r_ns(s_array,r_n_terms,M,&r_n_vals[i*r_n_terms],
	      &r_n_vals[(i-1)*r_n_terms],s,int_double(1.0-gap*i),zeta_terms,R);
  for(i=no_gaps-1;i>num;i--)
    make_r_ns_neg(s_array,r_n_terms,M,&r_n_vals[i*r_n_terms],
	      &r_n_vals[(i+1)*r_n_terms],s,int_double(1.0-gap*i),zeta_terms,R);
		  


#ifdef DEBUGGING
  print_errors(s,r_n_vals,r_n_terms,no_gaps,M);
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
					out_val=calc_r_n(N,s_array,x,&r_n_vals[gap_ptr],frac);
				}
				else // use upper row
				{
					gap_ptr=r_n_terms*(no_gaps-steps+1);
					x=int_double(steps)-(frac/gap+1.0);
					out_val=calc_r_n_neg(N,s_array,x,&r_n_vals[gap_ptr],frac);
				}
				// out_val is now ~ sum_{n=RN}^\infty (n+num/q)^(-s)
				
				out_val+=zeta(s,frac,0,RN-1);
				out_val*=q_minus_s_2;

				//fwrite(&out_val,sizeof(dcomplex),1,out_file);
				/*
				if((num==499))
				{
					print_int_complex(s);
					printf(" ");
					print_int_complex(out_val);
					printf("\n");
				}
				*/

		  }
  }
  return(num_fracs);
}

int main(int argc, char **argv)
{
#ifdef SSE
  MachineEstimate::Computation comp;
#endif
#ifdef INT_DOUBLE
_fpu_rndd;
#endif

double ds,ds1;

  dcomplex s,*s_array,*r_n_vals,*zeta_terms/* ,*n_s */;



  int r_n_terms,i,j,no_gaps,q_start,q_end,num_s,num_fracs,num_zeta,N,M,R;
  int_double mgap;
  double gap,s_delta,im_s;
  ifstream zeta_file;   // contains zeta(s)-1 for s=0.5,1.5....99.5,0.5+dsI....,0.5+2dsI.....
  FILE *out_file;
  fpos_t out_file_ptr;

  clock_t no_clicks;

  no_clicks=clock();

  if(argc!=10)
    print_usage();
  q_start=atoi(argv[1]);
  if(q_start<3)
    print_usage();
  q_end=atoi(argv[2]);
  if(q_end<q_start)
    print_usage();
  r_n_terms=atoi(argv[3]);
  if(r_n_terms<=0)
    print_usage();
  N=atoi(argv[4]);
  if((N<=0)||(N>r_n_terms))
	  print_usage();
  M=atoi(argv[5]);
  if((M<=0)||(M>=r_n_terms))
	  print_usage();
  zeta_file.open(argv[6]);
  zeta_file >> RN;
  if(!zeta_file.is_open())
    fatal_error("Couldn't open ifile. Exiting.\n");
  out_file=fopen(argv[7],"wb");
  if(!out_file)
    fatal_error("Couldn't open ifile. Exiting.\n");
  gap=atof(argv[8]);
  if((gap<=0.0)||(gap>=1.0))
    print_usage();
  // number of terms to use in sum from 1 (n+alpha)^(-s)
  R=atoi(argv[9]);
  if(R<=0)
	  print_usage();

  no_gaps=(int) ceil(1.0/gap);
  cout << "Running with " << no_gaps << " gaps." << endl;
  cout << "Lattice is Hur(s,alpha) less first " << RN << " terms." << endl;
  cout << "         and " << r_n_terms << " r_n_terms." << endl;
  cout << "    of which " << M << " are by Taylor." << endl;
  cout << "         and " << r_n_terms-M << " are by summing." << endl;
  cout << "       Using " << N << " terms for intermediate values." << endl;
  cout << "   and using " << R << " terms to sum." << endl;


  cout << "Allocating memory for s_array." << endl;
  if(!(s_array=(dcomplex *) _aligned_malloc(sizeof(dcomplex)*(r_n_terms)*(r_n_terms),16)))
    fatal_error("Can't allocate memory for s_array. Exiting.");
  
  
  cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=(dcomplex *) _aligned_malloc(sizeof(dcomplex)*(no_gaps+1)*r_n_terms,16)))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");


  cout << "Allocating memory for zeta_terms.\n" << endl;
  if(!(zeta_terms=(dcomplex *) _aligned_malloc(sizeof(dcomplex)*R,16)))
	  fatal_error("Can't allocate memory for zeta_terms. Exiting.");


  fwrite(&q_start,sizeof(int),1,out_file);
  fwrite(&q_end,sizeof(int),1,out_file);
  fgetpos(out_file,&out_file_ptr);
  fwrite(&num_fracs,sizeof(int),1,out_file); // garbage at this point, will be overwritten

  zeta_file >> num_zeta; // number of zeta terms per s
  if(num_zeta<r_n_terms)
    fatal_error("Not enough zeta values in zeta_file. Exiting.\n");
  zeta_file >> num_s;    // number of s's
  fwrite(&num_s,sizeof(int),1,out_file);
  zeta_file >> s_delta;      // step between each s.
                             // assumes delta is representable

  printf("Going to try and process %d values of s using step size %20.18e\n    ",
	 num_s,s_delta);

  im_s=0.0;
  s=int_complex(d_half,d_zero);
  for(i=0;i<num_s;i++)
  {
	  fwrite(&im_s,sizeof(double),1,out_file);
	  im_s+=s_delta;
  }
  im_s=0.0;
  for(i=0;i<num_s;i++)
    {
//		printf("s:0.5+i%10.8e \n",im_s);
		zeta_file >> ds >> ds1; // gamma(s/2) real and imag
		  zeta_file >> ds >> ds1; // gamma((s+1)/2) real and imag
		for(j=0;j<r_n_terms;j++)
	  {
		  zeta_file >> ds >> ds1;
		  r_n_vals[j]=dcomplex(int_double(ds,nextafter(ds,DBL_MAX)),
			  int_double(ds1,nextafter(ds1,DBL_MAX)));
		  r_n_vals[j+no_gaps*r_n_terms]=r_n_vals[j]+pow(RN,-s);
	  };
	  for(;j<num_zeta;j++)
		  zeta_file >> ds >> ds;    // throw away extra ones

	  num_fracs=calc_hurwitz1(s,r_n_terms,N,M,gap,s_array,
		  r_n_vals,no_gaps,q_start,q_end,out_file,zeta_terms,
		  R);
	  if(!num_fracs)
		  fatal_error("Error running Hurwitz routine. Exiting.");
      im_s+=s_delta;
	  s.imag=im_s;
    };

  
  fsetpos(out_file,&out_file_ptr);
  fwrite(&num_fracs,sizeof(int),1,out_file);
  zeta_file.close();
  fclose(out_file);

  printf("Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  return(0);
};


