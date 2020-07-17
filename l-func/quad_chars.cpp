/*
File: quad_chars.cpp

Created: 24th June 2009

Version: <v> = 1.0

Last Modified: 14th June 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: Takes Output from l-func-mpfi-1.0
                      num_s (int)
					  N (int)
					  rn (int)
					  NO_GAPS (int)

					  im_s (double)
					  gam_s (double)
					  gam_s_2 (double)
					  N*(NO_GAPS+1) h_rn (double)

Build instructions: 

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "1.0"
#define ERRORS
//#define DEBUGGING



#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <math.h>

using namespace std;

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

inline long long int gcd (long long int a, long long int b)
  /* Euclid algorithm gcd */
{
  long long int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};


// make a the smaller
inline int co_prime(long long int a, long long int b)
{
  return(gcd(a,b)==1);
};

// implement me!
inline int kronecker (long long int a, long long int n)
{
// we are assured that n>a>0, gcd(a,n)=1, n odd
	long long int m8,tmp;
	int k=1;
	bool odd_pow;
	while(true)
	{
		if(a==1)
			return(k);
		if((a&1)==0) // even
		{
			m8=(n&7);
			if((m8==1)||(m8==7))
			{
				a>>=1;
				while((a&1)==0)
					a>>=1;
			}
			else
			{
				a>>=1;
				odd_pow=true;
				while((a&1)==0)
				{
					a>>=1;
					odd_pow=!(odd_pow);
				}
				if(odd_pow)
					k=-k;
			}
		}
		if(a==1)
			return(k);
		tmp=a;
		if(((a&3)==1)||((n&3)==1))
		{
			a=n%a;
			n=tmp;
		}
		else
		{
			a=n%a;
			n=tmp;
			k=-k;
		}
	}
}

// only a vector here. last element used for taylor error
void create_s_array(double s, double *s_array, unsigned int n, double &gap)
{
  unsigned int i;
  s_array[0]=gap/2.0;
  for(i=1;i<n;i++)
	  s_array[i]=s_array[i-1]*(s+i)*gap/(i+1); // ditto and (i+1)
}

void print_r_n_vec(double *r_n_vec,unsigned int n)
{
  unsigned int i;

  for(i=0;i<n;i++)
    {
      cout << "r_n_vec[" << i << "] is ";
      printf("%10.8e",r_n_vec[i]);
      cout << endl;
    };
}

double calc_r_n(unsigned int r_n_terms, double *s_array, double &delta,
		  double *r_n_vec,double &alpha)
{
  unsigned int i;
  double res;
  /*
  double term,delta_n=delta;

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
		  delta_n*=delta;
    }
	*/
  
  // add smallest terms first. May bulge in the middle but...
  res=s_array[r_n_terms-2]*r_n_vec[r_n_terms-1]*delta;
  for(i=r_n_terms-2;i>0;i--)
	  res=(res+r_n_vec[i]*s_array[i-1])*delta;

  return(res+r_n_vec[0]);
  
  // return(res);
}

double calc_r_n_neg(unsigned int r_n_terms, double *s_array, double &delta,
		  double *r_n_vec,double &alpha)
{
  unsigned int i;
  double delta_n=delta; // will be negative
  double res,term;

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

inline void skip_double(FILE *infile)
{
	double d;
	fread(&d,sizeof(double),1,infile);
}


inline void skip_int_double(FILE *infile)
{
	double d[2];
	fread(d,sizeof(double),2,infile);
}

inline void skip_dcomplex(FILE *infile,unsigned int n)
{
	double d[4];
	unsigned int i;
	for(i=0;i<n;i++)
		fread(d,sizeof(double),4,infile);
}



int calc_hurwitz1(double &s,unsigned int r_n_terms,unsigned int file_N,double gap,double *s_array,
				  double *r_n_vals,unsigned int no_gaps,long long int q_start,
				  long long int q_end, FILE *in_file, FILE *out_file,unsigned int RN
				  )
{
  double out_val,res;
  unsigned int i,j,gap_ptr,steps;
  long long int q, num;
  double x,frac;
  double q_minus_s,minus_s;
  double dsteps,csteps;

  minus_s=-s;

  create_s_array(s,s_array,r_n_terms,gap); // end column used for error estimation
  gap_ptr=0;
  
  for(i=0;i<=no_gaps;i++)
  {

	  for(j=0;j<r_n_terms;j++)
	  {
  		  fread(&r_n_vals[gap_ptr++],sizeof(double),1,in_file); // left of real
		  skip_double(in_file); // right of real
		  skip_int_double(in_file); // imag
	  }
		  skip_dcomplex(in_file,file_N-r_n_terms);
  }

  for(q=q_start;q<=q_end;q++)
    {
		if((q&1)==0)
			continue;      // Jacobi not defined for even q
		res=0.0;
		q_minus_s=pow((double)q,-s);
		for(num=1;num<q;num++)
		  if(co_prime(num,q))
			  {
				frac=double(num)/q;
				dsteps=double(num)/(gap*q);
				csteps=ceil(dsteps);
				steps=(int) csteps;
				if((csteps-dsteps)>=0.5) // use lower row of R_ns
				{
					gap_ptr=r_n_terms*(no_gaps-steps);
					x=double(steps)-frac/gap;
					out_val=calc_r_n(r_n_terms,s_array,x,&r_n_vals[gap_ptr],frac);
				}
				else // use upper row
				{
					gap_ptr=r_n_terms*(no_gaps-steps+1);
					x=double(steps)-(frac/gap+1.0);
					out_val=calc_r_n_neg(r_n_terms,s_array,x,&r_n_vals[gap_ptr],frac);
				}
				// out_val is now ~ sum_{n=RN}^\infty (n+num/q)^(-s)
				out_val=out_val+pow(frac,minus_s);
				if(kronecker(num,q)==1)
					res+=out_val;
				else
					res-=out_val;
		  }
		  res*=q_minus_s;
		  printf("res=%10.8e\n",res);
  }
  return(1);
}

int main(int argc, char **argv)
{

  double s,*s_array,*r_n_vals;
  int no_gaps,num_s,num_fracs,N,rn,file_N;
  long long int q_start,q_end;
  double gap;
  FILE *in_file,*out_file;
  clock_t no_clicks;

  no_clicks=clock(); // start timing

  if(argc!=6)
    print_usage();
  q_start=_atoi64(argv[1]);
  if(q_start<3)
    print_usage();
  q_end=_atoi64(argv[2]);
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
  if(!(s_array=(double *) _aligned_malloc(sizeof(double)*N,16)))
    fatal_error("Can't allocate memory for s_array. Exiting.");
  
  
  cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=(double *) _aligned_malloc(sizeof(double)*(no_gaps+1)*N,16)))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");


  fwrite(&q_start,sizeof(int),1,out_file);
  fwrite(&q_end,sizeof(int),1,out_file);
  fwrite(&num_s,sizeof(int),1,out_file);
  printf("Going to try to process %d values of s.\n",num_s);

  s=0.5;

		skip_double(in_file); // im_s
		skip_dcomplex(in_file,2); // gammas

		num_fracs=calc_hurwitz1(s,N,file_N,gap,s_array,r_n_vals,no_gaps,q_start,q_end,in_file,out_file,rn);
		if(!num_fracs)
			fatal_error("Error running Hurwitz routine. Exiting.");

	fclose(in_file);
  fclose(out_file);

  printf("Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  printf("Fractions processed=%d.\n",num_fracs);
  return(0);
}


