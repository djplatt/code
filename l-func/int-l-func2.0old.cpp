/*

File: int-l-func2.0.cpp

Created: 23rd August 2008

Version: <v> = 2.0

Last Modified: 1st September 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: 

Build instructions: g++ -oint-l-func<v> int-l-func<v>.cpp
                    ../sse/IntervalsSSE2.cpp -I../sse -O3

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "1.0"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "IntervalComplex.h"  // the complex interval class

#define T_MAX (4000) /* maximum height */
#define T_POINTS (10000) /* number of points up T line to use */


#define SUCCESS (1)
#define FAILURE (0)
#define OUT_VEC_SIZE (256)
#define QS_PER_FILE (100)

#define debug printf("got to line %d\n",__LINE__);

using namespace std;

typedef int_complex dcomplex;
typedef MachineEstimate int_double;

void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: int-l-func%s (q-start) (q-end) (no. terms) <n> (ifname) (ofname) (gap)\n",VERSION);
	printf("  (q-start)   - integer >=3\n");
	printf("  (q-end)     - integer >=qstart\n");
	printf("  (no. terms) - integer > 0\n");
	printf("  (N)         - integer > 0\n");
	printf("  (ifname)    - file with s, zeta values\n");
	printf("  (ofname)    - output file name\n");
	printf("  (gap)       - 0.0 < gap < 1.0\n");
	exit(FAILURE);
};

void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(FAILURE);
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


dcomplex s_n(const dcomplex &s, const int_double &alpha, int N)
  /* return sum from 0 to N-1 of (n+alpha)^(-s) 
     N will probably end up being 0 or 1 so might
     revisit */
{
  int n;
  dcomplex res(0.0,0.0);
  dcomplex s0(-s.real,-s.imag);	
  for(n=0;n<N;n++)
    res+=pow(alpha+n,s0);
  //  cout << "s_n called with " << s << " " << alpha << " " << N << " returning " << res << endl;
  return(res);
};

void create_s_array(const dcomplex &s,dcomplex *s_array, unsigned int n)
{
  unsigned int i,del;
  double fact;

  //printf("in ceate_s_array\n");
  s_array[0]=s;
  //print_int_complex(s_array[0]);
  for(i=1;i<n;i++)
      s_array[i*n+i]=dcomplex(real(s)+i,imag(s));
  //debug;
  del=1;
  fact=2;
  while(del<n)
    {
      for(i=0;i<n-del;i++)
	{
	  //debug;
	  s_array[i*n+i+del]=s_array[i*n+i]*s_array[i*n+n+i+del]/fact;
	  //print_int_complex(s_array[i*n+i+del]);printf("\n");
	};
      del++;
      fact++;
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
};

void make_r_ns (dcomplex *s_array,
		unsigned int n, 
		double delta,
		dcomplex *this_r_n_vec,
		dcomplex *last_r_n_vec)
{
  unsigned int i,j;
  double *delta_n;
  
  delta_n=new double[n-1];

  /*
  cout << "Entered make r_n with n=" << n << " and delta =" << delta << endl;
  for(i=0;i<n;i++)
    cout << old_r_n_vec[i] << endl;
  cout << "***" << endl;
  */

  //cout << "In make_r_ns." << endl;
  
  //print_r_n_vec(last_r_n_vec,n);
  delta_n[0]=delta;
  for(i=1;i<n-1;i++)
    delta_n[i]=delta_n[i-1]*delta;
  for(i=0;i<n;i++)
    {
      this_r_n_vec[i]=last_r_n_vec[i];
      for(j=0;j<n-i-1;j++)
	{
	  //cout << "delta_n[" << j << "]= " << delta_n[j] << endl;
	  this_r_n_vec[i]+=s_array[i*(n-1)+j+i]*delta_n[j]*last_r_n_vec[j+i+1];
	  //cout << "this_r_n_vec[" << i << "] = " << this_r_n_vec[i] << endl;
	};
    };
  // give delta_n back
  //print_r_n_vec(this_r_n_vec,n);
};

dcomplex calc_r_n(unsigned int r_n_terms, dcomplex *s_array, const int_double &delta,
		  dcomplex *r_n_vec)
{
  unsigned int i;
  dcomplex res,term;
  int_double delta_n;

// printf("in calc_r_n:\n");
//  print_r_n_vec(r_n_vec,r_n_terms);
//  printf("delta = ");print_interval(delta);printf("\n");
  delta_n=delta;
  res=r_n_vec[0];

  /* could decide how many terms to take according to delta */

  for(i=0;i<r_n_terms-1;)
    {
      term=s_array[i]*delta_n;
      //printf("term = ");print_int_complex(term);printf("\n");
      i++;
      //printf("last_r_n_vec[%d] = ",i);print_int_complex(last_r_n_vec[i]);printf("\n");
      term*=r_n_vec[i];
      //printf("term = ");print_int_complex(term);printf("\n");
      res+=term;
      //printf("res = ");print_int_complex(res);printf("\n");
      delta_n*=delta;
      //printf("delta_n = ");print_int_double(delta_n);printf("\n");
    };
  return(res);
};


dcomplex make_out_val(const dcomplex &r_n,dcomplex *n_s,
		      unsigned int a, unsigned int q, unsigned int N)
{
  unsigned int n;
  dcomplex res;

  for(n=0;n<N;n++)
  {
      res+=n_s[n*q+a];
//      printf("n_s: ");print_int_complex(n_s[n*q+a]);printf("\n");
//      printf("res: ");print_int_complex(res);printf("\n");
  };

//  printf("a: %d, q: %d\nr_n: ",a,q);
//  print_int_complex(r_n);
//  printf("\nres: ");
//  print_int_complex(res);
//  printf("\nq^(-s): ");
//  print_int_complex(n_s[q]);
//  printf("\nr_n*q^(-s): ");
//  print_int_complex(r_n*n_s[q]);
//  printf("\nr_n*q^(-s)+res: ");
//  print_int_complex(r_n*n_s[q]+res);
//  printf("\n\n");

  return(n_s[q]*r_n+res);
};

int calc_hurwitz1(const dcomplex &s, unsigned int r_n_terms, unsigned int N,
		  dcomplex *zeta_vals,double gap,dcomplex *s_array,
		  dcomplex *r_n_vals,unsigned int no_gaps,
		  dcomplex *n_s, unsigned int q_start,
		  unsigned int q_end, ofstream *out_file)
{
  dcomplex out_val;
  unsigned int i,gap_ptr,q,num,steps;
  int_double x,frac;
  dcomplex s0=dcomplex(-s.real,-s.imag);

  //printf("in calc hurwitz with s = ");print_int_complex(s);printf("\n");

  create_s_array(s,s_array,r_n_terms-1);
  //out_file->write((char *) s_array,sizeof(dcomplex)*(r_n_terms-1)*(r_n_terms-1));
  /*
  for(i=0;i<(r_n_terms-1)*(r_n_terms-1);i++)
  {
      printf("s_array[%d] =",i);
      print_int_complex(s_array[i]);
      printf("\n");
  };
  */
  n_s[1]=int_complex(int_double(1.0,1.0),int_double(0.0,0.0));
  for(i=2;i<=N*q_end;i++)
  {
    n_s[i]=pow(i,s0);
    //printf("n_s[%d]= ",i);print_int_complex(n_s[i]);printf("\n");
  };
  //out_file->write((char *) n_s,sizeof(dcomplex)*(N*q_end+1));

  s0=s;

  for(i=0;i<r_n_terms;i++)
    {
      r_n_vals[i]=zeta_vals[i]-s_n(s0,1.0,N);
      s0+=1.0;
    };

  for(i=1;i<no_gaps;i++)
    make_r_ns(s_array,r_n_terms,gap,&r_n_vals[i*r_n_terms],
	      &r_n_vals[(i-1)*r_n_terms]);

  //out_file->write((char *) r_n_vals,sizeof(dcomplex)*r_n_terms*no_gaps); 

  //  print_r_n_vec(r_n_vals,r_n_terms*no_gaps);

  //  open_file(0,s);
  //file_no=1;

  for(q=q_start;q<=q_end;q++)
    {
		if((q&3)==2)
			continue;      // if q = 2 mod 4 then no primitive characters

      //gap_ptr=no_gaps*r_n_terms;
      //      x=1.0-gap*no_gaps;
      //if((q%QS_PER_FILE)==1)
      //	{
      //  flush_buffer();
      //  out_file.close();
      //  open_file(file_no++,s);
      //};
      for(num=1;num<q;num++)
		  if(co_prime(num,q))
			  {
				frac=num;
	    //	      printf("frac = ");print_int_double(frac);printf("\n");
				frac=frac/q;
	    //	      printf("frac = ");print_int_double(frac);printf("\n");
				steps=(int) ceil(num/(gap*q));
				gap_ptr=r_n_terms*(no_gaps-steps);
				x=steps*gap-frac;
	    
	    /*
	      while(frac>x)
	      {
	      x=x+gap;
	      gap_ptr-=r_n_terms;
	      };
	    */
	    //print_int_double(x-frac);printf("\n%d\n",gap_ptr);
				out_val=calc_r_n(r_n_terms,s_array,
					x,&r_n_vals[gap_ptr]);
	    //	      printf("out_val = ");print_int_complex(out_val);printf("\n");
				out_val=make_out_val(out_val,n_s,num,q,N);
				//printf("out_val = ");print_int_complex(out_val);printf("\n");
				//printf("fwrite returned %d\n",fwrite(&out_val,sizeof(dcomplex),1,out_file));
				out_file->write((char *) &out_val,sizeof(dcomplex));
	  };
    };
  return(SUCCESS);
};

int main(int argc, char **argv)
{
  MachineEstimate::Computation comp;

  double ds;

  // hardwired for now to circumvent new dcomplex[n] problem.
  // make sure runs with 3 <=1000 30 1 <f1> <f2> 0.125 
  dcomplex s,*zeta_vals,*s_array,*r_n_vals,*n_s;



  int r_n_terms,N,i,j,no_gaps,q_start,q_end,num_s,num_fracs,num_zeta;
  int_double s_delta;
  double gap;
  ifstream zeta_file;
  ofstream out_file;


  if(argc!=8)
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
  if(N<=0)
    print_usage();
  zeta_file.open(argv[5]);
  if(!zeta_file.is_open())
    fatal_error("Couldn't open ifile. Exiting.\n");
  out_file.open(argv[6],ios :: out | ios :: binary);
  if(!out_file.is_open())
    fatal_error("Couldn't open ifile. Exiting.\n");
  gap=atof(argv[7]);
  if((gap<=0.0)||(gap>=1.0))
    print_usage();
  
  printf("dcomplex: %d int_double: %d\n",sizeof(dcomplex),sizeof(int_double));
  if(!(zeta_vals= (dcomplex *) _aligned_malloc(sizeof(dcomplex)*r_n_terms,16)))
    fatal_error("Couldn't allocate memory for zeta_vals. Exiting.\n");

  cout << "Allocating memory for s_array." << endl;
  if(!(s_array=(dcomplex *) _aligned_malloc(sizeof(dcomplex)*(r_n_terms-1)*(r_n_terms-1),16)))
    fatal_error("Can't allocate memory for s_array. Exiting.");
  
  no_gaps=(int) ceil(1.0/gap);
  cout << "Running with " << no_gaps << " gaps." << endl;
  
  cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=(dcomplex *) _aligned_malloc(sizeof(dcomplex)*no_gaps*r_n_terms,16)))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");

  cout << "Allocating memory for n_s." << endl;
  if(!(n_s=(dcomplex *) _aligned_malloc(sizeof(dcomplex)*(N*q_end+1),16)))
    fatal_error("Can't allocate memory for n_s. Exiting.");
  

  out_file.write((char *) &q_start,sizeof(int));
  out_file.write((char *) &q_end,sizeof(int));

  num_fracs=0;
  for(i=q_start;i<=q_end;i++)
    if((i&3)!=2)
      for(j=1;j<i;j++)
		if(co_prime(i,j))
		  num_fracs++;

  printf("Going to generate %d values per s.\n",num_fracs);
  out_file.write((char *) &num_fracs,sizeof(int));

  zeta_file >> num_zeta; // number of zeta terms per s
  if(num_zeta<r_n_terms)
    fatal_error("Not enough zeta values in zeta_file. Exiting.\n");
  zeta_file >> num_s;    // number of s's
  zeta_file >> ds;      // step between each s.
  s_delta=int_double(ds);    // assumes delta representable

  printf("Going to try and process %d values of s using step size\n    ",
	 num_s);
  print_interval(s_delta);printf("\n");
  
  s.real=0.5;
  for(i=0;i<num_s;i++)
    {
      for(j=0;j<r_n_terms;j++)
	{
	  zeta_file >> ds;
	  zeta_vals[j].real=int_double(ds,nextafter(ds,DBL_MAX));
	  zeta_file >> ds;
	  zeta_vals[j].imag=int_double(ds,nextafter(ds,DBL_MAX));
	  //print_int_complex(zeta_vals[j]);printf("\n");
	};

      for(;j<num_zeta;j++)
	  zeta_file >> ds >> ds;    // throw away extra ones

      //out_file.write((char *) zeta_vals,sizeof(dcomplex)*r_n_terms);

      if(!calc_hurwitz1(s,r_n_terms,N,zeta_vals,gap,s_array,
			r_n_vals,no_gaps,n_s,q_start,q_end,&out_file))
	fatal_error("Error running Hurwitz routine. Exiting.");

      s.imag+=s_delta;
    };

  zeta_file.close();
  out_file.close();

  cout << "Enter a number to quit. ";
  cin >> i;
  return(SUCCESS);
};


