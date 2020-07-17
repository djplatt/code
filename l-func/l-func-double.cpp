/*

File: l-func-double.cpp

Created: 6th July 2008

Version: <v> = 1.0

Last Modified: 6th July 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: 

Build instructions: g++ -ol-func-double l-func-double.cpp -O3

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
#include <complex>
#include <time.h>

#define Q (10000) /* maximum conductor */

#define T_MAX (4000) /* maximum height */
#define T_POINTS (1000) /* number of points up T line to use */


#define SUCCESS (1)
#define FAILURE (0)
#define OUT_VEC_SIZE (256)
#define QS_PER_FILE (100)


using namespace std;

typedef complex<double> dcomplex;

typedef struct {double val; 
  unsigned int num; 
  unsigned int den;
  dcomplex r_n;} frac;

union out_type{double d;char c;} out_vec[OUT_VEC_SIZE];

void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: l-func%s (no. terms) <n> (fname) (gap)\n",VERSION);
	printf("  (no. terms) - integer > 0\n");
	printf("  (N)       - integer > 0\n");
	printf("  (fname)   - file with s, zeta values\n");
	printf("  (gap)     - 0.0 < gap < 1.0\n");
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

ofstream out_file;

dcomplex buffer[OUT_VEC_SIZE];
unsigned int buffer_ptr=0;

void flush_buffer()
{
  out_file.write((char *) &buffer, sizeof(dcomplex)*buffer_ptr);
  buffer_ptr=0;
};


void write_dcomplex(dcomplex z)
{
  buffer[buffer_ptr++]=z;
  if(buffer_ptr==OUT_VEC_SIZE)
    flush_buffer();
};

void close_file()
{
  if(buffer_ptr!=0)
    flush_buffer();
  out_file.close();
};

void open_file(unsigned int f_num,dcomplex s)
{
  char str[100];

  if(sprintf(str,"data//hurfile%04d%6.4f.dat",f_num,imag(s))<0)
    fatal_error("Unable to create filename. Exiting.");
  out_file.open(str,ios::out|ios::binary|ios::trunc);
  if(!out_file.is_open())
    fatal_error("Couldn't open file. Exiting.");
};

dcomplex pow(double x,dcomplex s)
  /* return x^s */
{
  double tlnx;
  double x_to_sigma;

  tlnx=imag(s)*log(x);
  x_to_sigma=pow(x,real(s));
  return(dcomplex(x_to_sigma*cos(tlnx),x_to_sigma*sin(tlnx)));
};

dcomplex s_n(dcomplex s, double alpha, int N)
  /* return sum from 0 to N of (n+alpha)^(-s) 
     N will probably end up being 0 or 1 so might
     revisit */
{
  int n;
  dcomplex res(0,0);

  for(n=0;n<N;n++)
    res+=pow(n+alpha,-s);
  //  cout << "s_n called with " << s << " " << alpha << " " << N << " returning " << res << endl;
  return(res);
};

/*
void create_s_array(dcomplex s,dcomplex *s_array,int n)
   fill the n*n array s_array with the taylor coefficients 
   s_array[i,j] <- (s+i)(s+i+1)...(s+j)/(j-i+1)! 
{
  int i,j;
  double fact;


  for(i=0;i<n;i++)
    s_array[i*n+i]=dcomplex(real(s)+i,imag(s));

  for(i=0;i<n-1;i++)
    for(j=i+1;j<n;j++)
      s_array[i*n+j]=s_array[i*n+j-1]*s_array[j*n+j];
  fact=2;
  for(i=0;i<n-1;i++)
    {
    for(j=0;j<n-i-1;j++)
      s_array[j*n+j+1+i]/=fact;
    fact=fact*(i+3);
    };
};
*/

void create_s_array(dcomplex s,dcomplex *s_array, unsigned int n)
{
  unsigned int i,del;
  double fact;

  for(i=0;i<n;i++)
    s_array[i*n+i]=dcomplex(real(s)+i,imag(s));

  del=1;
  fact=2;
  while(del<n)
    {
      for(i=0;i<n-del;i++)
	s_array[i*n+i+del]=s_array[i*n+i]*s_array[i*n+n+i+del]/fact;
      del++;
      fact++;
    };
  /*
  for(i=0;i<n;i++)
    for(del=i;del<n;del++)
      {
	cout << "Re s_array["<<i<<","<<del<<"] = ";
	printf("%20.18f\n",real(s_array[i*n+del]));
      };
  exit(0);
  */
};


void print_r_n_vec(dcomplex *r_n_vec,unsigned int n)
{
  unsigned int i;

  for(i=0;i<n;i++)
    cout << "r_n[" << i << "] is " << r_n_vec[i] << endl;
};



void make_r_ns (dcomplex *s_array,
		unsigned int n, 
		double delta,
		dcomplex *this_r_n_vec,
		dcomplex *last_r_n_vec,
		double *delta_n)
{
  unsigned int i,j;

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
  //print_r_n_vec(this_r_n_vec,n);
};



dcomplex calc_r_n(unsigned int r_n_terms, dcomplex *s_array, double delta,
		  dcomplex *last_r_n_vec)
{
  unsigned int i;
  dcomplex res,term;
  double delta_n;

  delta_n=delta;
  res=last_r_n_vec[0];

  /* could decide how many terms to take according to delta */

  for(i=0;i<r_n_terms-1;)
    {
      term=s_array[i]*delta_n;
      i++;
      term*=last_r_n_vec[i];
      res+=term;
      delta_n*=delta;
    };
  return(res);
};

dcomplex make_out_val(dcomplex r_n,dcomplex *n_s,
		      unsigned int a, unsigned int q, unsigned int N)
{
  unsigned int n;
  dcomplex res;

  for(n=0;n<N;n++)
    res+=n_s[n*q+a];

  return(r_n*n_s[q]+res);
};

int calc_hurwitz1(dcomplex s, unsigned int r_n_terms, unsigned int N,
		  dcomplex *zeta_vals,double gap)
{
  dcomplex s0,*s_array;
  unsigned int file_no,i,no_gaps,gap_ptr,q,num,num_fracs;
  dcomplex *r_n_vals,out_val,*n_s;
  double x,frac,dq,*delta_ns;

  no_gaps=(int) ceil(1.0/gap);

  cout << "Running with " << no_gaps << " gaps." << endl;
  cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=new dcomplex[no_gaps*r_n_terms]))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");

  cout << "Allocating memory for s_array." << endl;
  if(!(s_array=new dcomplex[(r_n_terms-1)*(r_n_terms-1)]))
    fatal_error("Can't allocate memory for s_array. Exiting.");

  create_s_array(s,s_array,r_n_terms-1);

  if(!(n_s=new dcomplex[N*Q+1]))
    fatal_error("Can't allocate memory for n_s. Exiting.");
if(!(delta_ns=new double[r_n_terms-1]))
	fatal_error("Can't allocate memory for delta_ns. Exiting.");

  n_s[1]=1;
  for(i=2;i<=N*Q;i++)
    n_s[i]=pow(i,-s);

  s0=s;

  for(i=0;i<r_n_terms;i++)
    {
      r_n_vals[i]=zeta_vals[i]-s_n(s0,1.0,N);
      //cout << s0 << " " << zeta_vals[i] << " " << s_n(s0,1.0,N) << endl;
      //cout << "r_n[" << i << "] is " << last_big_r_n_vec[i] << endl;
      s0=dcomplex(real(s0)+1,imag(s0));
    };

  for(i=1;i<no_gaps;i++)
    make_r_ns(s_array,r_n_terms,gap,&r_n_vals[i*r_n_terms],
	      &r_n_vals[(i-1)*r_n_terms],delta_ns);
 
  //print_r_n_vec(r_n_vals,r_n_terms*no_gaps);
  open_file(0,s);
  file_no=1;
  /*
  out_val=make_out_val(r_n_vals[0],n_s,1,1,N);
  write_dcomplex(out_val);
  cout << "1 1" << endl;
  */
  for(q=3;q<=Q;q++)
    {
      if((q&3)==2)
	continue;      // if q = 2 mod 4 then no primitive characters

      num_fracs=0;
      gap_ptr=no_gaps*r_n_terms;
      x=1.0-gap*no_gaps;
      dq=(double) q;
      if((q%QS_PER_FILE)==1)
	{
	  flush_buffer();
	  out_file.close();
	  open_file(file_no++,s);
	};
      for(num=1;num<q;)
	{
	  if(co_prime(num,q))
	    {
	      num_fracs++;
	      frac=(double) num/dq;
	      while(frac>=x)
		{
		  x=x+gap;
		  gap_ptr-=r_n_terms;
		};
	      out_val=calc_r_n(r_n_terms,s_array,
			       x-frac,&r_n_vals[gap_ptr]);
	      //cout << "out_val1= " << out_val << endl;
	      out_val=make_out_val(out_val,n_s,num,q,N);
	      //cout << "out_val2= " << out_val << endl;
	      //cout << out_val << endl;
	      write_dcomplex(out_val);
	    };
	  num++;
	  if(num==q)
	    {
	      //cout << num << " " << num_fracs << endl;
	      break;
	    };
	};
    };
  close_file();
	  

  return(SUCCESS);
};

 
int main(int argc, char **argv)
{
  double zeta_valr,zeta_vali;
  dcomplex s,*zeta_vals;
  unsigned int r_n_terms,N,i;
//  char *fname;
  double gap;
  ifstream zeta_file;

  clock_t no_clicks;

  no_clicks=clock();


  if(argc!=5)
    print_usage();
  r_n_terms=atoi(argv[1]);
  if(r_n_terms<=0)
    print_usage();
  N=atoi(argv[2]);
  if(N<=0)
    print_usage();
  zeta_file.open(argv[3]);
  if(!zeta_file.is_open())
    fatal_error("Couldn't open zeta file. Exiting.\n");
  gap=atof(argv[4]);
  if((gap<=0.0)||(gap>=1.0))
    print_usage();
  if(!(zeta_vals=new dcomplex[r_n_terms]))
    fatal_error("Couldn't allocate memory for zeta_vals. Exiting.\n");
  zeta_file >> zeta_valr;
      if(zeta_file.eof())
	fatal_error("Insufficient zeta values given. Exiting.\n");
  zeta_file >> zeta_vali;
  s=dcomplex(zeta_valr,zeta_vali);
  for(i=0;i<r_n_terms;i++)
    {
      if(zeta_file.eof())
	fatal_error("Insufficient zeta values given. Exiting.\n");
      zeta_file >> zeta_valr;
      if(zeta_file.eof())
	fatal_error("Insufficient zeta values given. Exiting.\n");
      zeta_file >> zeta_vali;
      zeta_vals[i]=dcomplex(zeta_valr,zeta_vali);
 
    };
  zeta_file.close();


  printf("num terms = %20d\nN       = %20d\n",
	 r_n_terms,N);




  cout << "Entering calc_hurwitz1." << endl;
  if(!calc_hurwitz1(s,r_n_terms,N,zeta_vals,gap))
	  fatal_error("Error running Hurwitz routine. Exiting.");

   printf("Time Elapsed = %8.5f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));

  return(SUCCESS);
};

