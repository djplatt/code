/*

File: l-func2.0.cpp

Created: 2nd July 2008

Version: <v> = 2.0

Last Modified: 2nd Juky 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: Might get away with 1/2 of the fracs array
                      if its symmetric about 1/2, but that part is
		      pretty quick anyway.

Build instructions: g++ -ol-func<v> l-func<v>.cpp -O2

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "2.0"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>

#define Q (10000) /* maximum conductor */
//#define GAP ((double) 0.01)

#define T_MAX (100) /* maximum height */
#define T_POINTS (500) /* number of points up T line to use */


#define SUCCESS (1)
#define FAILURE (0)
#define OUT_VEC_SIZE (64)

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
	printf("  (N)       - integer >= 0\n");
	printf("  (fname)   - file with s, zeta values\n");
	printf("  (gap)     - 0.0 < gap < 1.0\n");
	exit(FAILURE);
};

void fatal_error(char *error_string)
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

inline int co_prime(unsigned int a, unsigned short int b)
{
  return(gcd(a,b)==1);
};

/* vector to hold the Farey Fractions */
frac fracs[Q<<1];

inline int _frac_comp (frac a, frac b)
/* comparison function for qsort 
   note no two fractions can't be equal so
   return(0) isn't a possibility. See
   qsort documentation (stdlib.h) */

{
  if(b.val-a.val>0)
    return(1);
  return(-1);
};

int frac_comp(const void *a,const void *b)
{
  return(_frac_comp(*(frac*)a,*(frac*)b));
};  

unsigned int a,b,c,d;

void init_more_fracs()
{
  a=1;
  b=1;
  c=Q-1;
  d=Q;
};

ofstream *files;

bool open_files(unsigned int num_files)
{
  unsigned int i;
  char str[20];

  for(i=0;i<num_files;i++)
    {
      if(sprintf(str,"hurfile%04d.dat",i)<0)
	fatal_error("Unable to create filenames. Exiting.");
      files[i].open(str,ios::out|ios::binary|ios::trunc);
      files[i].write(&out_vec,sizeof(double)*OUT_VEC_SIZE);
      // if(!files[i].is_open()); // this no work!
      //fatal_error("Couldn't open file. Exiting.");
    };
  return(true);
};

bool close_files(unsigned int num_files)
{
  unsigned int i;
  for(i=0;i<num_files;i++)
    files[i].close();
  return(true);
};


int more_fracs()
  /* put the next Q fractions into fracs. Return Q or the number of fracs
     actualy constructed. */
{
  unsigned int e,f,k,ptr;

  ptr=1;
  fracs[0].num=a;
  fracs[0].den=b;
  fracs[0].val=(double) a/(double) b;
  while((a>0)&&(ptr<Q))
    {
      k=(Q+b)/d;
      e=k*c-a;
      f=k*d-b;
      a=c;
      b=d;
      c=e;
      d=f;
      fracs[ptr].num=a;
      fracs[ptr].den=b;
      fracs[ptr++].val=(double) a/(double) b;
    };

  //  cout << "Exiting more_fracs." << endl;

  return(ptr);
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


void create_s_array(dcomplex s,dcomplex *s_array,int n)
  /* fill the n*n array s_array with the taylor coefficients */
  /* s_array[i,j] <- (s+i)(s+i+1)...(s+j)/(j-i+1)! */
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

dcomplex *this_r_n_vec,*last_r_n_vec;


void make_r_ns (dcomplex *s_array,
		unsigned int n, 
		double delta)
{
  int i,j;
  double delta_n[n-1];

  /*
  cout << "Entered make r_n with n=" << n << " and delta =" << delta << endl;
  for(i=0;i<n;i++)
    cout << old_r_n_vec[i] << endl;
  cout << "***" << endl;
  */

  //cout << "In make_r_ns." << endl;
  delta_n[0]=delta;
  for(i=1;i<n-1;i++)
    delta_n[i]=delta_n[i-1]*delta;
  for(i=0;i<n;i++)
    {
      this_r_n_vec[i]=last_r_n_vec[i];
      for(j=0;j<n-i-1;j++)
	this_r_n_vec[i]+=s_array[i*(n-1)+j+i]*delta_n[j]*last_r_n_vec[j+i+1];
    };
  /*
  for(i=0;i<n;i++)
    cout << new_r_n_vec[i] << endl;
  cout << "Exiting make_r_ns." << endl;
  */
};


void print_r_n_vec(dcomplex *r_n_vec,unsigned int n)
{
  unsigned int i;

  for(i=0;i<n;i++)
    cout << "r_n[" << i << "] is " << r_n_vec[i] << endl;
};


dcomplex calc_r_n(unsigned int r_n_terms, dcomplex *s_array, double delta)
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

double GAP;

int calc_hurwitz(unsigned int r_n_terms,
		 dcomplex *s_array,
		 double *last_p_over_Q)
{
  unsigned int num_fracs,farey_ptr,i;
  double delta;

  num_fracs=more_fracs();
  /*
  */

  if(fracs[num_fracs-1].val==0.0)
    num_fracs--;

  fracs[0].r_n=fracs[Q-1].r_n;

  for(farey_ptr=1;farey_ptr<num_fracs;farey_ptr++)
    {
      delta=last_p_over_Q[0]-fracs[farey_ptr].val;
      if(delta<GAP)
	fracs[farey_ptr].r_n=calc_r_n(r_n_terms,s_array,delta);
      else
      {
	/*
	cout << "Making new r_n's " << fracs[farey_ptr].num << endl;
	*/
	make_r_ns(s_array,r_n_terms,delta);
	
	//printf("this %d last %d\n",this_r_n_vec,last_r_n_vec);
	swap(this_r_n_vec,last_r_n_vec);

	fracs[farey_ptr].r_n=last_r_n_vec[0];
	last_p_over_Q[0]=fracs[farey_ptr].val;
      };
    };
  /*
  cout << fracs[0].num << "/" << fracs[0].den << " ";
  printf("%20.18f\n",real(fracs[0].r_n));
  cout << fracs[1].num << "/" << fracs[1].den << " ";
  printf("%20.18f\n",real(fracs[1].r_n));
  cout << fracs[num_fracs-2].num << "/" << fracs[num_fracs-2].den << " ";
  printf("%20.18f\n",real(fracs[num_fracs-2].r_n));
  */
  //cout << fracs[num_fracs-1].num << "/" << fracs[num_fracs-1].den << " ";
  //printf("%20.18f\n",real(fracs[num_fracs-1].r_n));
  

  return(num_fracs);

};


int calc_hurwitz1(dcomplex s, unsigned int r_n_terms, unsigned int N,
		  dcomplex *zeta_vals)
{
  dcomplex s0,*s_array;
  double *last_p_over_Q;
  unsigned int seg_no,i,j,num_fracs,num_files;

  num_files=(int) ceil((double Q)/1000.0);

  if(num_files>9999)
    fatal_error("Can't open more than 9999 file. Exiting.");

  if(!(files=new ofstream[num_files]))
    fatal_error("Can't allocate memory for files. Exiting.");

  if(!open_files(num_files))
    fatal_error("Failed to open files. Exiting.");

  if(!(last_p_over_Q=new double))
    fatal_error("Cant allocate memory for last_p_over_Q. Exiting.");

  last_p_over_Q[0]=1.0;

  if(!(s_array=new dcomplex[(r_n_terms-1)*(r_n_terms-1)]))
    fatal_error("Can't allocate memory for s_array. Exiting.");
  create_s_array(s,s_array,r_n_terms-1);


  /*
  for(j=0;j<r_n_terms-1;j++)
    for(i=0;i<=j;i++)
      cout << "s_array[" << i << "," << 
	j << "] is " << s_array[i*(r_n_terms-1)+j] << endl;
  */

  if(!(last_r_n_vec=new dcomplex[r_n_terms]))
    fatal_error("Couldn't allocate memory for last_r_n_vec. Exiting");
  if(!(this_r_n_vec=new dcomplex[r_n_terms]))
    fatal_error("Couldn't allocate memory for this_r_n_vec. Exiting");

  s0=s;
  for(i=0;i<r_n_terms;i++)
    {
      last_r_n_vec[i]=zeta_vals[i]-s_n(s0,1.0,N);
      //cout << s0 << " " << zeta_vals[i] << " " << s_n(s0,1.0,N) << endl;
      //cout << "r_n[" << i << "] is " << last_big_r_n_vec[i] << endl;
      real(s0)++;
    };
  fracs[Q-1].r_n=last_r_n_vec[0];

  init_more_fracs();
  while(a>0)
    num_fracs=calc_hurwitz(r_n_terms,s_array,last_p_over_Q);

  printf("Value of %d/%d was %20.18f\n",fracs[num_fracs-1].num,
	 fracs[num_fracs-1].den,
	 real(fracs[num_fracs-1].r_n));
  num_fracs--;
  printf("Value of %d/%d was %20.18f\n",fracs[num_fracs-1].num,
	 fracs[num_fracs-1].den,
	 real(fracs[num_fracs-1].r_n));
  return(SUCCESS);
};

 
int main(int argc, char **argv)
{
  double zeta_val;
  dcomplex s,*zeta_vals;
  unsigned int r_n_terms,N,i;
  char *fname;
  ifstream zeta_file;

  if(argc!=5)
    print_usage();
  r_n_terms=atoi(argv[1]);
  if(r_n_terms<=0)
    print_usage();
  N=atoi(argv[2]);
  if(N<0)
    print_usage();
  zeta_file.open(argv[3]);
  if(!zeta_file.is_open())
    fatal_error("Couldn't open zeta file. Exiting.\n");
  GAP=atof(argv[4]);
  if((GAP<=0.0)||(GAP>=1.0))
    print_usage();
  if(!(zeta_vals=new dcomplex[r_n_terms]))
    fatal_error("Couldn't allocate memory for zeta_vals. Exiting.\n");
  zeta_file >> zeta_val;
  real(s)=zeta_val;
      if(zeta_file.eof())
	fatal_error("Insufficient zeta values given. Exiting.\n");
  zeta_file >> zeta_val;
  imag(s)=zeta_val;
  for(i=0;i<r_n_terms;i++)
    {
      if(zeta_file.eof())
	fatal_error("Insufficient zeta values given. Exiting.\n");
      zeta_file >> zeta_val;
      real(zeta_vals[i])=zeta_val;
      if(zeta_file.eof())
	fatal_error("Insufficient zeta values given. Exiting.\n");
      zeta_file >> zeta_val;
      imag(zeta_vals[i])=zeta_val;
 
    };
  zeta_file.close();


  printf("num terms = %20d\nN       = %20d\n",
	 r_n_terms,N);




  cout << "Entering calc_hurwitz1." << endl;
  if(!calc_hurwitz1(s,r_n_terms,N,zeta_vals))
    fatal_error("Error running Hurwitz routine. Exiting.");


  return(SUCCESS);
};

