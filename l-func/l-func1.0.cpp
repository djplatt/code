/*

File: l-func1.0.cpp

Created: 27th June 2008

Version: <v> = 1.0

Last Modified: 1st Juky 2008

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

#define VERSION "1.0"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>

#define Q (10000) /* maximum conductor use 1000 for debug */
#define NUM_FRACS (30397486) /* number of distint Farey fractions
				 0 < p/q <= 1 with q <= Q  use 304192         */
#define T_MAX 100 /* maximum height */
#define T_POINTS 500 /* number of points up T line to use */


#define SUCCESS (1)
#define FAILURE (0)

using namespace std;

typedef complex<double> dcomplex;

typedef struct {double val; 
  unsigned short int num; 
  unsigned short int den;
  dcomplex z;} frac;

void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: l-func%s (terms 1) (gap) (terms 2) <n> <fname>\n",VERSION);
	printf("  (terms 1) - integer > (terms 2)\n");
	printf("  (gap)     - double  > 0.0 < 1.0\n");
	printf("  (terms 2) - integer > 0\n");
	printf("  (N)       - integer >= 0\n");
	printf("  (fname)   - file with s, zeta values\n");
	exit(FAILURE);
};

void fatal_error(char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(FAILURE);
};

inline unsigned short int gcd (unsigned short int a, unsigned short int b)
  /* Euclid algorithm gcd */
{
  unsigned short int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

inline int co_prime(unsigned short int a, unsigned short int b)
{
  return(gcd(a,b)==1);
};

/* vector to holf the Farey Fractions */
frac fracs[NUM_FRACS];

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

int make_fracs()
  /* build a vector of Farey fractions p/q such that
     (p,q)=1, p<=q, p>0, q<=Q
     sort them so that fracs[0].val is 1/1 */
{
  unsigned short int num,den;
  unsigned int ptr,ptr2;

  cout << "Creating Farey fractions." << endl;
  
  fracs[0].val=1;
  fracs[0].den=1;
  fracs[0].num=1;
  ptr=1;
  for(den=1;den<=Q;den++)
    for(num=1;num<den;num++)
      if(co_prime(den,num))
	{
	  if(ptr==NUM_FRACS)
	    fatal_error(
	      "Generated more Farey fractions than expected. Exiting");
	  fracs[ptr].val=((double) num)/((double) den);
	  fracs[ptr].num=num;
	  fracs[ptr++].den=den;
	};


  cout << "Found " << ptr << " distinct Farey fractions." << endl;

  if(ptr!=NUM_FRACS)
    fatal_error("Didn't find enough fractions. Exiting.");

  cout << "Sorting." << endl;

  qsort(fracs,ptr,sizeof(frac),frac_comp);
  /*
  cout << fracs[0].val << endl;
  for(ptr2=1;ptr2<ptr;ptr2=ptr2<<1)
    printf("%20.18f\n",fracs[ptr2].val);
  */

  return(SUCCESS);
};

/* find_farey: returns the index of the lowest Farey fraction that
   equals or exceeds x. If x>1 then not defined, but probably returns
   the index of 1/1. This routine will need thinking about when we do
   an interval version. Probably just work on (say) left points. */
int find_farey(double x)
{
  int low,high,mid;

  if(x<=1.0/Q)
    return(NUM_FRACS-1);
  low=NUM_FRACS-1;
  high=0;
  mid=(low+high)>>1;
  while(low-mid>1)
    {
      if(fracs[mid].val>x)
	high=mid;
      else
	low=mid;
      mid=(low+high)>>1;
    };
  if(fracs[mid].val>=x)
    return(mid);
  return(high);
};

/*
dcomplex zeta(dcomplex s)
  // must do this properly. since t<=100 Euler Maclaurin summation?
{
  //  cout << s << endl;
  if(s==dcomplex(0.5))
    return(dcomplex(-1.4603545088095868129));
  if(s==dcomplex(1.5))
    return(dcomplex(2.6123753486854883433));
  if(s==dcomplex(2.5))
    return(dcomplex(1.3414872572509171798));
  if(s==dcomplex(3.5))
    return(dcomplex(1.1267338673170566464));
  if(s==dcomplex(4.5))
    return(dcomplex(1.0547075107614542640));
  if(s==dcomplex(5.5))
    return(dcomplex(1.0252045799546856946));
  if(s==dcomplex(6.5))
    return(dcomplex(1.0120058998885247961));
  if(s==dcomplex(7.5))
    return(dcomplex(1.0058267275365228077));
  if(s==dcomplex(8.5))
    return(dcomplex(1.0028592508824156277));
  if(s==dcomplex(9.5))
    return(dcomplex(1.0014125906121736623));
  fatal_error("Cant to generic zeta's yet!. Exiting.");
  return(dcomplex());
};
*/

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

void make_r_ns (dcomplex *new_r_n_vec, 
		dcomplex *old_r_n_vec,
		dcomplex *s_array,
		unsigned int s_array_size,
		unsigned int n, 
		double delta)
{
  int i,j;
  double delta_n[n];

  /* 
  cout << "Entered make r_n with n=" << n << " and delta =" << delta << endl;
  for(i=0;i<n;i++)
    cout << old_r_n_vec[i] << endl;
  cout << "***" << endl;
  */  

  delta_n[0]=1;
  for(i=1;i<n;i++)
    delta_n[i]=delta_n[i-1]*delta;
  for(i=0;i<n;i++)
    {
      new_r_n_vec[i]=old_r_n_vec[i];
      for(j=0;j<n-i-1;j++)
	new_r_n_vec[i]+=s_array[i*s_array_size+j]*delta_n[j+1]*old_r_n_vec[j+1];
    };
  /*
  for(i=0;i<n;i++)
    cout << new_r_n_vec[i] << endl;
  cout << "Exiting make_r_ns." << endl;
  */
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

int calc_hurwitz1(dcomplex s, unsigned int big_r_n_terms, double gap, 
		  unsigned int little_r_n_terms,unsigned int N,
		  dcomplex *zeta_vals)
{
  double delta,last_val;
  int i,j,farey_ptr;
  dcomplex s0,*this_r_n_vec,*last_r_n_vec;
  dcomplex *this_big_r_n_vec,*last_big_r_n_vec,*s_array;

  s0=s;

  if(!(s_array=new dcomplex[(big_r_n_terms-1)*(big_r_n_terms-1)]))
     fatal_error("Can't allocate memory for s_array. Exiting.");
  create_s_array(s,s_array,big_r_n_terms-1);


  
  for(j=0;j<big_r_n_terms-1;j++)
    for(i=0;i<=j;i++)
      cout << "s_array[" << i << "," << 
	j << "] is " << s_array[i*(big_r_n_terms-1)+j] << endl;
  

  if(!(last_big_r_n_vec=new dcomplex[big_r_n_terms]))
    fatal_error("Couldn't allocate memory for r_n terms. Exiting");
  if(!(this_big_r_n_vec=new dcomplex[big_r_n_terms]))
    fatal_error("Couldn't allocate memory for r_n terms. Exiting");
  if(!(last_r_n_vec=new dcomplex[little_r_n_terms]))
    fatal_error("Couldn't allocate memory for r_n terms. Exiting");
  if(!(this_r_n_vec=new dcomplex[little_r_n_terms]))
    fatal_error("Couldn't allocate memory for r_n terms. Exiting");

  for(i=0;i<big_r_n_terms;i++)
    {
      last_big_r_n_vec[i]=zeta_vals[i]-s_n(s0,1.0,N);
      //cout << s0 << " " << zeta_vals[i] << " " << s_n(s0,1.0,N) << endl;
      //cout << "r_n[" << i << "] is " << last_big_r_n_vec[i] << endl;
      real(s0)++;
    };
  fracs[0].z=zeta_vals[0];

  for(i=0;i<little_r_n_terms;i++)
    last_r_n_vec[i]=last_big_r_n_vec[i];

  i=0;
  farey_ptr=1;
  last_val=1.0;

  while(farey_ptr<NUM_FRACS)
    {
      delta=last_val-fracs[farey_ptr].val;
      if(delta>gap)
	{
	  i++;
	  while(delta>gap)
	    {
	      /* next fraction too far away so generate big_r_n_terms
		 r_n's from last big one */
	      
	      make_r_ns(this_big_r_n_vec,
			last_big_r_n_vec,
			s_array,big_r_n_terms-1,
			big_r_n_terms,
			gap);
	      swap(this_big_r_n_vec,last_big_r_n_vec);
	      delta=delta-gap;
	    };
	  if(delta!=0.0)
	    make_r_ns(this_big_r_n_vec,
		      last_big_r_n_vec,
		      s_array,big_r_n_terms-1,
		      big_r_n_terms,
		      delta);

	  fracs[farey_ptr].z=s_n(s,fracs[farey_ptr].val,N)+this_big_r_n_vec[0];

	  /*
	  printf("%20.18f\n",real(fracs[farey_ptr].z));
	  exit(0);
	  */

	  swap(this_big_r_n_vec,last_big_r_n_vec);


	  for(j=0;j<little_r_n_terms;j++)
	    last_r_n_vec[j]=last_big_r_n_vec[j];

	  last_val=fracs[farey_ptr++].val;
	}
      else
	{
	  /* generate little_r_n_terms r_n's from previous fraction

	  cout << "Before make_r_n " << endl;
	  cout << "this_r_n_vec" << endl;
	  for(j=0;j<little_r_n_terms;j++)
	    cout << this_r_n_vec[j] << " ";
	  cout << endl;
	  cout << "last_r_n_vec" << endl;
	  for(j=0;j<little_r_n_terms;j++)
	    cout << last_r_n_vec[j] << " ";
	  cout << endl;
	  */

	  make_r_ns(this_r_n_vec,
		    last_r_n_vec,
		    s_array,big_r_n_terms-1,
		    little_r_n_terms,
		    delta);

	  /*
	  cout << "After make_r_n " << endl;
	  cout << "this_r_n_vec" << endl;
	  for(j=0;j<little_r_n_terms;j++)
	    cout << this_r_n_vec[j] << " ";
	  cout << endl;
	  cout << "last_r_n_vec" << endl;
	  for(j=0;j<little_r_n_terms;j++)
	    cout << last_r_n_vec[j] << " ";
	  cout << endl;
	  */


	  fracs[farey_ptr].z=s_n(s,fracs[farey_ptr].val,N)+this_r_n_vec[0];

	  /*
	  printf("Z= %20.18f\n",real(fracs[farey_ptr].z));
	  */

	  swap(this_r_n_vec,last_r_n_vec);

	  /*
	  cout << "After swap " << endl;
	  cout << "this_r_n_vec" << endl;
	  for(j=0;j<little_r_n_terms;j++)
	    cout << this_r_n_vec[j] << " ";
	  cout << endl;
	  cout << "last_r_n_vec" << endl;
	  for(j=0;j<little_r_n_terms;j++)
	    cout << last_r_n_vec[j] << " ";
	  cout << endl;
	  exit(0);
	  */

	  last_val=fracs[farey_ptr++].val;
	};
    };

  cout << "I had to use " << big_r_n_terms << " terms " << i << " times."  << endl;

  for(i=0;i<30;i++)
    {
      cout << "Z(" << s << "," << fracs[i].num << "/" << fracs[i].den << "=";
      printf("%20.18f\n",real(fracs[i].z));
    };
  for(i=NUM_FRACS-11;i<NUM_FRACS;i++)
    {
      cout << "Z(" << s << "," << fracs[i].num << "/" << fracs[i].den << "=";
      printf("%20.18f\n",real(fracs[i].z));
    };

  return(SUCCESS);
};

int main(int argc, char **argv)
{
  double gap,zeta_val;
  dcomplex s,*zeta_vals;
  unsigned int big_r_n_terms,little_r_n_terms,N,i;
  char *fname;
  ifstream zeta_file;

  if(argc!=6)
    print_usage();
  gap=atof(argv[2]);
  if((gap<=0.0)||(gap>=1.0))
     print_usage();
  big_r_n_terms=atoi(argv[1]);
  if(big_r_n_terms<=0)
    print_usage();
  little_r_n_terms =atoi(argv[3]);
  if(little_r_n_terms<=0)
    print_usage();
  if(little_r_n_terms>big_r_n_terms)
    print_usage();
  N=atoi(argv[4]);
  if(N<0)
    print_usage();
  (zeta_file.open(argv[5]));
  if(!zeta_file.is_open())
    fatal_error("Couldn't open zeta file. Exiting.\n");
  if(!(zeta_vals=new dcomplex[big_r_n_terms]))
    fatal_error("Couldn't allocate memory for zeta_vals. Exiting.\n");
  zeta_file >> zeta_val;
  real(s)=zeta_val;
      if(zeta_file.eof())
	fatal_error("Insufficient zeta values given. Exiting.\n");
  zeta_file >> zeta_val;
  imag(s)=zeta_val;
  for(i=0;i<big_r_n_terms;i++)
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


  printf("terms 1 = %20d\ngap     = %20.18f\nterms 2 = %20d\nN       = %20d\n",
	 big_r_n_terms,gap,little_r_n_terms,N);


  if(!make_fracs())
    fatal_error("Error making Farey Fractions. Exiting.");


  cout << "Entering calc_hurwitz1." << endl;
  if(!calc_hurwitz1(s,big_r_n_terms,gap,
		    little_r_n_terms,N,zeta_vals))
    fatal_error("Error running Hurwitz routine. Exiting.");


  return(SUCCESS);
};

