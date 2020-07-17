/*

File: l-func-gmp.cpp

Created: 30th September 2008

Version: <v> = 1.0

Last Modified: 30th September 2008

Dialect: C++

Requires: gmpxx gmp

Implementation notes: 

Build instructions: g++ -ol-func-gmp l-func-gmp.cpp -O3 -lgmpxx -lgmp

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "1.0"
#define _CRT_SECURE_NO_WARNINGS  // shut Microsoft up

#ifdef __GNUC__
#define _aligned_malloc(x,y) malloc(x)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <complex>
#include "/cygdrive/c/cygwin/home/dave/mpfrcpp/mpreal.h"

#define PREC 100
#define debug printf("got to line %d\n",__LINE__);

using namespace std;
using namespace mpfr;


typedef complex<mpreal> dcomplex;
typedef dcomplex int_complex;

typedef mpreal int_double;

void print_int_complex(dcomplex z)
{printf("%20.18e+%20.18ei",double(real(z)),double(imag(z)));}

void print_interval(int_double x)
{printf("%20.18e",double(x));}

#define double(x) mpfr_get_d(x,GMP_RNDN)



// use sincos opcode
dcomplex pow (int_double x, dcomplex z)
{
	int_double tlnx,xtosigma;

	tlnx=log(x)*imag(z);
	xtosigma=pow(x,real(z));
	return(dcomplex(xtosigma*cos(tlnx),xtosigma*sin(tlnx)));
}

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
}

inline int gcd (unsigned int a, unsigned int b)
  /* Euclid algorithm gcd */
  /* would binary GCD be faster? */
{
  unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    }
  return(b);
}

#define co_prime(a,b) (gcd(a,b)==1)

void create_s_array(dcomplex &s,dcomplex *s_array, unsigned int n, int_double gap)
{
  unsigned int i,del;
  int_double fact=2.0;

  s_array[0]=s*gap;
  for(i=1;i<n;i++)
      s_array[i*n+i]=dcomplex(real(s)+(int_double)i,imag(s))*gap;

  del=1;
  while(del<n)
    {
      for(i=0;i<n-del;i++)
		  s_array[i*n+i+del]=s_array[i*n+i]*s_array[i*n+n+i+del]/fact;
      del++;
      fact=fact+1.0;
    }
}

void print_r_n_vec(dcomplex *r_n_vec,unsigned int n)
{
  unsigned int i;

  for(i=0;i<n;i++)
    {
      cout << "r_n_vec[" << i << "] is ";
      print_int_complex(r_n_vec[i]);
      cout << endl;
    }
}


// returns approximation to Zeta(s,alpha)-1/alpha
// not used in favour of (faster?) inline version
// within make_rns
dcomplex zeta(dcomplex s,int_double alpha, unsigned int ZETA_TERMS)
{
	dcomplex tot=dcomplex(0.0,0.0);
	dcomplex delta;
	unsigned int i;

	for(i=1;i<=ZETA_TERMS;i++)
		tot+=pow(i+alpha,-s);
	return(tot);
}

// put r_n(s+n,beta-alpha) into this_r_n_vec[0..n-1]
// where r_n(s+n,beta) is in last_r_n_vec[0..n-1]
// use Taylor expansion for M terms, series for the rest.
// there must be at least one "rest"
void make_r_ns (dcomplex *s_array,
		unsigned int n,
		unsigned int M,
		dcomplex *this_r_n_vec,
		dcomplex *last_r_n_vec,
		dcomplex s,
		int_double alpha,
		dcomplex *zeta_terms,
		unsigned int ZETA_TERMS)
{
  unsigned int i,j;
  dcomplex s_m;
 
  // use Taylor for the terms 0..M-1
  // s_array contains s(s+1)(s+2)...(s+n-1)*gap^n/n!
  for(i=0;i<M;i++)
    {
      this_r_n_vec[i]=last_r_n_vec[i];
      for(j=0;j<n-i-1;j++)
		  this_r_n_vec[i]+=s_array[i*(n-1)+j+i]*last_r_n_vec[j+i+1];
    }

  // use series for the rest
  // r_n(s,alpha)=sum n=1..ZETA_TERMS (n+alpha)^(-s)
  s_m=s+dcomplex(M,0);
  for(i=0;i<ZETA_TERMS;i++)
	  zeta_terms[i]=pow(1+i+alpha,-s_m);
  this_r_n_vec[M]=zeta_terms[0];
  for(i=1;i<ZETA_TERMS;i++)
	  this_r_n_vec[M]+=zeta_terms[i];

  for(i=M+1;i<n;i++)
  {
	  for(j=0;j<ZETA_TERMS;j++)
		  zeta_terms[j]/=((int_double)j+1.0+alpha);
	  this_r_n_vec[i]=zeta_terms[0];
	  for(j=1;j<ZETA_TERMS;j++)
		  this_r_n_vec[i]+=zeta_terms[j];
  }
  /* the above does this, but quicker
     is it less accurate?
    for(i=M;i<n;i++)
	     this_r_n_vec[i]=zeta(s+dcomplex(i,0),alpha,ZETA_TERMS);
*/
}

// use r_n_terms of the r_n_vec for alpha to calculate H(s,alpha-delta*gap)
dcomplex calc_r_n(unsigned int r_n_terms, dcomplex *s_array, const int_double &delta,
		  dcomplex *r_n_vec)
{
  unsigned int i;
  dcomplex res,term;

  res=r_n_vec[0];

  for(i=0;i<r_n_terms-1;)
    {
      term=s_array[i]*pow(delta,(int_double) i+1);
      i++;
      term*=r_n_vec[i];
      res+=term;
    }
  return(res);
}

// for a given s, calculate all the values of
// H(s,a/q)*q^(-s) where (a,q)=1 and a<q
int calc_hurwitz1(dcomplex &s, unsigned int r_n_terms,
				  unsigned int N, unsigned int M,
				  int_double gap,dcomplex *s_array,
				  dcomplex *r_n_vals,unsigned int no_gaps,
				  unsigned int q_start,
				  unsigned int q_end, FILE *out_file,
				  dcomplex *zeta_terms,
				  unsigned int ZETA_TERMS)
{
  dcomplex out_val,minus_s=-s;
  complex<double> z;
  unsigned int i,gap_ptr,q,num,num_fracs;
  int_double x,frac,steps;
  dcomplex q_minus_s,num_minus_s;

  create_s_array(s,s_array,r_n_terms-1,gap);

  // the first row of r_n_vals is set to r_n[i]=Zeta(s+i)-1 by main
  for(i=1;i<no_gaps;i++)
    make_r_ns(s_array,r_n_terms,M,&r_n_vals[i*r_n_terms],
	      &r_n_vals[(i-1)*r_n_terms],s,1.0-i*gap,zeta_terms,ZETA_TERMS);

  num_fracs=0;
  for(q=q_start;q<=q_end;q++)
  {
	  if((q&3)==2)
		  continue;      // if q = 2 mod 4 then no primitive characters
	  q_minus_s=pow(q,minus_s);
	  for(num=1;num<q;num++)
		  if(co_prime(num,q))  // q>num so this way round?
		  {
			  num_fracs++;
			  num_minus_s=pow(num,minus_s);
			  frac=int_double(num)/q;
			  steps=ceil(num/(gap*q));
			  gap_ptr=r_n_terms*(no_gaps-steps);
			  x=steps-frac/gap;
			  out_val=calc_r_n(N,s_array,x,&r_n_vals[gap_ptr]);
			  out_val*=q_minus_s;
			  out_val+=num_minus_s;
			  if(num==1)
			  {
				  print_int_complex(out_val);
			      printf("\n");
			  }
			  z=complex<double>(double(real(out_val)),double(imag(out_val)));
			  fwrite(&z,sizeof(complex<double>),1,out_file);
		  }
  }
  return(num_fracs);
}

int main(int argc, char **argv)
{
	double ds,ds1,gap;
	double x;
	dcomplex s,*s_array,*r_n_vals,*zeta_terms;
	int r_n_terms,N,M,i,j,no_gaps,q_start,q_end,num_s,num_fracs,num_zeta,
		ZETA_TERMS;
	int_double s_delta,*gap_ns;
	int_double im_s;
	ifstream zeta_file;   // contains zeta(s)-1 for s=0.5,1.5....99.5,0.5+dsI....,0.5+2dsI.....
	FILE *out_file;
	fpos_t out_file_ptr;
	clock_t no_clicks; // to collect timings.

	mpfr_set_default_prec(PREC);
	no_clicks=clock(); // start timing
	if(argc!=10)       // check correct number of arguments
		print_usage();
	q_start=atoi(argv[1]);
	if(q_start<3)
		print_usage();
	q_end=atoi(argv[2]);
	if(q_end<q_start)
		print_usage();
	r_n_terms=atoi(argv[3]);
	if(r_n_terms<=0) // number of lattice points per step
		print_usage();
	N=atoi(argv[4]); // how many lattice points used to calculate intermediates
	if((N<=0)||(N>r_n_terms))
		print_usage();
	M=atoi(argv[5]); // how many lattice points to calculate using Taylor
	// rather than simple series.
	if((M<=0)||(M>=r_n_terms))
		print_usage();
	zeta_file.open(argv[6]); 
	if(!zeta_file.is_open())
		fatal_error("Couldn't open ifile. Exiting.\n");
	out_file=fopen(argv[7],"wb");
	if(!out_file)
		fatal_error("Couldn't open ifile. Exiting.\n");
	gap=atof(argv[8]);
	if((gap<=0.0)||(gap>=1.0))
		print_usage();
	ZETA_TERMS=atoi(argv[9]);
	if(ZETA_TERMS<=0)
		print_usage();
	if(sizeof(long double)==8)
		printf("Warning, long double is merely double.\n");
	no_gaps=ceil(1.0/gap);
	cout << "Running with " << no_gaps << " gaps." << endl;
	cout << "         and " << ZETA_TERMS << " zeta terms." << endl;

	cout << "Allocating memory for s_array." << endl;
	if(!(s_array=new dcomplex[(r_n_terms-1)*(r_n_terms-1)]))
		fatal_error("Can't allocate memory for s_array. Exiting.");

	cout << "Allocating memory for R_N vectors." << endl;
	if(!(r_n_vals=new dcomplex[no_gaps*r_n_terms]))
		fatal_error("Can't allocate memory for r_n_vals. Exiting.");

	cout << "Allocating memory for gap_ns." << endl;
	if(!(gap_ns=new int_double[r_n_terms-1]))
		fatal_error("Can't allocate memory for gap_ns. Exiting.");

	cout << "Allocating memory for zeta_terms." << endl;
	if(!(zeta_terms=new dcomplex[ZETA_TERMS]))
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
	zeta_file >> ds;      // step between each s.
	s_delta=int_double(ds);    // assumes delta representable

	printf("Going to try to process %d values of s using step size\n    ",
		num_s);
	print_interval(s_delta);printf("\n");

	im_s=0.0;
	s=int_complex(0.5,im_s);
	for(i=0;i<num_s;i++)
	{
		x=(double) im_s;
		fwrite(&x,sizeof(double),1,out_file);	
		for(j=0;j<r_n_terms;j++)
		{
			zeta_file >> ds >> ds1;
			r_n_vals[j]=dcomplex(ds,ds1);
		}
		for(;j<num_zeta;j++)
			zeta_file >> ds >> ds;    // throw away extra ones
		num_fracs=calc_hurwitz1(s,r_n_terms,N,M,gap,s_array,
			r_n_vals,no_gaps,q_start,q_end,out_file,zeta_terms,
			ZETA_TERMS);
		if(!num_fracs)
			fatal_error("Error running Hurwitz routine. Exiting.");
		im_s+=s_delta;
		s=dcomplex(0.5,im_s);
	}

	fsetpos(out_file,&out_file_ptr);
	fwrite(&num_fracs,sizeof(int),1,out_file);
	zeta_file.close();
	fclose(out_file);

	printf("Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
	printf("Wrote %d fractions.\n",num_fracs);
	return(0);
}
