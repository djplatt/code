/*
  File: f_hat.cpp

  Created: 10th Feb 2010

  Version: <v> = 1.0

  Last Modified: 10th February 2010

  1.0 Initial implementation

  Dialect: C++


  Implementation notes: 

  Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

  Uses DFT to create F_hat(x) for x = [2*Pi*n0/B,2*Pi*n1/B) for all
  primitive characters mod q.

  These files are converted by low_q_conv.cpp, then read by f.cpp

  By: DJ Platt
  Bristol University

  Copyright 2010.

  This work is funded by the UK ESPRC. */

#define VERSION "1.0"
//#define PRINT
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft.h"
#include "f_defs.h"

void print_usage()
/* called when wrong arguments passed via command line */
{
  printf("Usage: int-l-func%s (q) (n0) (n1) (N) (eta) (facs_file) (ofname)\n",VERSION);
  printf("  (q)         - integer>=3,q&3!=2\n");
  printf("  (n0)        - integer>0.0\n");
  printf("  (n1)        - integer>n0\n");
  printf("  (N)         - integer>=2*n1\n");
  printf("  (eta_o)     - float in [0.0,1.0)\n");
  printf("  (eta_e)     - float in [0.0,1.0)\n");
  printf("  (facs_file) - file with factors data.\n");
  printf("  (ofname)    - output file.\n");
  exit(1);
}


void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};

int_double d_log_pi,d_log2,d_small;
int_complex c_small;

#define SMALL ((double) 1e-307)
#define ln_SMALL ((double) 706.89362355)  // exp(-ln_SMALL)<=SMALL
#define MAX_LOGS (10000)

int_double logns[MAX_LOGS];

void setup()
{
  d_log_pi=log(d_pi);
  d_log2=log(int_double(2));
  double small=1e-307; //NB using nextafter(0.0) is very slow. denormalisation problem?
  d_small=int_double(-small,small);
  c_small=int_complex(d_small,d_small);
  for(int n=1;n<MAX_LOGS;n++)
    logns[n]=log(int_double(n));
}

#define TARGET_ACC (40)
// error from truncating sum of F_hat to M terms
inline int_complex F_hat_e_err(const int_double &x,const unsigned int q, const int_double &logq, unsigned int *M, const int_double &delta, const int_double &logdelta, int_double &sindelta)
{
  int_double a,b;
  int_double lambda=d_pi*exp(x+x)*sindelta/q;
  b=d_log2-logq*0.25+x*0.5;
  if(b.left>TARGET_ACC)
    M[0]=1;
  else
    {
      a=sqrt((TARGET_ACC-b)/lambda);
      M[0]=(unsigned int) a.left;
      M[0]++;
    }
  //printf("M=%d\n",M[0]);
  a=exp(b-lambda*(M[0]+1)*(M[0]+1));
  int_double r=exp(-lambda*(3+2*M[0]));
  b=a/(d_one-r);
  b.left=b.right;
  //print_int_double_str("Err_e=",b);
  //exit(0);
  return(int_complex(b,b));
}

inline int_complex F_hat_o_err(const int_double &x,const unsigned int q, const int_double &logq, unsigned int *M, const int_double &delta, const int_double &logdelta,int_double &sindelta)
{
  int_double a,b;
  int_double lambda=d_pi*exp(x+x)*sindelta/q;
  b=d_log2-logq*0.75+x*1.5;
  if(b.left>TARGET_ACC)
    M[0]=1;
  else
    {
      a=sqrt((TARGET_ACC-b)/lambda);
      M[0]=(unsigned int) a.left;
      M[0]++;
    }
  //printf("M=%d\n",M[0]);
  a=exp(b-lambda*(M[0]+1)*(M[0]+1))*(M[0]+1);
  //print_int_double_str("a=",a);
  int_double r=exp(-lambda*(3+2*M[0]))*(M[0]+2);
  r/=M[0]+1;
  //print_int_double_str("r=",r);
  b=a/(d_one-r);
  b.left=b.right;
  //print_int_double_str("Err_o=",b);
  return(int_complex(b,b));
}

inline int_complex F_hat_e_term(const int_complex &outer_exp, const int_complex &exp_2_u_x, const unsigned int n)
{
  int_complex inner_exp=outer_exp-exp_2_u_x*n*n;
  if(inner_exp.real.right>ln_SMALL)
    return(c_small);
  return(exp(inner_exp));
}


inline int_complex F_hat_o_term(const int_complex &outer_exp, const int_complex &exp_2_u_x, const unsigned int n)
{
  int_double logn;
  if(n<MAX_LOGS)
    logn=logns[n];
  else
    logn=log(int_double(n));
  int_complex inner_exp=outer_exp-exp_2_u_x*n*n+logn;
  if(inner_exp.real.right>ln_SMALL)
    return(c_small);
  return(exp(inner_exp));
}

int_double two_pi_by_B,delta,logdelta,sindelta,logq;


void make_l_even(unsigned int q, 
		unsigned int num_s, 
		factor *factors,
		unsigned int *q_n,
		unsigned int *a_n,
		unsigned int *offset_vec,
		int_complex *omegas,
		FILE *out_file,
		const double B,
		unsigned int n0,
		const double &eta)
{
  unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS];
  int_complex z;
  int dims[MAX_FACS];
  int no_dims;
  bool power_2,primitive,neg_one;
  unsigned int pr=factors[q].pr,n_prims;
  unsigned int i,j,k,offset,s_done,M;
  int_double x;
  int_complex u_x,exp_2_u_x,outer_exp,err;

  //  printf("Processing Q=%d\n",q);


  j=0;
  for(i=1;i<q;i++)
    if(co_prime(i,q))
      a_n[j++]=i;
  if(pr)  // q has a primitive root so nice and easy 1-d FFT
    {
      fill_q_ns(q_n,pr,phi_q,q);
      init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
      // don't actually need to do this, covered by odd
      prep_omegas(omegas,q,phi_q,q_n,a,bs,b_star_conjs,conv_sizes[0],b_spare);
      s_done=0;
      while(s_done<num_s)
	{
	  x=two_pi_by_B*(n0+s_done);
	  //print_int_double_str("x=",x);
	  u_x=int_complex(int_double(x),d_pi*eta/4.0);
	  exp_2_u_x=exp(u_x*2)*d_pi/q; // Pi/q*exp(2*U(x))
	  outer_exp=u_x*0.5+d_log2-logq*0.25;
	  //print_int_complex_str("u(x)=",u_x);
	  //print_int_complex_str("outer_exp=",outer_exp);
	  err=F_hat_e_err(x,q,logq,&M,delta,logdelta,sindelta);
	  //printf("%d\n",M);
	  for(i=1;i<q;i++)
	    if(co_prime(i,q))
	      {
		Fft_vec[q_n[i]]=c_zero;
		for(k=i;k<=M;k+=q)
		  Fft_vec[q_n[i]]+=F_hat_e_term(outer_exp,exp_2_u_x,k);
	      }
	  /*
	  for(int n=0;n<phi_q;n++)
	    if(!contains_zero(Fft_vec[n]))
	      printf("before abs_error= %d %d\n",-abs_error(Fft_vec[n].real),-abs_error(Fft_vec[n].imag));
	  */
	  bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare); 
	  /*
	    for(int n=0;n<phi_q;n++)
	    if(!contains_zero(Fft_vec[n]))
	      printf("after abs_error= %d %d\n",-abs_error(Fft_vec[n].real),-abs_error(Fft_vec[n].imag));
	  */
	  for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
	    if(prim_p(i,factors[q].primes[0]))        
	      {
		neg_one=(i&1);
		if(!neg_one)
		  {
		    if(s_done==0) // first time through, so finish off omega values
		      {
			omegas[i]=finish_omega(omegas[i],false);
			//print_int_complex_str("omega=",omegas[i]);
			fwrite(&i,sizeof(unsigned int),1,out_file);
			fwrite(&omegas[i],sizeof(int_complex),1,out_file);
			fwrite(&neg_one,sizeof(bool),1,out_file);
		      }
		    z=Fft_vec[i]*omegas[i];
		    z+=err;
		    fwrite(&z,sizeof(int_complex),1,out_file);
		    //print_int_complex_str("z=",z);
		  }
	      }
	  s_done++;
	}
    }
  else // q doesn't have a primitive root
    {
      no_dims=factors[q].num_facs;
      fac=factors[q].facs[0];        // the first p^n
      power_2=(factors[fac].pr==0);  // no conductor => p=2, n>=3
      if(power_2)
	{
	  no_dims++;
	  fill_q_ns_2s(q_n,fac);    // use the {-1,1}X{5} trick
	  for(i=1;i<factors[q].num_facs;i++) // move all the dimensions up one
	    dims[i+1]=factors[factors[q].facs[i]].phi;
	  dims[1]=factors[q].facs[0]/4;
	  dims[0]=2;                         // slip in a two
	}
      else
	{
	  fill_q_ns(q_n,factors[fac].pr,factors[fac].phi,fac); // use the generator
	  for(i=0;i<factors[q].num_facs;i++)
	    dims[i]=factors[factors[q].facs[i]].phi;
	}
      offset=fac;
      for(i=1;i<factors[q].num_facs;i++)  // do the rest on the factors, all will have generators
	{
	  fac=factors[q].facs[i];      
	  pr=factors[fac].pr;
	  fill_q_ns(&q_n[offset],pr,factors[fac].phi,fac);  // use the generator
	  offset+=fac;
	}
      s_done=0;
      make_offsets(q,q_n,offset_vec,factors,a_n);    // reverse q_n so we know how to populate fft vector
      // don't need to do this either
      prep_omegas_nd(offset_vec,q,omegas,no_dims,dims,phi_q);
      while(s_done<num_s)
	{
	  x=two_pi_by_B*(n0+s_done);
	  u_x=int_complex(int_double(x),d_pi*eta/4.0);
	  exp_2_u_x=exp(u_x*2)*d_pi/q; // Pi/q*exp(2*U(x))
	  outer_exp=u_x*0.5+d_log2-logq*0.25;
	  err=F_hat_e_err(x,q,logq,&M,delta,logdelta,sindelta);
	  for(i=0,j=1;j<q;j++)
	    {
	      if(co_prime(j,q))
		{
		  Fft_vec[offset_vec[i]]=c_zero;
		  for(k=j;k<=M;k+=q)
		    Fft_vec[offset_vec[i]]+=F_hat_e_term(outer_exp,exp_2_u_x,k);
		  i++;
		}
	    }
	  do_nd_fft(Fft_vec,no_dims,dims,phi_q);

	  for(i=0;i<factors[q].num_facs;i++)
	    coords[i]=0;

	  for(i=0;i<phi_q;i++)
	    {
	      primitive=true;
	      for(j=0;j<factors[q].num_facs;j++)
		if(coords[j]%factors[q].primes[j]==0)
		  {
		    primitive=false;
		    break;
		  }

	      if(primitive)
		{
		  neg_one=neg_one_p(coords,factors[q].num_facs);
		  if(power_2&&(i<(phi_q>>1)))
		    neg_one=!neg_one;
		  if(!neg_one)
		    {
		      if(s_done==0) // first time for this im_s and chi
			{
			  omegas[i]=finish_omega(omegas[i],neg_one);
			  fwrite(&i,sizeof(unsigned int),1,out_file);
			  fwrite(&omegas[i],sizeof(int_complex),1,out_file);
			  fwrite(&neg_one,sizeof(bool),1,out_file);
			}
		      z=Fft_vec[i]*omegas[i];
		      z+=err;
		      fwrite(&z,sizeof(int_complex),1,out_file);
		    }
		}

	      j=factors[q].num_facs-1;
	      coords[j]++;
	      while(coords[j]==factors[factors[q].facs[j]].phi)
		{
		  coords[j]=0;
		  if(j==0)
		    break;
		  j--;
		  coords[j]++;
		}
	    }
	  s_done++;
	}
    }
}

unsigned int num_odds=0;

void make_l_odd(unsigned int q, 
		unsigned int num_s, 
		factor *factors,
		unsigned int *q_n,
		unsigned int *a_n,
		unsigned int *offset_vec,
		int_complex *omegas,
		FILE *out_file,
		const double B,
		unsigned int n0,
		const double &eta)
{
  unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS];
  int_complex z;
  int dims[MAX_FACS];
  int no_dims;
  bool power_2,primitive,neg_one;
  unsigned int pr=factors[q].pr,n_prims;
  unsigned int i,j,k,offset,s_done,M;
  two_pi_by_B=d_two_pi/B;
  delta=d_pi_2*(1.0-eta);
  logdelta=log(delta);
  int_double x;
  sin_cos(delta,&sindelta,&x);
  logq=log(int_double(q));
  int_complex u_x,exp_2_u_x,outer_exp,err;

  //  printf("Processing Q=%d\n",q);


  j=0;
  for(i=1;i<q;i++)
    if(co_prime(i,q))
      a_n[j++]=i;
  if(pr)  // q has a primitive root so nice and easy 1-d FFT
    {
      fill_q_ns(q_n,pr,phi_q,q);
      init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
      prep_omegas(omegas,q,phi_q,q_n,a,bs,b_star_conjs,conv_sizes[0],b_spare);
      s_done=0;
      while(s_done<num_s)
	{
	  x=two_pi_by_B*(n0+s_done);
	  //print_int_double_str("x=",x);
	  u_x=int_complex(int_double(x),d_pi*eta/4.0);
	  exp_2_u_x=exp(u_x*2)*d_pi/q; // Pi/q*exp(2*U(x))
	  outer_exp=u_x*1.5+d_log2-logq*0.75;
	  //print_int_complex_str("u(x)=",u_x);
	  //print_int_complex_str("outer_exp=",outer_exp);
	  err=F_hat_o_err(x,q,logq,&M,delta,logdelta,sindelta);
	  //printf("M=%d\n",M);
	  for(i=1;i<q;i++)
	    if(co_prime(i,q))
	      {
		Fft_vec[q_n[i]]=c_zero;
		for(k=i;k<=M;k+=q)
		  Fft_vec[q_n[i]]+=F_hat_o_term(outer_exp,exp_2_u_x,k);
	      }
	  //for(int n=0;n<phi_q;n++)
	  //  {
	  //    printf("abs error = %d,%d\n",-abs_error(Fft_vec[n].real),-abs_error(Fft_vec[n].imag));
	  //    print_int_complex_str("Pre-FFT",Fft_vec[n]);
	  //  }
	  bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare); 
	  /*
	    for(int n=0;n<phi_q;n++)
	    {
	      printf("abs error = %d,%d\n",-abs_error(Fft_vec[n].real),-abs_error(Fft_vec[n].imag));
	      //print_int_complex_str("Post-FFT",Fft_vec[n]);
	    }
	  exit(0);
	  */
	  for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
	    if(prim_p(i,factors[q].primes[0]))        
	      {
		neg_one=(i&1);
		if(neg_one)
		  {

		    if(s_done==0) // first time through, so finish off omega values
		      {
			num_odds++;
			omegas[i]=finish_omega(omegas[i],true);
			//print_int_complex_str("omega=",omegas[i]);
			fwrite(&i,sizeof(unsigned int),1,out_file);
			fwrite(&omegas[i],sizeof(int_complex),1,out_file);
			fwrite(&neg_one,sizeof(bool),1,out_file);
		      }
		    z=Fft_vec[i]*omegas[i];
		    z+=err;
		    fwrite(&z,sizeof(int_complex),1,out_file);
		    //if(i==1) print_int_complex_str("z=",z);
		  }
	      }
	  s_done++;
	}
    }
  else // q doesn't have a primitive root
    {
      no_dims=factors[q].num_facs;
      fac=factors[q].facs[0];        // the first p^n
      power_2=(factors[fac].pr==0);  // no conductor => p=2, n>=3
      if(power_2)
	{
	  no_dims++;
	  fill_q_ns_2s(q_n,fac);    // use the {-1,1}X{5} trick
	  for(i=1;i<factors[q].num_facs;i++) // move all the dimensions up one
	    dims[i+1]=factors[factors[q].facs[i]].phi;
	  dims[1]=factors[q].facs[0]/4;
	  dims[0]=2;                         // slip in a two
	}
      else
	{
	  fill_q_ns(q_n,factors[fac].pr,factors[fac].phi,fac); // use the generator
	  for(i=0;i<factors[q].num_facs;i++)
	    dims[i]=factors[factors[q].facs[i]].phi;
	}
      offset=fac;
      for(i=1;i<factors[q].num_facs;i++)  // do the rest on the factors, all will have generators
	{
	  fac=factors[q].facs[i];      
	  pr=factors[fac].pr;
	  fill_q_ns(&q_n[offset],pr,factors[fac].phi,fac);  // use the generator
	  offset+=fac;
	}
      s_done=0;
      make_offsets(q,q_n,offset_vec,factors,a_n);    // reverse q_n so we know how to populate fft vector
      prep_omegas_nd(offset_vec,q,omegas,no_dims,dims,phi_q);
      while(s_done<num_s)
	{
	  x=two_pi_by_B*(n0+s_done);
	  u_x=int_complex(int_double(x),d_pi*eta/4.0);
	  exp_2_u_x=exp(u_x*2)*d_pi/q; // Pi/q*exp(2*U(x))
	  outer_exp=u_x*1.5+d_log2-logq*0.75;
	  err=F_hat_o_err(x,q,logq,&M,delta,logdelta,sindelta);
	  for(i=0,j=1;j<q;j++)
	    {
	      if(co_prime(j,q))
		{
		  Fft_vec[offset_vec[i]]=c_zero;
		  for(k=j;k<=M;k+=q)
		    Fft_vec[offset_vec[i]]+=F_hat_o_term(outer_exp,exp_2_u_x,k);
		  i++;
		}
	    }
	  do_nd_fft(Fft_vec,no_dims,dims,phi_q);

	  for(i=0;i<factors[q].num_facs;i++)
	    coords[i]=0;

	  for(i=0;i<phi_q;i++)
	    {
	      primitive=true;
	      for(j=0;j<factors[q].num_facs;j++)
		if(coords[j]%factors[q].primes[j]==0)
		  {
		    primitive=false;
		    break;
		  }

	      if(primitive)
		{
		  neg_one=neg_one_p(coords,factors[q].num_facs);
		  if(power_2&&(i<(phi_q>>1)))
		    neg_one=!neg_one;
		  if(neg_one)
		    {
		      if(s_done==0) // first time for this im_s and chi
			{
			  num_odds++;
			  //printf("index %i is primitive\n",i);
			  omegas[i]=finish_omega(omegas[i],neg_one);
			  fwrite(&i,sizeof(unsigned int),1,out_file);
			  fwrite(&omegas[i],sizeof(int_complex),1,out_file);
			  fwrite(&neg_one,sizeof(bool),1,out_file);
			}
		      z=Fft_vec[i]*omegas[i];
		      z+=err;
		      fwrite(&z,sizeof(int_complex),1,out_file);
		    }
		}

	      j=factors[q].num_facs-1;
	      coords[j]++;
	      while(coords[j]==factors[factors[q].facs[j]].phi)
		{
		  coords[j]=0;
		  if(j==0)
		    break;
		  j--;
		  coords[j]++;
		}
	    }
	  s_done++;
	}
    }
}


int main(int argc, char **argv)
{

  _fpu_rndd();

  int_complex s,gam,*omegas;
  unsigned int *a_n,*q_n,*offset_vec,num_fracs,max_num_prims,n_prims,i;
  unsigned int N,q,num_s,phi_q,n0,n1;
  FILE *out_file;
  ifstream facs_file;
  factor *factors;
  clock_t no_clicks;
  double eta_odd,eta_even,B;

  no_clicks=clock(); // start timing
  setup();
  //printf("At program start.\n");

  //printf("argc=%d\n",argc);
  if(argc!=9)
    print_usage();
  
  q=atoi(argv[1]);
  if((q<3)||(q>MAX_Q)||((q&3)==2))
    print_usage();
  //printf("q=%d\n",q);

  n0=atoi(argv[2]);
  //printf("n0=%d\n",n0);
  if(n0<0)
    print_usage();

  n1=atoi(argv[3]);
  //printf("n1=%d\n",n1);
  if(n1<=n0)
    print_usage();

  N=atoi(argv[4]);
  //printf("N=%d\n",N);
  if((N<n1)||((N&63)!=0))
    print_usage();
  B=N*one_over_A;
  //printf("B=%f\n",B);

  eta_odd=atof(argv[5]);
  if((eta_odd<0.0)||(eta_odd>=1.0))
    print_usage();
  //printf("eta_odd=%f\n",eta_odd);

  eta_even=atof(argv[6]);
  if((eta_even<0.0)||(eta_even>=1.0))
    print_usage();
  //printf("eta_even=%f\n",eta_even);

  facs_file.open(argv[7]);
  if(!facs_file.is_open())
    fatal_error("Couldnt open factors file. Exiting.\n");

  if(!(factors=(factor *) calloc(q+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  if(!read_factors(factors,q,facs_file))
    fatal_error("Error reading factor file. Exiting.\n");

  facs_file.close();

  phi_q=factors[q].phi;
  
  out_file=fopen(argv[8],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");


  fwrite(&q,sizeof(unsigned int),1,out_file);
  num_s=n1-n0;
  //printf("num_s=%d\n",num_s);
  fwrite(&num_s,sizeof(int),1,out_file);
  n_prims=num_prims(q,factors); // no of prim characters mod q
  //printf("num_prims=%d\n",n_prims);
  fwrite(&n_prims,sizeof(unsigned int),1,out_file);
  fwrite(&n0,sizeof(int),1,out_file);
  fwrite(&N,sizeof(int),1,out_file);
  fwrite(&eta_odd,sizeof(double),1,out_file);
  fwrite(&eta_even,sizeof(double),1,out_file);
  long file_pos=ftell(out_file);
  fwrite(&num_odds,sizeof(unsigned int),1,out_file);


  if(!(omegas=(int_complex *) _aligned_malloc(phi_q*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  //printf("Initialising ws...\n");
  init_ws(_ws);

  //printf("Starting...\n");
  make_l_odd(q,num_s,factors,q_n,offset_vec,a_n,omegas,out_file,B,n0,eta_odd);
  make_l_even(q,num_s,factors,q_n,offset_vec,a_n,omegas,out_file,B,n0,eta_even);

  fseek(out_file,file_pos,SEEK_SET);
  fwrite(&num_odds,sizeof(unsigned int),1,out_file);

  fclose(out_file);

  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  printf("f_hat successful completion on q=%d n=[%d,%d)\n",q,n0,n1);
  return(0);
}
