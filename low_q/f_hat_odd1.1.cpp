/*
  File: f_hat_odd.cpp

  Created: 10th Feb 2010

  Version: <v> = 1.1

  Last Modified: 8th March 2010

  1.0 Initial implementation
  1.1 Output now pre-sorted for f_odd.cpp

  Dialect: C++


  Implementation notes: 

  Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

  Uses DFT to create F_hat_odd(x) for x = [2*Pi*n0/B,2*Pi*n1/B) for all
  primitive characters mod q.

  By: DJ Platt
  Bristol University

  Copyright 2010.

  This work is funded by the UK ESPRC. */

#define VERSION "1.1"
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
  printf("Usage: f_hat_odd%s (ifname) (facs_file) (ofname)\n",VERSION);
  printf("  (ifname)    - output from F_hat_odd_terms.c\n");
  printf("  (facs_file) - file with factors data.\n");
  printf("  (ofname)    - input to f.cpp (via spec file)\n");
  exit(1);
}


void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};

unsigned int conj_j (unsigned int j, unsigned int num_chi, unsigned int q)
{
	if(q&7) // q not divisible by 8
	  return(num_chi-j-1);
	if(j<(num_chi>>1))
	  return((num_chi>>1)-j-1);
	return(num_chi+(num_chi>>1)-j-1);
}


int_double d_log_pi,d_log2,d_small;
int_complex c_small;

#define MAX_LOGS (10000)

int_double logns[MAX_LOGS];

void setup()
{
  d_log_pi=log(d_pi);
  d_log2=log(int_double(2));

  d_small=exp(int_double(LN_SMALL));
  d_small.left=-d_small.right;
  c_small=int_complex(d_small,d_small);
  logns[1]=d_zero;
  for(int n=2;n<MAX_LOGS;n++)
    logns[n]=log(int_double(n));
}

#define TARGET_ACC (40)
// error from truncating sum of F_hat to M terms
inline int_complex F_hat_o_err(const int_double &x,const unsigned int q, const int_double &logq, unsigned int *M, const int_double &delta, const int_double &logdelta,int_double &pi_sindelta_by_q)
{
  int_double a,b;
  int_double lambda=exp(x+x)*pi_sindelta_by_q;
  b=d_log2-logq*0.75+x*1.5;
  //print_int_double_str("b=",b);
  //print_int_double_str("lambda=",lambda);
  if(b.left>TARGET_ACC)
    M[0]=1;
  else
    {
      a=sqrt((TARGET_ACC-b)/lambda);
      M[0]=(unsigned int) a.left;
      M[0]++;
    }
  //print_int_double_str("(m+1)^2*lambda=",lambda*(M[0]+1)*(M[0]+1));
  a=exp(b-lambda*(M[0]+1)*(M[0]+1))*(M[0]+1);
  //print_int_double_str("a=",a);
  int_double r=exp(-lambda*(3+2*M[0]))*(M[0]+2);
  r/=M[0]+1;
  r=d_one-r;
  //print_int_double_str("r=",r);
  b=a/r;
  b.left=b.right;
  //print_int_double_str("Err_o=",b);
  return(int_complex(b,b));
}

inline int_complex F_hat_o_term(const int_complex &outer_exp, const int_complex &exp_2_u_x, const unsigned int n)
{
  unsigned int n_sqr=n*n;
  int_double logn;
  if(n<MAX_LOGS)
    logn=logns[n];
  else
    logn=log(int_double(n));
  //printf("in F_hat_o_term with n=%d\n",n);
  //print_int_complex_str("outer_exp=",outer_exp);
  //print_int_complex_str("exp(2u(x))",exp_2_u_x);
  //exit(0);
  int_complex inner_exp=outer_exp-exp_2_u_x*n_sqr+logn;
  //print_int_complex_str("inner_exp=",inner_exp);
  //printf("abs_error=%d %d\n",-abs_error(inner_exp.real),-abs_error(inner_exp.imag));
  if(-inner_exp.real.right<LN_SMALL)
    return(c_small);
  //print_int_complex_str("F_hat_o_term returning",exp(inner_exp));exit(0);
  return(exp(inner_exp));
}

unsigned int *prims,*indices;
int_complex *f_omegas,*zs;
bool *neg_ones;

void make_l_odd(unsigned int q, 
		unsigned int num_s, 
		factor *factors,
		unsigned int *q_n,
		unsigned int *a_n,
		unsigned int *offset_vec,
		int_complex *omegas,
		const double B,
		unsigned int n0,
		const double &eta,
		FILE *infile)
{
  unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS],M;
  int dims[MAX_FACS];
  int no_dims;
  bool power_2,primitive,neg_one,first=true;
  unsigned int pr=factors[q].pr,prim_ptr,neg_ptr;
  unsigned int i,j,k,offset,s_done;
  int_double x,log2,err;
  int_double two_pi_by_B,delta,logdelta,pi_sindelta_by_q,logq,eta_pi_by_4;
  
  two_pi_by_B=read_int_double(infile);
  delta=read_int_double(infile);
  logdelta=read_int_double(infile);
  eta_pi_by_4=read_int_double(infile);
  pi_sindelta_by_q=read_int_double(infile);
  log2=read_int_double(infile);
  logq=read_int_double(infile);

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
	  fread(&M,sizeof(unsigned int),1,infile);
	  err=read_int_double(infile);
	  for(i=1;i<min(q,M+1);i++)
	    if(co_prime(i,q))
	      {
		Fft_vec[q_n[i]]=read_int_complex(infile);
		//printf("Pre FFT abs_error  = %d,%d\n",-abs_error(Fft_vec[q_n[i]].real),-abs_error(Fft_vec[q_n[i]].real));
	      }
	  for(;i<q;i++)
	    if(co_prime(i,q))
	      Fft_vec[q_n[i]]=c_zero;
	  bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare); 
	  prim_ptr=0;
	  neg_ptr=0;
	  for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
	    if(prim_p(i,factors[q].primes[0]))
	      {
		neg_one=(i&1);
		neg_ones[prim_ptr]=neg_one;
		if(neg_one)
		  {
		    if(s_done==0) // first time through, so finish off omega values
		      {
			f_omegas[prim_ptr]=finish_omega(omegas[i],true);
			indices[prim_ptr]=i;
			prims[prim_ptr]=neg_ptr;
		      }
		    zs[neg_ptr*num_s+s_done]=Fft_vec[i]*f_omegas[prim_ptr]+int_complex(err,err);
		    neg_ptr++;
		  }
		prim_ptr++;
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
	  fread(&M,sizeof(unsigned int),1,infile);
	  err=read_int_double(infile);
	  for(i=0,j=1;j<min(q,M+1);j++)
	    if(co_prime(j,q))
		Fft_vec[offset_vec[i++]]=read_int_complex(infile);
	  for(;j<q;j++)
	    if(co_prime(j,q))
	      Fft_vec[offset_vec[i++]]=c_zero;
	  
	  do_nd_fft(Fft_vec,no_dims,dims,phi_q);

	  for(i=0;i<factors[q].num_facs;i++)
	    coords[i]=0;
	  prim_ptr=0;
	  neg_ptr=0;
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
		  neg_ones[prim_ptr]=neg_one;
		  if(neg_one)
		    {
		      prims[prim_ptr]=neg_ptr;
		      if(s_done==0) // first time for this im_s and chi
			{
			  //printf("index %i is primitive\n",i);
			  f_omegas[prim_ptr]=finish_omega(omegas[i],neg_one);
			  indices[prim_ptr]=i;
			}
		      zs[neg_ptr*num_s+s_done]=Fft_vec[i]*f_omegas[prim_ptr]+int_complex(err,err);
		      neg_ptr++;
		    }
		  prim_ptr++;
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

  int_complex gam,*omegas;
  unsigned int *a_n,*q_n,*offset_vec,num_fracs,n_prims,i;
  unsigned int N,q,num_s,phi_q,n0,n1,prim,prim1;
  FILE *infile,*out_file;
  ifstream facs_file;
  factor *factors;
  clock_t no_clicks;
  double eta_odd,B;
  bool real_p;

  no_clicks=clock(); // start timing
  setup();
  //printf("At program start.\n");

  //printf("argc=%d\n",argc);
  if(argc!=4)
    print_usage();

  
  infile=fopen(argv[1],"rb");
  if(!infile)
    fatal_error("Couldn't open infile. Exiting.\n");

  facs_file.open(argv[2]);
  if(!facs_file.is_open())
    fatal_error("Couldnt open factors file. Exiting.\n");

  
  out_file=fopen(argv[3],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");

  fread(&q,sizeof(unsigned int),1,infile);  

  if(!(factors=(factor *) calloc(q+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  if(!read_factors(factors,q,facs_file))
    fatal_error("Error reading factor file. Exiting.\n");

  facs_file.close();

  fread(&eta_odd,sizeof(double),1,infile);  
  fread(&n0,sizeof(unsigned int),1,infile);  
  fread(&n1,sizeof(unsigned int),1,infile);  
  fread(&N,sizeof(unsigned int),1,infile);  

  B=N*one_over_A;
  //printf("B=%f\n",B);
  phi_q=factors[q].phi;

  num_s=n1-n0;
  //printf("num_s=%d\n",num_s);
  n_prims=num_prims(q,factors); // no of prim characters mod q
  if(!(prims=(unsigned int *) malloc(n_prims*sizeof(unsigned int))))
    fatal_error("Couldn't allocate memory for prims.\n");
  if(!(indices=(unsigned int *) malloc(n_prims*sizeof(unsigned int))))
    fatal_error("Couldn't allocate memory for indices.\n");
  if(!(f_omegas=(int_complex *) _aligned_malloc(n_prims*sizeof(int_complex),16)))
    fatal_error("Couldn't allocate memory for omegas.\n");
  if(!(neg_ones=(bool *) malloc(n_prims*sizeof(bool))))
    fatal_error("Couldn't allocate memory for neg_ones.\n");
  if(!(zs=(int_complex *) _aligned_malloc(((n_prims+1)/2)*num_s*sizeof(int_complex),16)))
    fatal_error("Couldn't allocate memory for zs.\n");

  //printf("num_prims=%d\n",n_prims);

  fwrite(&q,sizeof(unsigned int),1,out_file);
  fwrite(&n_prims,sizeof(unsigned int),1,out_file);
  fwrite(&N,sizeof(int),1,out_file);
  fwrite(&eta_odd,sizeof(double),1,out_file);
  fwrite(&num_s,sizeof(int),1,out_file);
  fwrite(&n0,sizeof(int),1,out_file);


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
  make_l_odd(q,num_s,factors,q_n,offset_vec,a_n,omegas,B,n0,eta_odd,infile);

  for(prim=0;prim<n_prims;prim++)
    {
      prim1=conj_j(prim,n_prims,q);
      real_p=(prim1==prim);
      if(prim1<prim)
	continue;
      //printf("Saving q:%d ind1:%d ind2:%d\n",q,prim,prim1);
      fwrite(&neg_ones[prim],sizeof(bool),1,out_file);
      fwrite(&real_p,sizeof(bool),1,out_file);
      if(!neg_ones[prim])
	continue;
      //printf("prims[%d]=%d\n",prim,prims[prim]);
      fwrite(&indices[prim],sizeof(unsigned int),1,out_file);
      fwrite(&f_omegas[prim],sizeof(int_complex),1,out_file);
      //print_int_complex_str("saving",zs[prims[prim]*num_s+1]);
      fwrite(&zs[prims[prim]*num_s],sizeof(int_complex),num_s,out_file);
      if(!real_p) // it has a conjugate so save it
	{
	  //printf("Saving composite ind:%d\n",prim1);
	  //printf("prims[%d]=%d\n",prim1,prims[prim1]);

	  fwrite(&indices[prim1],sizeof(unsigned int),1,out_file);
	  fwrite(&f_omegas[prim1],sizeof(int_complex),1,out_file);
	  //print_int_complex_str("saving",zs[prims[prim1]*num_s+1]);
	  fwrite(&zs[prims[prim1]*num_s],sizeof(int_complex),num_s,out_file);
	}
      //else
      //printf("\n");
    }


  fclose(out_file);
  fclose(infile);

  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  printf("f_hat_odd1.1 successful completion on q=%d n=[%d,%d)\n",q,n0,n1);
  return(0);
}
