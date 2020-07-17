/*
  File: safe.cpp

  Created: 10th Mar 2010

  Version: <v> = 1.0

  Last Modified: 10th March 2010

  1.0 Initial implementation

  Dialect: C++


  Implementation notes: 

  Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

  Use the smoothed approximate functional equation to calculate F(t) for large t

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



#include "../includes/int_double11.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft.h"
#include "f_defs.h"

void print_usage()
/* called when wrong arguments passed via command line */
{
  printf("Usage: safe%s (q) (t0) (N) (facs_file) (ofname)\n",VERSION);
  printf("  (q)    - modulus != 2 (4)\n");
  printf("  (t0)   - 1st t value\n");
  printf("  (n)    - number of t values to output.\n");
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

int_complex G(const int_complex &s, const int_complex &z)
{
  return(c_zero);
}


unsigned int *indices;
int_complex *f_omegas,*zs;
bool *neg_ones;

void safe_even()
{};

void safe_odd(unsigned int q, 
	      unsigned int num_s, 
	      factor *factors,
	      unsigned int *q_n,
	      unsigned int *a_n,
	      unsigned int *offset_vec,
	      int_complex *omegas,
	      const double &t0,
	      const int_double &pibyq)
{
  unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS],M;
  int dims[MAX_FACS];
  int no_dims;
  bool power_2,primitive,neg_one,first=true;
  unsigned int pr=factors[q].pr,prim_ptr,neg_ptr;
  unsigned int i,j,k,offset,s_done;
  int_complex s=int_complex(int_double(0.5),int_double(t0-one_over_A));
  int_complex delta,sby2,s1by2,saby2,sa1by2;
  sby2.real=int_double(0.25);
  s1by2.real=sby2.real;
  saby2.real=int_double(0.75);
  sa1by2.real=saby2.real;

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
	  s.imag+=one_over_A;
	  delta=calc_delta(s.imag);
	  M=calc_M(s.imag); // make M big enough to control truncation of sum
	  d2=delta*delta;
	  pid2byq=pibyq*d2;
	  piby2dq=pibyq/d2;
	  sby2.imag=s.imag/2;
	  s1by2.imag=sby2.imag;
	  saby2.imag=-sby2.imag;
	  sa1by2.imag=saby2.imag;
	  for(i=1;i<q;i++)
	    if(co_prime(i,q))
	      {
		Fft_vec[q_n[i]]=c_zero;
		for(j=i;j<=M;j+=q)
		  Fft_vec[q_n[i]]=delta*j*G(saby2,pid2byq*j*j);
	      }
	  bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare); 
	  prim_ptr=0;
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
		      }
		    zs[prim_ptr*num_s+s_done]=Fft_vec[i]*f_omegas[prim_ptr];
		  }
		prim_ptr++;
	      }
	  for(i=1;i<q;i++)
	    if(co_prime(i,q))
	      {
		Fft_vec[q_n[i]]=c_zero;
		for(j=i;j<=M;j+=q)
		  Fft_vec[q_n[q-i]]=j*G(sa1by2,pibyd2q*j*j);
	      }
	  bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare); 
	  prim_ptr=0;
	  for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
	    if(prim_p(i,factors[q].primes[0]))
	      {
		neg_one=(i&1);
		if(neg_one)
		  zs[prim_ptr*num_s+s_done]+=Fft_vec[i]*f_omegas[prim_ptr]/d2;
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
  int q1,num_s1;
  double t0;
  unsigned int *a_n,*q_n,*offset_vec,num_fracs,n_prims,i;
  unsigned int N,q,num_s,phi_q,n0,n1,prim,prim1;
  FILE *infile,*out_file;
  ifstream facs_file;
  factor *factors;
  clock_t no_clicks;
  double eta_odd,B;
  bool real_p;
  int_double pibyq;

  no_clicks=clock(); // start timing
  setup();
  //printf("At program start.\n");

  //printf("argc=%d\n",argc);
  if(argc!=6)
    print_usage();

  q1=atoi(argv[1]);
  if((q1<3)||((q1&3)==2))
    fatal_error("q must be >=3 and != 2 mod 4.Exiting.\n");
  q=q1;
  pibyq=d_pi/q;

  t0=atof(argv[2]);
  if(t0<0.0)
    fatal_error("t0 must be >=0.0. Exiting.\n");

  num_s1=atoi(argv[3]);
  if(num_s1<=0)
    fatal_error("N must be > 0. Exiting\n");
  num_s=num_s1;

  facs_file.open(argv[4]);
  if(!facs_file.is_open())
    fatal_error("Couldnt open factors file. Exiting.\n");

  
  out_file=fopen(argv[5],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");

  if(!(factors=(factor *) calloc(q+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  if(!read_factors(factors,q,facs_file))
    fatal_error("Error reading factor file. Exiting.\n");

  facs_file.close();

  phi_q=factors[q].phi;

  n_prims=num_prims(q,factors); // no of prim characters mod q
  if(!(prims=(unsigned int *) malloc(n_prims*sizeof(unsigned int))))
    fatal_error("Couldn't allocate memory for prims.\n");
  if(!(indices=(unsigned int *) malloc(n_prims*sizeof(unsigned int))))
    fatal_error("Couldn't allocate memory for indices.\n");
  if(!(f_omegas=(int_complex *) _aligned_malloc(n_prims*sizeof(int_complex),16)))
    fatal_error("Couldn't allocate memory for omegas.\n");
  if(!(neg_ones=(bool *) malloc(n_prims*sizeof(bool))))
    fatal_error("Couldn't allocate memory for neg_ones.\n");
  if(!(zs=(int_complex *) _aligned_malloc((n_prims*num_s*sizeof(int_complex),16))))
    fatal_error("Couldn't allocate memory for zs.\n");
  //printf("num_prims=%d\n",n_prims);

  if(!(omegas=(int_complex *) _aligned_malloc(phi_q*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  //printf("Initialising ws...\n");
  init_ws(_ws);

  //printf("Starting...\n");
  safe_odd(q,num_s,factors,q_n,offset_vec,a_n,omegas,t0,pibyq);
  safe_even();
  for(prim=0;prim<n_prims;prim++)
    {
      prim1=conj_j(prim,n_prims,q);
      real_p=(prim1==prim);
      if(prim1<prim)
	continue;
      //printf("Saving q:%d ind1:%d ind2:%d\n",q,prim,prim1);
      fwrite(&neg_ones[prim],sizeof(bool),1,out_file);
      fwrite(&real_p,sizeof(bool),1,out_file);
      fwrite(&indices[prim],sizeof(unsigned int),1,out_file);
      fwrite(&f_omegas[prim],sizeof(int_complex),1,out_file);
      //print_int_complex_str("saving",zs[prims[prim]*num_s+1]);
      fwrite(&zs[prim*num_s],sizeof(int_complex),num_s,out_file);
      if(!real_p) // it has a conjugate so save it
	{
	  //printf("Saving composite ind:%d\n",prim1);
	  //printf("prims[%d]=%d\n",prim1,prims[prim1]);

	  fwrite(&indices[prim1],sizeof(unsigned int),1,out_file);
	  fwrite(&f_omegas[prim1],sizeof(int_complex),1,out_file);
	  //print_int_complex_str("saving",zs[prims[prim1]*num_s+1]);
	  fwrite(&zs[prim1*num_s],sizeof(int_complex),num_s,out_file);
	}
      //else
      //printf("\n");
    }


  fclose(out_file);
  fclose(infile);

  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
     printf("safe successful completion on q=%d t0=%e N=%d\n",q,t0,num_s);
  return(0);
}
