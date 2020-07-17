/*
  File: f_even_odd.cpp

  Created: 10th Feb 2010

  Version: <v> = 1.1

  Last Modified: 8th March 2010

  1.0 Initial implementation
  1.1 Output now pre-sorted for f_even.cpp

  Dialect: C++


  Implementation notes: 

  Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

  Uses DFT to create F_hat_even(x) for x = [2*Pi*n0/B,2*Pi*n1/B) for all
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
  printf("Usage: f_hat_even%s (ifname) (facs_file) (ofname)\n",VERSION);
  printf("  (ifname)    - output from F_hat_odd_terms.c\n");
  printf("  (ofname)    - input to f.cpp (via spec file)\n");
  exit(1);
}


void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};

inline int phi(int p, int p_n)
{
  return(p_n-p_n/p);
}

long unsigned int pow_mod(long unsigned int a, long unsigned int b, long unsigned int m)
{
  long unsigned int a_pw=a,pw=b,res=1;
  while(true)
    {
      //printf("pw=%ld a_pw=%ld res=%ld\n",pw,a_pw,res);
      if(pw&1)
	res=(res*a_pw)%m;
      pw>>=1;
      if(pw==0)
	return(res);
      a_pw=(a_pw*a_pw)%m;
    }
}
      
unsigned long int pr(unsigned long int i, factor *factors)
{
  int phi=factors[i].phi;
  for(int p=2;p<i;p++)
    {
      if(gcd(p,i)!=1)
	continue;
      bool good=true;
      for(int j=0;j<factors[phi].num_facs;j++)
	{
	  if(pow_mod(p,phi/factors[phi].primes[j],i)!=1)
	    continue;
	  good=false;
	  break;
	}
      if(good)
	return(p);
    }
}

bool make_factors(factor *factors, int q_end)
{
  printf("Making factor database.\n");
  //printf("pow_mod(7,15,99)=%ld\n",pow_mod(7,19,99));
  int *primes,max_p=floor(sqrt(q_end));
  primes=(int *) malloc(sizeof(int)*(q_end+1));
  for(int i=2;i<=q_end;i++)
    primes[i]=0;
  for(int i=4;i<=q_end;i+=2)
    primes[i]=2;
  for(int i=3;i<=max_p;i+=2)
    if(primes[i]==0)
      for(int j=i*i;j<=q_end;j+=i)
	if(primes[j]==0)
	  primes[j]=i;
  printf("Prime sieve completed.\n");
  // now each entry primes[i] is 0 if i prime, else = smallest prime factor
  factors[3].num_facs=1;
  factors[3].primes[0]=3;
  factors[3].facs[0]=3;
  factors[4].num_facs=1;
  factors[4].primes[0]=2;
  factors[4].facs[0]=4;

  for(int f=5;f<=q_end;f++)
    if(primes[f]==0) // a prime
      {
	factors[f].num_facs=1;
	factors[f].primes[0]=f;
	factors[f].facs[0]=f;
      }
    else
      {
	factors[f].primes[0]=primes[f];
	int n=f/primes[f];
	if(factors[n].primes[0]==primes[f])
	  {
	    factors[f].num_facs=factors[n].num_facs;
	    factors[f].facs[0]=primes[f]*factors[n].facs[0];
	    for(int i=1;i<factors[n].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i];
		factors[f].facs[i]=factors[n].facs[i];
	      }
	  }
	else
	  {
	    factors[f].num_facs=factors[n].num_facs+1;
	    factors[f].facs[0]=primes[f];
	    factors[f].primes[0]=primes[f];
	    for(int i=1;i<factors[f].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i-1];
		factors[f].facs[i]=factors[n].facs[i-1];
	      }
	  }
      }
  free(primes);
  printf("Factors computed.\n");
  // now calculate phi(f)
  for(int i=3;i<=q_end;i++)
    {
      factors[i].phi=1;
      for(int j=0;j<factors[i].num_facs;j++)
	factors[i].phi*=phi(factors[i].primes[j],factors[i].facs[j]);
    }
  printf("phi computed.\n");

  //now do the prim roots
  factors[3].pr=2;
  factors[4].pr=3;
  //long unsigned int max_pr=3;
  for(int i=5;i<=q_end;i++)
    {
      if(((factors[i].num_facs==1)&&(factors[i].primes[0]!=2))|| // p^m, p an odd prime
	 ((factors[i].num_facs==2)&&(factors[i].facs[0]==2)))    // 2p^m
	{
	  factors[i].pr=pr(i,factors);
	  /*
	  if(factors[i].pr>max_pr)
	    {
	      max_pr=factors[i].pr;
	      printf("New largest pr=%lu for modulus %ld\n",max_pr,i);
	    }
	  */
	}
      else
	factors[i].pr=0;
    }
  printf("pr's computed.\n");
  /*  
  for(int i=3;i<50;i++)
    {
      printf("%ld has %ld factors, phi(%ld)=%ld, pr(%ld)=%ld, factors are",i,factors[i].num_facs,i,factors[i].phi,i,factors[i].pr);
      for(int j=0;j<factors[i].num_facs;j++)
	printf(" %ld %ld",factors[i].primes[j],factors[i].facs[j]);
      printf("\n");
    }
  */
  return(true);
}

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

void setup()
{
  d_log_pi=log(d_pi);
  d_log2=log(int_double(2));

  d_small=exp(int_double(LN_SMALL));
  d_small.left=-d_small.right;
  c_small=int_complex(d_small,d_small);
}

unsigned int *prims,*indices;
int_complex *f_omegas,*zs;
bool *neg_ones;

void make_l_even(unsigned int q, 
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
  //printf("Processing Q=%d\n",q);


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
		if(!neg_one)
		  {
		    if(s_done==0) // first time through, so finish off omega values
		      {
			f_omegas[prim_ptr]=finish_omega(omegas[i],neg_one);
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
      //printf("prepping omegas.\n");
      prep_omegas_nd(offset_vec,q,omegas,no_dims,dims,phi_q);
      //printf("omegas prepped.\n");
      while(s_done<num_s)
	{
	  x=two_pi_by_B*(n0+s_done);
	  fread(&M,sizeof(unsigned int),1,infile);
	  //printf("M=%lu\n",M);
	  err=read_int_double(infile);
	  for(i=0,j=1;j<min(q,M+1);j++)
	    if(co_prime(j,q))
	      {
		//printf("Placing value into Fft_vec[%u].\n",offset_vec[i]);
		Fft_vec[offset_vec[i++]]=read_int_complex(infile);
	      }
	  //printf("Non zeros entries in.\n");
	  for(;j<q;j++)
	    if(co_prime(j,q))
	      Fft_vec[offset_vec[i++]]=c_zero;
	  //printf("Doing nd FFT.\n");
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
		  if(!neg_one)
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
  double eta_even,B;
  bool real_p;

  no_clicks=clock(); // start timing
  setup();
  //printf("At program start.\n");

  //printf("argc=%d\n",argc);
  if(argc!=3)
    print_usage();

  
  infile=fopen(argv[1],"rb");
  if(!infile)
    fatal_error("Couldn't open infile. Exiting.\n");

  out_file=fopen(argv[2],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");

  fread(&q,sizeof(unsigned int),1,infile);
  QT=q*H;

  if(!(factors=(factor *) calloc(q+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  if(!make_factors(factors,q))
    fatal_error("Error reading factor file. Exiting.\n");


  fread(&eta_even,sizeof(double),1,infile);  
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
  fwrite(&eta_even,sizeof(double),1,out_file);
  fwrite(&num_s,sizeof(int),1,out_file);
  fwrite(&n0,sizeof(int),1,out_file);

  if(!(Fft_vec=(int_complex *) _aligned_malloc(phi_q*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(omegas=(int_complex *) _aligned_malloc(phi_q*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(phi_q,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  printf("Initialising ws...\n");
  init_ws(_ws);

  printf("Starting...\n");
  make_l_even(q,num_s,factors,q_n,offset_vec,a_n,omegas,B,n0,eta_even,infile);

  for(prim=0;prim<n_prims;prim++)
    {
      prim1=conj_j(prim,n_prims,q);
      real_p=(prim1==prim);
      if(prim1<prim)
	continue;
      //printf("Saving q:%d ind1:%d ind2:%d\n",q,prim,prim1);
      fwrite(&neg_ones[prim],sizeof(bool),1,out_file);
      fwrite(&real_p,sizeof(bool),1,out_file);
      if(neg_ones[prim])
	continue;
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
  printf("f_hat_even1.1 successful completion on q=%d n=[%d,%d)\n",q,n0,n1);
  return(0);
}
