/*
File: int-l-func-half1.0.cpp

Created: 5th November 2009

Version: <v> = 1.0
 
Dialect: C++

Requires: -lrt

Implementation notes: Takes Output from l-func-mpfi-1.1
            num_s (int)
            N (int)
            rn (int) must be 1 in this implementation
            NO_GAPS (int)

            im_s (double)
            gam_s (dcomplex)
            gam_s_2 (dcomplex)
            N*(NO_GAPS+1) h_rn (dcomplex)

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "1.0"


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <semaphore.h>
using namespace std;

#include "../includes/int_double11.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft.h"


void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: int-l-func%s (q-start) (q-end) (ifname) (ofname) (N) (semaphore name)\n",VERSION);
  printf("  (q-start)   - integer >=3\n");
  printf("  (q-end)     - integer >=(q-start)\n");
  printf("  (ifname)    - file with lattice values.\n");
  printf("  (facs_file) - file with factors data.\n");
  printf("  (ofname)    - output file.\n");
  printf("  (N)         - number of Taylor terms to use.\n");
  exit(1);
}


void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};


void prep_omegas_nd1(unsigned int *offsets,unsigned int q, unsigned int n_dims,int *dims,
					unsigned int phi_q)
{
	unsigned int i,w_ptr=0;

	for(i=0;i<n_dims;i++)
		{
			ws_ptr[i]=w_ptr;
			if((dims[i]==2)||(dims[i]==4))
				continue;

			if(bin_power_2(dims[i]))
			{
				w_ptr+=dims[i];
				continue;
			}

			if(dims[i]<=MAX_SIMPLE_DFT)
				continue;

			init_bluestein_fft(dims[i],&conv_sizes[i],&bs[w_ptr],&b_star_conjs[w_ptr]);
			w_ptr+=conv_sizes[i];
		}
}

#define HUR_HALF_TERMS (10)
#define HUR_HALF_K (10)
#define HUR_HALF_ERR ((double) (2*HUR_HALF_K-0.5)*(2*HUR_HALF_K-0.5))

int_double s_array[HUR_HALF_K];

inline void init_hur_half()
{
}

inline int_complex hur_half(const int_double &alpha)
{
	int i;
  int_double err,n_alpha;
  int_double res=d_zero,n_alpha_s,term;
  double s1=0.5;
  int_double s_array[MAX_BERNOULLI_K];

  n_alpha=alpha+HUR_HALF_TERMS;
  n_alpha_s=pow(n_alpha,0.5);

  for(i=0;i<HUR_HALF_TERMS;i++)
    {
      res=res+pow(alpha+i,-0.5);
      //      if((i%1000)==0)
      //	printf("relative error in res is %d %d\n",rel_error(res.real),rel_error(res.imag));
    }

  //debug;
  res=res-n_alpha_s*2;

  n_alpha_s=n_alpha_s/n_alpha; // n_alpha_s=(n+alpha)^(-s)
		
  res=res+n_alpha_s/2;

  s_array[0]=n_alpha_s/(n_alpha*2);

  for(i=1;i<HUR_HALF_K;i++)
    {
      s1=s1+1;
      s_array[i]=s_array[i-1]*s1/n_alpha;
      s1=s1+1;
      s_array[i]=s_array[i]*s1/n_alpha;
    }
  
  for(i=0;i<HUR_HALF_K;i++)
    {
//		print_int_double_str("",s_array[i]);
//		print_int_double_str("",h_bernoulli[i]);
      term=s_array[i]*h_bernoulli[i];
      res=res+term;
    }

  err=term/HUR_HALF_ERR;
  if(err.left<=err.right)
    err.right=err.left;
  else
    err.left=err.right;
//  print_int_double_str("last term=",term);
//  print_int_double_str("error term is ",err);
  res+=err;
//  print_int_double_str("alpha=",alpha);
//  print_int_double_str("Zeta(0.5,alpha)=",res);
  return(res);

	//return(hurwitz(c_half,alpha));
}

void make_l(unsigned int q,
			factor *factors,
			unsigned int *q_n,
			unsigned int *a_n,
			unsigned int *offset_vec,
			int_complex *omegas)
{
	unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS];
	int_complex omega,omega_a,z,z1,q_comb;
	int dims[MAX_FACS];
	int no_dims;
	bool power_2,primitive,neg_one;
	unsigned int pr=factors[q].pr;
	unsigned int i,j,offset,s_done;
	long double lnq=ln_q(q);
    if((q%100)==0)
	   printf("Processing Q=%d\n",q);

	//fwrite(&q,sizeof(unsigned int),1,out_file);
	//n_prims=num_prims(q,factors); // no of prim characters mod q
	//fwrite(&n_prims,sizeof(unsigned int),1,out_file);

	//q_minus_half=pow(int_double(q),-0.5);
	j=0;
	for(i=1;i<q;i++)
		if(co_prime(i,q))
			a_n[j++]=i;
	if(pr)  // q has a primitive root so nice and easy 1-d FFT
	{
		//printf("has pr\n");
		fill_q_ns(q_n,pr,phi_q,q);
		init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
		//prep_omegas(omegas,q,phi_q,q_n,a,bs,b_star_conjs,conv_sizes[0],b_spare);
		//for(i=0;i<phi_q;i++)
		//  print_int_complex_str("prepped omegas[]= ",omegas[i]);
		for(i=1;i<q;i++)
			if(co_prime(i,q))
			{
				Fft_vec[q_n[i]]=hur_half(int_double(i)/q);
				//print_int_complex_str("Hurwitz returned ",Fft_vec[q_n[i]]);
			}
		bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare);

		// now fft_vec contains L(chi,s)*q^(s)
		for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
			if(prim_p(i,factors[q].primes[0]))        
			{
				if(contains_zero(Fft_vec[i].real)&&(contains_zero(Fft_vec[i].imag)))
				{
					printf("Might have a zero at q=%d index=%d.\n",q,i);
					print_int_complex_str("L(1/2,chi)=",Fft_vec[i]);
				}
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
		//debug;
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
		prep_omegas_nd1(offset_vec,q,no_dims,dims,phi_q);

		for(i=0,j=1;i<phi_q;i++)
		{
			while(!co_prime(j,q)) j++;
			//printf("Zeta(1/2,%d/%d)=",j,q);print_int_complex(real_hur(0.5,int_double(j)/q));
			Fft_vec[offset_vec[i]]=hur_half(int_double(j++)/q);
			Fft_vec[offset_vec[i]].imag=d_zero;
		}
		//debug;
		do_nd_fft(Fft_vec,no_dims,dims,phi_q);

		for(i=0;i<factors[q].num_facs;i++)
			coords[i]=0;
		//debug;
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
					if(contains_zero(Fft_vec[i].real)&&(contains_zero(Fft_vec[i].imag)))
					{
						printf("Might have a zero at q=%d index=%d.\n",q,i);
					print_int_complex_str("L(1/2,chi)=",Fft_vec[i]);
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
	}
}

int main(int argc, char **argv)
{

  _fpu_rndd();
  set_h_bernoulli();

  int_complex s,gam,*s_array,*r_n_vals,*omegas;
  im_s *im_s_vec;
  unsigned int *a_n,*q_n,*offset_vec,num_fracs,q,max_num_prims,n_prim,i;
  int no_gaps,q_start,q_end,num_s,N,rn,file_N;
  int_double mgap,im_s_2_mod;
  double gap;
  FILE *in_file,*out_file;
  ifstream facs_file;
  factor *factors;
  clock_t no_clicks,no_clicks1;

  no_clicks=clock(); // start timing
  //printf("At program start.\n");

  if(argc!=4)
    print_usage();
  q_start=atoi(argv[1]);
  if(q_start<3)
    print_usage();
  q_end=atoi(argv[2]);
  if(q_end<q_start)
    print_usage();
  facs_file.open(argv[3]);
    if(!facs_file.is_open())
      fatal_error("Couldnt open factors file. Exiting.\n");

  if(!(factors=(factor *) calloc(q_end+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  if(!read_factors(factors,q_end,facs_file))
    fatal_error("Error reading factor file. Exiting.\n");

  facs_file.close();
  if(!(omegas=(int_complex *) _aligned_malloc((q_end-1)*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  init_ws(_ws);
  for(q=q_start;q<=q_end;q++)
  {
    if((q&3)!=2)
      make_l(q,factors,q_n,offset_vec,a_n,omegas);
  }

  fclose(out_file);

  //printf("FFT Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks1)/((double) CLOCKS_PER_SEC));
  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  //printf("Lambda value straddled zero %d times\n",temp_count);
  printf("l-func successful completion on q=%d-%d file %s\n",q_start,q_end,argv[3]);
  return(0);
}
