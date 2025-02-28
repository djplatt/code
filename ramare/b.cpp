/*
File: b.cpp

Created: 16th November

Version: <v> = Original

Dialect: C++

Requires: -lrt

Build instructions: g++ -fomit-frame-pointer -frounding-math -finline-functions -lcrlibm

By: DJ Platt
    Heilbronn Institute
    Bristol University

Copyright 2012,2013,2014.

Compute L'(0,chi)/L(0,chi) for all odd primitive characters 
Compute lim s->0 L'(0,chi)/L(0,chi)-1/s for all even primitivecharacters
of modulus q-start to q-end inclusive.
Sum over all characters for each q.
Divide by phi_2(q) (number of primitive characters)

The author is a Heilbronn Research Fellow */

#define VERSION ""

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <semaphore.h>
#include <assert.h>
//#include <algorithm>
using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/hurwitz1.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft-half.h"
//#include "L1.h"

void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: b%s (q-start) (q-end) (ofname) (num_facs)\n",VERSION);
  printf("  (q-start)   - integer >=3\n");
  printf("  (q-end)     - integer >=(q-start)\n");
  printf("  (ofname)    - output file.\n");
  printf("  (num_facs)  - do q with this # of prime factors (or all q if num_facs == 0).\n");
  exit(1);
}


void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};

// returns the index of the character conjugate to j
// return j if it is real
// num_chi is number of primitive characters
inline int conj_j ( int j,  int num_chi,  int q)
{
	if(q&7) // q not divisible by 8
		return(num_chi-j-1);
	if(j<(num_chi>>1))
		return((num_chi>>1)-j-1);
	return(num_chi+(num_chi>>1)-j-1);
}


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

void init_fft_nd(unsigned int *offsets,unsigned int q,unsigned int n_dims,int *dims,
					unsigned int phi_q)
{
  unsigned int n,n1,i,w_ptr=0;
  int_double theta=int_double(2.0)/q,q_minus_half;
  /*
  q_minus_half=d_one/sqrt(int_double(q));//exp(log(int_double(q))*(-0.5));
  n1=0;
  for(n=1;n<q;n++)
    if(co_prime(n,q))
      {
	sin_cospi(theta*n,&omegas[offsets[n1]].imag,&omegas[offsets[n1]].real);
	n1++;
      }
  */
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
      //printf("initialising bluestein dim %ld\n",i);
      init_bluestein_fft(dims[i],&conv_sizes[i],&bs[w_ptr],&b_star_conjs[w_ptr]);
      //printf("bluestein initialised.\n");
      w_ptr+=conv_sizes[i];
    }
  //printf("doing %ld dimension fft\n",n_dims);
  //do_nd_fft(omegas,n_dims,dims,phi_q);
  //printf("%ld dimension fft finished\n",n_dims);  

  //for(n=0;n<phi_q;n++)
  //omegas[n]*=q_minus_half;
}

// compute b=L(0,chi)/L'(0,chi) for odd characters
void b(unsigned int q, 
      factor *factors,
      unsigned int *q_n,
      unsigned int *a_n,
      unsigned int *offset_vec,
      FILE *out_file)
{
  unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS];
  int_complex z,z1,q_comb;
  //int_double q_minus_half=sqrt(d_one/q);
  int dims[MAX_FACS];
  int no_dims;
  bool power_2,primitive,neg_one;
  unsigned int pr=factors[q].pr,n_prims;
  unsigned int i,j,k,offset,s_done;
  int_double log_q=log(int_double(q)),x;

  //printf("Processing Q=%d\n",q);

  //fwrite(&q,sizeof(unsigned int),1,out_file);
  n_prims=num_prims(q,factors); // no of prim characters mod q
  //fwrite(&n_prims,sizeof(unsigned int),1,out_file);

  j=0;
  for(i=1;i<q;i++)
    if(co_prime(i,q))
      a_n[j++]=i;
  if(pr)  // q has a primitive root so nice and easy 1-d FFT
  {
    fill_q_ns(q_n,pr,phi_q,q);
    init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
    for(i=1;i<q;i++)
      if(co_prime(i,q))
	Fft_vec[q_n[i]]=-int_double(i)/q;// +0.5
      
      

      bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare);

      // now Fft_vec contains sum chi(a) a/q
      for(i=0;i<phi_q;i++)
	Fft_copy[i]=Fft_vec[i];
      // now copied to Fft_copy
      //printf("Setting up lngammas.\n");
      for(i=1;i<q;i++)
	if(co_prime(i,q))
	  Fft_vec[q_n[i]]=lngamma(int_double(i)/q);//-0.5*d_ln_two_pi;
      /*
      for(i=0;i<phi_q;i++)
	{
	  printf("Zeta(1/2,a/%ld)=",q);
	  print_int_complex_str("",Fft_vec[i]);
	}
      */
      bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare);
      for(i=0;i<phi_q;i++)
	Fft_copy1[i]=Fft_vec[i]; // contains sum chi(a) log gamma (a/q)

      for(i=1;i<q;i++)
	if(co_prime(i,q))
	  Fft_vec[q_n[i]]=hurwitz0(int_double(i)/q);
      bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare);

      //Fft_vec contains sum chi(a) zeta''(0,a/q)
      int_complex res;
      int_double all_sum=0.0;
      for(i=1,k=0;i<phi_q;i++)  // i=0 corresponds to principal chi
	if(prim_p(i,factors[q].primes[0]))//&&(conj_j(k,n_prims,q)>=k)) // only do primitive conjugates
	  {
	    if(i&1) // odd
	      {
		res=Fft_copy1[i]/Fft_copy[i]-d_ln_two_pi-d_gamma;
		//odd_sum+=norm(res);
	      }
	    else // even
	      {
		Fft_copy1[i].real*=2.0;
		Fft_copy1[i].imag*=2.0;
		res=Fft_vec[i]/Fft_copy1[i]-d_ln_two_pi-d_gamma;
		//even_sum+=norm(res);
	      }
	    all_sum+=norm(res);
	    k++;
	  }
      printf("q= %d phi_2= %d %20.17e %20.17e ",q,k,all_sum.left,-all_sum.right);
      all_sum/=k;
      printf("%20.17e %20.17e\n",all_sum.left,-all_sum.right);

  }
  
  else // q doesn't have a primitive root
    {
      //printf("%ld does not have a generator.\n",q);
      no_dims=factors[q].num_facs;
      fac=factors[q].facs[0];        // the first p^n
      power_2=(factors[fac].pr==0);  // no conductor => p=2, n>=3
      if(power_2)
	{
	  //printf("%ld contains a power of 2^n, n>=3\n");
	  no_dims++;
	  fill_q_ns_2s(q_n,fac);    // use the {-1,1}X{5} trick
	  for(i=1;i<factors[q].num_facs;i++) // move all the dimensions up one
	    dims[i+1]=factors[factors[q].facs[i]].phi;
	  dims[1]=factors[q].facs[0]/4;
	  dims[0]=2;                         // slip in a two
	}
      else
	{
	  //printf("Filling q_ns for factor %ld.\n",fac);
	  fill_q_ns(q_n,factors[fac].pr,factors[fac].phi,fac); // use the generator
	  //printf("q_ns filled.\n");
	  for(i=0;i<factors[q].num_facs;i++)
	    dims[i]=factors[factors[q].facs[i]].phi;
	}
      offset=fac;
      for(i=1;i<factors[q].num_facs;i++)  // do the rest on the factors, all will have generators
	{
	  fac=factors[q].facs[i];      
	  pr=factors[fac].pr;
	  //printf("Filling q_ns for factor %ld.\n",fac);
	  fill_q_ns(&q_n[offset],pr,factors[fac].phi,fac);  // use the generator
	  //printf("q_ns filled.\n");
	  offset+=fac;
	}
      s_done=0;
      //printf("making offsets.\n");
      make_offsets(q,q_n,offset_vec,factors,a_n);    // reverse q_n so we know how to populate fft vector
      //printf("prepping omegas.\n");
      init_fft_nd(offset_vec,q,no_dims,dims,phi_q);
      for(j=1,i=0;j<q;j++)
	if(co_prime(j,q))
	  Fft_vec[offset_vec[i++]]=-int_double(j)/q; // +0.5
      do_nd_fft(Fft_vec,no_dims,dims,phi_q);
      for(i=0;i<phi_q;i++)
	Fft_copy[i]=Fft_vec[i];
      
      for(j=1,i=0;j<q;j++)
	if(co_prime(j,q))
	  Fft_vec[offset_vec[i++]]=lngamma(int_double(j)/q);//-0.5*d_ln_two_pi;
      do_nd_fft(Fft_vec,no_dims,dims,phi_q);
      for(i=0;i<phi_q;i++)
	Fft_copy1[i]=Fft_vec[i];
      
      for(j=1,i=0;j<q;j++)
	if(co_prime(j,q))
	  Fft_vec[offset_vec[i++]]=hurwitz0(int_double(j)/q);
      do_nd_fft(Fft_vec,no_dims,dims,phi_q);
      
      
      for(i=0;i<factors[q].num_facs;i++)
	coords[i]=0;
      
      int_complex res;
      int_double all_sum=0.0;
      
      for(i=0,k=0;i<phi_q;i++)
	{
	  primitive=true;
	  for(j=0;j<factors[q].num_facs;j++)
	    if(coords[j]%factors[q].primes[j]==0)
	      {
		primitive=false;
		break;
	      }
	  
	  if(primitive)//&&(conj_j(k,n_prims,q)>=k))
	    {
	      neg_one=neg_one_p(coords,factors[q].num_facs);
	      if(power_2&&(i<(phi_q>>1)))
		neg_one=!neg_one;
	      if(neg_one)
		{
		  //printf("character %d is odd. ",i);
		  res=Fft_copy1[i]/Fft_copy[i]-d_ln_two_pi-d_gamma;
		  //print_int_complex_str("",res);
		  //odd_sum+=norm(res);
		}
	      else // even
		{
		  //printf("character %d is even. ",i);
		  Fft_copy1[i].real*=2.0;
		  Fft_copy1[i].imag*=2.0;
		  res=Fft_vec[i]/Fft_copy1[i]-d_ln_two_pi-d_gamma;
		  //print_int_complex_str("",res);
		  //even_sum+=norm(res);
		}
	      all_sum+=norm(res);
	      k++; 
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
      printf("q= %d phi_2= %d %20.17e %20.17e ",q,k,all_sum.left,-all_sum.right);
      all_sum/=k;
      printf("%20.17e %20.17e\n",all_sum.left,-all_sum.right);
    }
}

int main(int argc, char **argv)
{

  _fpu_rndd();

  unsigned int *a_n,*q_n,*offset_vec,q,i;
  int q_start,q_end;
  FILE *out_file;
  factor *factors;

  if(argc!=5)
    print_usage();
  q_start=atoi(argv[1]);
  if(q_start<3)
    q_start=3;
  q_end=atoi(argv[2]);
  if(q_end<q_start)
    q_start=q_end;
  
  int num_facs=atoi(argv[4]);

  if(q_end>MAX_Q)
    {
      printf("q_end (%ld) exceeds MAX_Q (%ld). Exiting.\n",q_end,MAX_Q);
      exit(0);
    }
  

  out_file=fopen(argv[3],"wb");
  if(!out_file)
    {
      perror("Error opening out_file. ");
      printf("Failed to open file |%s|\n",argv[3]);
      fatal_error("Couldn't open in_file. Exiting.\n");
    }
  //fwrite(&q_start,sizeof(int),1,out_file);
  //fwrite(&q_end,sizeof(int),1,out_file);

  if(!(factors=(factor *) calloc(q_end+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  //printf("Making factors.\n");
  if(!make_factors(factors,q_end))
    fatal_error("Error creating factors. Exiting.\n");

  printf("allocated memory\n");
  
  //sem_wait(lfuncsem);
  
  //printf("remember to uncomment sem_wait etc.\n");

  if(!(q_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  printf("initialising ws\n");
  init_ws();
  printf("ws initialised\n");
  for(q=q_start;q<=q_end;q++)
    if((q&3)!=2)
      if((num_facs==0)||(factors[q].num_facs==num_facs))
	b(q,factors,q_n,offset_vec,a_n,out_file);

  //fclose(out_file);

  //printf("FFT Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks1)/((double) CLOCKS_PER_SEC));
  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  //printf("Lambda value straddled zero %d times\n",temp_count);
  printf("l-func successful completion on q=%d-%d file %s\n",q_start,q_end,argv[3]);


  return(0);
}
