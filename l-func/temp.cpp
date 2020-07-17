/*
File: int-l-func5.10.cpp

Created: 10th Jan 2009

Version: <v> = 5.10

Last Modified: 18th June 2012

5.0 changed summation code in calc_rn and calc_rn_neg
5.1 now multiply by (q/pi)^(it/2)
5.2 moved (q/pi)^(it/2) and q^(-s) into FFT
5.3 changed calc of q_comb for accuracy, introduced im_s_vec[].pi_minus_it_2
5.4 using semaphores to access lattice files
5.5 read lattice file in one go using read
5.6 improved calculation of hurwitz values
5.7 now uses crlibm 
5.8 makes its own factors, no longer uses semaphores
5.9 now includes int_fft1.2.h
    (this computes DFT of b once only for Bluestein convolution)
    thus no longer need b_spare
5.10 Exploit fact that phi_q is always even to halve FFT lengths.
     now includes int_fft1.3.h

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

Build instructions: g++ -O1 -mfpmath=387 -fomit-frame-pointer -frounding-math -finline-functions

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "5.10"


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <semaphore.h>
//#include <assert.h>
using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft1.3.h"
/*
extern "C" {
#include "../includes/nit.h"
}
*/
#include "../includes/upsamdefs.h"
//#include "../includes/qit_struct.h"
//#include "../includes/qit.h"

void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: int-l-func%s (q) (ifname) (ofname) (N)\n",VERSION);
  printf("  (q)         - integer >=3\n");
  printf("  (ifname)    - file with lattice values.\n");
  printf("  (ofname)    - output file.\n");
  printf("  (N)         - number of Taylor terms to use.\n");
  //  printf("  (qit file)  - file with qit values.\n");
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

// only a vector here. last element used for taylor error
void create_s_array(const int_complex &s,int_complex *s_array, unsigned int n, const double gap)
{
  unsigned int i;
  s_array[0].real=int_double(gap/2); // 0.5*gap is representable
  s_array[0].imag=s.imag*gap;
  for(i=1;i<n;i++)
    {
      s_array[i]=s_array[i-1]*(s+i)*gap;
      s_array[i]/=i+1;
    }
}

void print_s_array(int_complex *s_array, unsigned int n)
{
  for(int i=0;i<n;i++)
    {
      printf("s_array[%d]=",i);
      print_int_complex_str("",s_array[i]);
    }
}

void print_r_n_vec(int_complex *r_n_vec,unsigned int n)
{
  unsigned int i;

  for(i=0;i<n;i++)
    {
      cout << "r_n_vec[" << i << "] is ";
      print_int_complex(r_n_vec[i]);
      cout << endl;
    };
}

inline int_complex calc_r_n(unsigned int r_n_terms, int_complex *s_array, int_double &delta,
      int_complex *r_n_vec)
{
	unsigned int i=r_n_terms-2;
	int_complex res=r_n_vec[i+1]*s_array[i];
	for(;i>0;i--)
	  res=res*delta+r_n_vec[i]*s_array[i-1];
	res=res*delta+r_n_vec[0];
	return(res);
}

#define PRINT_ERROR_STEP 16
void print_errors(int_complex &s,int_complex *r_n_vals, unsigned int N, unsigned int no_gaps,unsigned int M)
{
  unsigned int i,j;
  int err;
  printf("Rel Error in r_n_vals at Im(s)=");
  print_int_double(s.imag);
  printf("\n");
  for(i=0;i<=PRINT_ERROR_STEP;i++)
  {
    printf("At%2d/%2d ",i,PRINT_ERROR_STEP);
    for(j=0;j<N;j++)
    {
      if(j==M)
        printf("|");
      err=rel_error(r_n_vals[N*i*no_gaps/PRINT_ERROR_STEP+j].real);
      printf("%3d",err);
    }
    printf("\n");
  }
}

double *lattice_buff;
int lattice_buff_ptr=0;

inline double read_double()
{
  return(lattice_buff[lattice_buff_ptr++]);
}

inline int_double read_int_double()
{
  int_double ans;
  ans.left=lattice_buff[lattice_buff_ptr++];
  ans.right=-lattice_buff[lattice_buff_ptr++];
  return(ans);
}

inline int_complex read_int_complex()
{
  int_complex ans;
  ans.real.left=lattice_buff[lattice_buff_ptr++];
  ans.real.right=-lattice_buff[lattice_buff_ptr++];
  ans.imag.left=lattice_buff[lattice_buff_ptr++];
  ans.imag.right=-lattice_buff[lattice_buff_ptr++];
  return(ans);
}

inline void skip_int_complex(unsigned int n)
{
  lattice_buff_ptr+=n*sizeof(int_complex);
}

int_complex taylor_error(int_complex &s, unsigned int N,double gap,int_complex *s_array)
{
  int_double a,r,err;

  // set a to max mod of H_RN(s+N,alpha) with alpha in [0,1]
    a=int_double((double)N+0.5);
    a=a/((double)N-0.5);
    // now multiply by s(s+1)..(s+N)*(gap/2)^N/N! assuming a/q is exactly half way
    a=a*sqrt(norm(s_array[N-1]))/(1<<N);
    r=sqrt(norm(s+N))*gap;
    r=r/(N+N+2);
  err=a/(r-1.0); // this is -ve what we want but never mind
  if(err.right>=0.0)  // interval in (-infty,0.0]
    err.right=err.left;
  else
  {
    if(err.left<=0.0) // interval in [0.0,infty)
      err.left=err.right;
    else
    {
      if(err.right>=err.left)
        err.left=-err.right;
      else
        err.right=-err.left;
    }
  }
  return(int_complex(err,err));
}

unsigned int hur_ptr=0,num_fracs;
int_complex *hur_vec;
/*
inline int_complex calc_qit(unsigned int n_dt, qit_t *qit)
{
  int_complex res=c_one;
  //printf("n_dt=%d\n",n_dt);
  for(int bit=0;bit<NUM_QIT_BITS;bit++,n_dt>>=1)
    if(n_dt&1)
      {
	res*=qit[0].bits[bit];
	//printf("bit %d is set\n",bit);
	//print_int_complex_str("res is now",res);
      }
  return(res);
}

inline int_complex q_it(unsigned int q, double t, qit_t *qit)
{
  double dn=t/(one_over_two_B/2.0);
  int n=dn;
  if((dn-(double) n)!=0.0)
    {
      printf("q_it called with non-integral multiple of 2B. Exiting.\n");
      exit(1);
    }
  return(calc_qit(n,qit));
}
*/

inline int_complex q_it1(int_double lnq, double t)
{
  int_complex res;
  int_double theta=lnq*t;
  sin_cos(theta,&res.imag,&res.real);
  //print_int_double_str("qit of ",lnq);
  //printf("with t=%f\n",t);
  //print_int_complex_str("returning ",res);
  return(res);
}

void make_l(unsigned int q, 
      unsigned int num_s, 
      factor *factors,
      unsigned int *q_n,
      im_s *im_s_vec,
      unsigned int *a_n,
      unsigned int *offset_vec,
      int_complex *omegas,
	    FILE *out_file)
//	    qit_t *qits)
{
  unsigned int phi_q=factors[q].phi,pq2=phi_q/2,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS];
  unsigned int my_hur_ptr,my_hur_ptr1;
  int_complex omega,omega_a,z,z1,q_comb;
  int_double q_minus_half,lnq=log(int_double(q));
  int dims[MAX_FACS];
  int no_dims;
  bool power_2,primitive,neg_one;
  unsigned int pr=factors[q].pr,n_prims;
  unsigned int i,j,offset,s_done;

  //printf("Processing Q=%d\n",q);

  fwrite(&q,sizeof(unsigned int),1,out_file);
  n_prims=num_prims(q,factors); // no of prim characters mod q
  fwrite(&n_prims,sizeof(unsigned int),1,out_file);

  q_minus_half=d_one/sqrt(int_double(q)); // pow(int_double(q),-0.5);
  j=0;
  for(i=1;i<q;i++)
    if(co_prime(i,q))
      a_n[j++]=i;
  if(pr)  // q has a primitive root so nice and easy 1-d FFT
  {
    fill_q_ns(q_n,pr,phi_q,q);
    init_bluestein_fft(pq2,conv_sizes,bs,b_star_conjs);
    int_complex *twiddles=(int_complex *)_aligned_malloc(sizeof(int_complex)*pq2,16);
    if(!twiddles)
      fatal_error("Failed to allocate memory for twiddles. Exiting.");
    // this is messy because I mixed +/- transform between simple and bluestein. Doh!
    twiddles[0]=c_one;
    if(pq2<=MAX_SIMPLE_DFT)
      for(i=1;i<pq2;i++)
	sin_cospi(int_double(i)/pq2,&twiddles[i].imag,&twiddles[i].real);
    else
      for(i=1;i<pq2;i++)
	sin_cospi(-int_double(i)/pq2,&twiddles[i].imag,&twiddles[i].real);
    //for(i=0;i<pq2;i++)
    //print_int_complex_str("twiddles[]=",twiddles[i]);

    prep_omegas_new(omegas,q,phi_q,q_n,a,bs,b_star_conjs,conv_sizes[0],twiddles);
    //printf("omegas prepped.\n");
    //init_bluestein_fft(pq2,conv_sizes,bs,b_star_conjs);
    s_done=0;
    my_hur_ptr1=hur_ptr;
    my_hur_ptr=hur_ptr;
//    printf("my_hur_ptr %ld\n",my_hur_ptr);
    while(s_done<num_s)
    {
      q_comb=im_s_vec[s_done].pi_minus_it_2*q_minus_half*conj(q_it1(lnq,im_s_vec[s_done].im_s/2.0));/*n_it(q,im_s_vec[s_done].im_s/2.0));*/
      //q_it(q,im_s_vec[s_done].im_s/2.0,qits+q));
      //printf("copying hurwitz values.\n");
      for(i=1;i<q;i++)
        if(co_prime(i,q))
	  {
	    //printf("Setting %d'th element of Fft_vec to ",q_n[i]);
	    //print_int_complex_str("",hur_vec[my_hur_ptr]);
	    Fft_vec[q_n[i]]=hur_vec[my_hur_ptr++];
	  }
      /*
      printf("Pre FFT.\n");
      int r_err=-1000,a_err=-1000,sum_r_err=0,sum_a_err=0;
      for(long unsigned int j=0;j<phi_q;j++)
	{
	  int this_r_err=rel_error(Fft_vec[j].real);
	  int this_a_err=abs_error(Fft_vec[j].real);
	  sum_r_err+=this_r_err;
	  sum_a_err+=this_a_err;
	  if(this_r_err>r_err) r_err=this_r_err;
	  if(this_a_err>a_err) a_err=this_a_err;
	}
       printf("worst case relative error 10^{%d} absolute error 10^{%d}\n",r_err,a_err);
       printf("average relative error 10^{%g} absolute error 10^{%g}\n",sum_r_err/(double) phi_q,sum_a_err/(double) phi_q);
      */
      //for(i=0;i<phi_q;i++)
      //print_int_complex_str("Fft_vec[]=",Fft_vec[i]);
      /*
      for(i=0;i<20;i++)
	{
	  printf("abs error before FFT[i]=%d\n",i,abs_error(Fft_vec[i].real));
	  print_int_complex_str("FFT[i]=",Fft_vec[i]);
	}
      exit(0);
      */
      bluestein_fft(2,Fft_vec,pq2,a,bs,b_star_conjs,conv_sizes[0]); // evens
      bluestein_fft(2,Fft_vec+1,pq2,a,bs,b_star_conjs,conv_sizes[0]); // odds
      /*
      printf("Post FFT.\n");
      r_err=-1000;a_err=-1000;sum_r_err=0;sum_a_err=0;
      for(long unsigned int j=0;j<phi_q;j++)
	{
	  int this_r_err=rel_error(Fft_vec[j].real);
	  int this_a_err=abs_error(Fft_vec[j].real);
	  sum_r_err+=this_r_err;
	  sum_a_err+=this_a_err;
	  if(this_r_err>r_err) r_err=this_r_err;
	  if(this_a_err>a_err) a_err=this_a_err;
	}
       printf("worst case relative error 10^{%d} absolute error 10^{%d}\n",r_err,a_err);
       printf("average relative error 10^{%g} absolute error 10^{%g}\n",sum_r_err/(double) phi_q,sum_a_err/(double) phi_q);
      */
      for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
        if(prim_p(i,factors[q].primes[0]))        
        {
          neg_one=(i&1);
	  int omega_index;
	  if(i<pq2)
	    omega_index=i+i;
	  else
	    omega_index=(i-pq2)*2+1;

          if(s_done==0) // first time through, so finish off omega values
          {
            omegas[omega_index]=finish_omega(omegas[omega_index],neg_one);
	    
	    //if(phi_q>MAX_SIMPLE_DFT) omegas[i]=conj(omegas[i]);

	    //print_int_complex_str("finished omega[]= ",omegas[i]);
            fwrite(&i,sizeof(unsigned int),1,out_file);
            fwrite(&omegas[omega_index],sizeof(int_complex),1,out_file);
            fwrite(&neg_one,sizeof(bool),1,out_file);
          }
	  if(i<pq2)
	    z=Fft_vec[i+i]+twiddles[i]*Fft_vec[i+i+1];
	  else
	    z=Fft_vec[i+i-phi_q]-twiddles[i-pq2]*Fft_vec[i+i-phi_q+1];
	  //print_int_complex_str("z=",z);
	  z*=q_comb;
	  z*=omegas[omega_index];
          if(neg_one) // chi(-1)=-1
            z*=im_s_vec[s_done].lambda_s_a;
          else
            z*=im_s_vec[s_done].lambda_s;
	  //printf("abs error in z=%d\n",abs_error(z.real));
          if(!contains_zero(z.imag))
          {
            printf("Imaginary part of z does not contain zero.\n");
	    print_int_complex_str("Fft_vec[]= ",Fft_vec[i]);
	    print_int_complex_str("z= ",z);
	    printf("im_s=%20.18e\n",im_s_vec[s_done].im_s);
	    printf("q=%d\n",q);
	    print_int_complex_str("omega ",omegas[omega_index]);
	    print_int_complex_str("q_comb ",q_comb);
	    print_int_complex_str("hur_vec[-2]=",hur_vec[my_hur_ptr-2]);
	    print_int_complex_str("hur_vec[-2]=",hur_vec[my_hur_ptr-1]);
	    fatal_error("Exiting.");
          }
	  /*
	  if((z.real.left<0)&&(z.real.right<0))
	    {
	      print_int_complex_str("fft         =",Fft_vec[i]);
	      print_int_complex_str("omega       =",omegas[i]);
	      print_int_complex_str("q^(-it/2)/2 =",q_minus_it_2);
	      print_int_complex_str("pi^(-it/2)  =",im_s_vec[s_done].pi_minus_it_2);
	      print_int_double_str("q^(-1/2)     =",q_minus_half);
	      if(neg_one)
		print_int_complex_str("Gamma        =",im_s_vec[s_done].lambda_s_a)
	      else
		print_int_complex_str("Gamma        =",im_s_vec[s_done].lambda_s)
	      print_int_complex_str("z            =",z);
	     
	      temp_count++;
	    }
	  */
	  //printf("absolute error in Re(z)=%d\n",abs_error(z.real));
          fwrite(&z.real,sizeof(int_double),1,out_file);
        }
        s_done++;
	my_hur_ptr1+=num_fracs;
	my_hur_ptr=my_hur_ptr1;
    }
    free(twiddles);
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
    my_hur_ptr1=hur_ptr;
    my_hur_ptr=hur_ptr;
//    printf("my_hur_ptr %ld\n",my_hur_ptr);
    while(s_done<num_s)
    {
      //sin_cos(-log(int_double(q))*(double) (im_s_vec[s_done].im_s/2.0),&q_minus_it_2.imag,&q_minus_it_2.real);
      //q_comb=q_minus_it(q,im_s_vec[s_done].im_s/2)*im_s_vec[s_done].pi_minus_it_2*pow(int_double(q),-0.5);
      //print_int_complex_str("q^it=",q_it(q,im_s_vec[s_done].im_s/2.0));
      q_comb=im_s_vec[s_done].pi_minus_it_2*q_minus_half*conj(q_it1(lnq,im_s_vec[s_done].im_s/2.0));//q_it(q,im_s_vec[s_done].im_s/2.0,qits+q));
      for(i=0;i<phi_q;i++)
	Fft_vec[offset_vec[i]]=hur_vec[my_hur_ptr++];

      do_nd_fft(Fft_vec,no_dims,dims,phi_q);
      //printf("FFT Done\n");
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
	      z=Fft_vec[i];
	      neg_one=neg_one_p(coords,factors[q].num_facs);
	      if(power_2&&(i<(phi_q>>1)))
		neg_one=!neg_one;
	      
	      if(s_done==0) // first time for this im_s and chi
		{
		  //printf("doing omegas\n");
		  omegas[i]=finish_omega(omegas[i],neg_one);
		  //printf("omegas done\n");
		  fwrite(&i,sizeof(unsigned int),1,out_file);
		  fwrite(&omegas[i],sizeof(int_complex),1,out_file);
		  fwrite(&neg_one,sizeof(bool),1,out_file);
		}
	      
	      z*=q_comb;
	      if(neg_one)
		z*=omegas[i]*im_s_vec[s_done].lambda_s_a;
	      else
		z*=omegas[i]*im_s_vec[s_done].lambda_s;

	      if(!contains_zero(z.imag))
		{
		  printf("Non Zero Imag Part %d %d %d\n",q,i,j);
		  print_int_complex_str("z=",z);
		  print_int_complex_str("Pi^(-it/2)=",im_s_vec[s_done].pi_minus_it_2);
		  print_int_double_str("q^(-0.5)=",q_minus_half);
		  //print_int_complex_str("q^(it)=",q_it(q,im_s_vec[s_done].im_s,qits+q));
		  //print_int_complex_str("q^(it/2)=",q_it(q,im_s_vec[s_done].im_s/2.0,qits+q));
		  print_int_complex_str("q_comb=",q_comb);
		  if(neg_one)
		    print_int_complex_str("Gamma((s+1)/2)*exp(Pi*t/4)=",im_s_vec[s_done].lambda_s_a)
		  else
		    print_int_complex_str("Gamma(s/2)*exp(Pi*t/4)=",im_s_vec[s_done].lambda_s)
		  print_int_complex_str("omega=",omegas[i]);
		  printf("Im(s)=%30.28e\n",im_s_vec[s_done].im_s);
		  fatal_error("Exiting.");
		}
	      /*
	      if((i==101)&&(s_done==123))
		{
		  print_int_complex_str("Fft=",Fft_vec[i]);
		  print_int_complex_str("q_comb=",q_comb);
		  print_int_complex_str("Pi^(-it/2)=",im_s_vec[s_done].pi_minus_it_2);
		  print_int_double_str("q^(-0.5)=",q_minus_half);
		  print_int_complex_str("q^(it)=",q_it(q,im_s_vec[s_done].im_s));
		  print_int_complex_str("q^(it/2)=",sqrt(q_it(q,im_s_vec[s_done].im_s)));
		  print_int_double_str("x=",z.real);
		}
	      */	      
	      fwrite(&z.real,sizeof(int_double),1,out_file);
	      
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
      my_hur_ptr1+=num_fracs;
      my_hur_ptr=my_hur_ptr1;
    }
  }
}

int calc_hurwitz1(int_complex &s,unsigned int r_n_terms,unsigned int file_N,double gap,int_complex *s_array,
          int_complex *r_n_vals,unsigned int no_gaps,unsigned int q_start,
		  unsigned int q_end, FILE *out_file)//, qit_t *qits)
{
  int_complex out_val;
  int_complex taylor_err;
  unsigned int i,j,gap_ptr,q,num,steps;
  int_double x,frac,lnq;
  int_complex q_s,first_term;
  double dsteps,csteps,dfloor;

  //printf("in calc_hurwitz1\n");
  create_s_array(s,s_array,r_n_terms,gap); // end column used for error estimation
  taylor_err=taylor_error(s,r_n_terms,gap,s_array);
  //printf("setup completed.\n"); 
  gap_ptr=0;
  
  for(i=0;i<=no_gaps;i++)
  {

    for(j=0;j<r_n_terms;j++)
      r_n_vals[gap_ptr++]=read_int_complex();
    skip_int_complex(file_N-r_n_terms);
  }

  num_fracs=0;
  __builtin_prefetch(s_array,0,3);
  for(q=q_start;q<=q_end;q++)
    {
    if((q&3)==2)
      continue;      // if q = 2 mod 4 then no primitive characters
 
    lnq=log(int_double(q));
    //q_s=n_it(q,s.imag.left)*sqrt(int_double(q));//pow(int_double(q),s.real.left);
    q_s=q_it1(lnq,s.imag.left)//q_it(q,s.imag.left,qits+q)
      *sqrt(int_double(q));

    for(num=1;num<q;num++)
      if(co_prime(num,q))
	{
	  //printf("Doing num=%ld\n",num);
	  __builtin_prefetch(&hur_vec[hur_ptr],1,3);  // going to write

	  // figure out nearest row in lattice
	  dsteps=double(num)/(gap*q);
	  dfloor=floor(dsteps);
	  if((dsteps-dfloor)>0.5)
	    {
	      csteps=dfloor+1.0;
	      steps=(int) csteps;
	      gap_ptr=r_n_terms*(no_gaps-steps);
	    }
	  else
	    {
	      csteps=dfloor;
	      steps=(int) csteps;
	      gap_ptr=r_n_terms*(no_gaps-steps);
	    }
	  __builtin_prefetch(&r_n_vals[gap_ptr],0,0); // going to read, not needed for long
	  num_fracs++;
	  frac=int_double(num)/q;
	  x=int_double(steps)-frac/gap;
	  //print_int_complex_str("n^it=",n_it(num,s.imag.left));
	  //print_int_complex_str("q^s=",q_s);
	  out_val=taylor_err+conj(q_it1(log(int_double(num)),s.imag.left))//q_it(num,s.imag.left,qits+num))
	    *q_s/sqrt(int_double(num));//pow(int_double(num),-0.5)*
	  //printf("out_val abs error = %d, %d.\n",abs_error(out_val.real),abs_error(out_val.imag));
	  out_val+=calc_r_n(r_n_terms,s_array,x,&r_n_vals[gap_ptr]);
	  //printf("out_val abs error = %d, %d.\n",abs_error(out_val.real),abs_error(out_val.imag));
	  //exit(0);
	  // out_val now contains Zeta(s,num/q)
	  hur_vec[hur_ptr++]=out_val;
	}
    }
  return(num_fracs);
}

int main(int argc, char **argv)
{

  // initialsie the int_double package
  _fpu_rndd();

  int_complex s,gam,*s_array,*r_n_vals,*omegas;
  im_s *im_s_vec;
  unsigned int *a_n,*q_n,*offset_vec,num_fracs,q,max_num_prims,i;
  int no_gaps,num_s,N,rn,file_N;
  int_double mgap,im_s_2_mod;
  double gap;
  FILE *in_file,*out_file;//,*qit_file;
  //ifstream facs_file;
  factor *factors;
  clock_t no_clicks,no_clicks1;

  no_clicks=clock(); // start timing
  //printf("At program start.\n");

  printf("Comamnd line:-");
  for(i=0;i<argc;i++) printf(" %s",argv[i]);
  printf("\n");

  if(argc!=5)
    print_usage();
  q=atoi(argv[1]);
  if(q>MAX_Q)
    fatal_error("q exceeds MAX_Q. Exiting.\n");

  if((q<3)||(q>MAX_Q)||((q&3)==2))
    print_usage();
  in_file=fopen(argv[2],"rb");
  if(!in_file)
    {
      perror("Error opening in_file. ");
      printf("Failed to open file |%s|\n",argv[2]);
      fatal_error("Couldn't open in_file. Exiting.\n");
    }


  if(!(factors=(factor *) calloc(q+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  if(!make_factors(factors,q))
    fatal_error("Error reading factor file. Exiting.\n");

  //facs_file.close();
  /*
  qit_file=fopen(argv[5],"rb");
  if(!qit_file)
    {
      perror("Error opening qit_file. ");
      printf("Failed to open file |%s|\n",argv[5]);
      fatal_error("Couldn't open qit file. Exiting.\n");
    }

  qit_t *qits=(qit_t *) malloc(sizeof(qit_t)*MAX_Q+1);
  if(!qits)
    fatal_error("Error allocating memory for qits. Exiting.\n");

  for(i=0;i<=MAX_Q;i++)
    if(fread(qits+i,sizeof(qit_t),1,qit_file)!=1)
      fatal_error("Problem reading qit file. Exiting.\n");
  fclose(qit_file);
  //print_int_complex_str("100^(5*32/128)=",qits[100].bits[5]);
  */
  N=atoi(argv[4]);
  /*
  sem_t *lfuncsem=sem_open(argv[7],0);
  if(lfuncsem==SEM_FAILED)
    {
      perror("");
      fatal_error("Can't open semaphore /lfuncsem. Exiting.\n"); 
    }
  */
  fread(&num_s,sizeof(int),1,in_file);
  fread(&file_N,sizeof(int),1,in_file);
  if((N>file_N)||(N<=0))
    fatal_error("N<=0 or exceeds N used to save lattice file. Exiting.\n");
  fread(&rn,sizeof(int),1,in_file);
  if(rn!=1)
    fatal_error("RN must be 1. Exiting.\n");
  out_file=fopen(argv[3],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");
  fread(&no_gaps,sizeof(int),1,in_file);

  gap=1.0/no_gaps; // this is exact, assuming no_gaps a reasonable power of 2
  //printf("gap set to %10.8e\n",gap);


  //  cout << "Allocating memory for s_array." << endl;
  if(!(s_array=(int_complex *) _aligned_malloc(sizeof(int_complex)*N,16)))
    fatal_error("Can't allocate memory for s_array. Exiting.");
  
  
  //  cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=(int_complex *) _aligned_malloc(sizeof(int_complex)*(no_gaps+1)*N,16)))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");

  //  cout << "Allocating memory for im_s_vec." << endl;
  if(!(im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*num_s,16)))
    fatal_error("Can't allocate memory for im_s_vec. Exiting.");

  if(!(Fft_vec=(int_complex *) _aligned_malloc(sizeof(int_complex)*q,16)))
    fatal_error("Can't allocate memory for Fft_vec. Exiting.");

  fwrite(&q,sizeof(int),1,out_file);
  fwrite(&q,sizeof(int),1,out_file); // start and end q the same

  num_fracs=factors[q].phi;
  max_num_prims=num_prims(q,factors);

  fwrite(&max_num_prims,sizeof(int),1,out_file);
  fwrite(&num_s,sizeof(int),1,out_file);

  s.real=int_double(0.5); // on the 1/2 line

  if(!(hur_vec=(int_complex *) _aligned_malloc(num_fracs*sizeof(int_complex)*num_s,16)))
    fatal_error("Failed to allocate memory for hur_vec. Exiting.");

  // try and force everything resident
  for(i=0;i<num_fracs*num_s;i+=getpagesize()/sizeof(int_complex))
    hur_vec[i]=c_zero;
  if(!(lattice_buff=(double *) malloc((13+4*(no_gaps+1)*file_N)*num_s*sizeof(double))))
    fatal_error("Failed to allocate memory for lattice_buff. Exiting.\n");
  //printf("allocated memory\n");
  /*
  sem_wait(lfuncsem);
  */
  //printf("remember to uncomment sem_wait etc.\n");
  fread(lattice_buff,sizeof(double),(13+4*(no_gaps+1)*file_N)*num_s,in_file);
  //printf("read lattice file\n");
  fclose(in_file);
  /*
  sem_post(lfuncsem);
  */
  for(i=0;i<num_s;i++)
    {
      im_s_vec[i].im_s=read_double();
      //printf("im_s=%10.8e\n",im_s_vec[i].im_s);
      s.imag=int_double(im_s_vec[i].im_s);
      im_s_vec[i].lambda_s=read_int_complex(); // GAMMA(s/2)*exp(Pi t/4)
      im_s_vec[i].lambda_s_a=read_int_complex(); // GAMMA((s+1)/2)*exp(Pi t/4)
      im_s_vec[i].pi_minus_it_2=read_int_complex();
      //print_int_complex_str("gamma(s/2)*exp(Pi t/4)=",im_s_vec[i].lambda_s);exit(0);
      num_fracs=calc_hurwitz1(s,N,file_N,gap,s_array,r_n_vals,no_gaps,q,q,out_file);//,qits);
      //printf("calc_hurwitz1 completed.\n");
      if(!num_fracs)
	fatal_error("Error running Hurwitz routine. Exiting.");      
    }

  free(lattice_buff);
  _aligned_free(r_n_vals);
  _aligned_free(s_array);
  
  //printf("Done %ld fracs and %ld s vals from Im(s)=%10.8e, Q=%ld.\n",num_fracs,num_s,im_s_vec[0].im_s,q);
  
  printf("Hur Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  no_clicks1=clock(); // start timing again
  
  //fwrite(im_s_vec,sizeof(im_s),num_s,out_file);

  if(!(omegas=(int_complex *) _aligned_malloc((q-1)*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(q-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(q-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(q-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  init_ws(_ws);
  hur_ptr=0;

  make_l(q,num_s,factors,q_n,im_s_vec,offset_vec,a_n,omegas,out_file);//,qits);

  fclose(out_file);

  printf("FFT Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks1)/((double) CLOCKS_PER_SEC));
  printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  //printf("Lambda value straddled zero %d times\n",temp_count);
  printf("l-func successful completion on q=%d-%d file %s\n",q,q,argv[3]);
  return(0);
}
