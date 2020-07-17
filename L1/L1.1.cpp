/*
File: L1.1.cpp

Created: 28th June 2011

Version: <v> = 1.1

Dialect: C++

Requires: -lrt

Implementation notes:
            num_s (int) = 1
            N (int) = 5 no of Taylor terms
            rn (int) must be 2
            NO_GAPS (int)

            N*(NO_GAPS+1) h_rn (dcomplex)

Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

By: DJ Platt
    Bristol University

Copyright 2008.

lattice file contains -gamma-Psi(3-i/4096),-Psi'(3-i/4096)/1!,-Psi''(3-i/4096)/2!,-Psi'''(3-i/4096)/3!,-Psi''''(3-i/4096)/4!
for i=0..4096


This work is funded by the UK ESPRC. */
#define LINUX
#define VERSION "1.1"
#define TAY_ERR ((double) 1.6e-20) // Take 5 taylor terms and RN=2

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
#include "../includes/im_s.h"
#include "../includes/int-fft-half.h"
//#include "L1.h"

void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: int-l-func%s (q-start) (q-end) (ifname) (ofname)\n",VERSION);
  printf("  (q-start)   - integer >=3\n");
  printf("  (q-end)     - integer >=(q-start)\n");
  printf("  (ifname)    - file with lattice values.\n");
  printf("  (ofname)    - output file.\n");
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

inline int_double calc_r_n(unsigned int r_n_terms, int_double &delta,
      int_double *r_n_vec)
{
  unsigned int j,i=r_n_terms-2;
  int_double res=r_n_vec[i+1];
  //print_int_double_str("in calc_r_n with delta=",delta);
  //for(j=0;j<r_n_terms;j++)
  //print_int_double_str("r_n[j]=",r_n_vec[j]);
	for(;i>0;i--)
	  {
	    res=res*delta+r_n_vec[i];
	    //print_int_double_str("next term=",r_n_vec[i]);
	    //print_int_double_str("res=",res);
	  }
	res=res*delta+r_n_vec[0];
	//print_int_double_str("last term=",r_n_vec[0]);
	//print_int_double_str("res=",res);

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
inline void skip_int_double(unsigned int n)
{
  lattice_buff_ptr+=n*sizeof(int_double);
}

void init_fft_nd(unsigned int *offsets,unsigned int q, int_complex *omegas,unsigned int n_dims,int *dims,
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


unsigned int hur_ptr=0,num_fracs;
int_complex *hur_vec;

//unsigned int temp_count=0;

int_double big_so_far=d_zero;

void make_l(unsigned int q, 
      unsigned int num_s, 
      factor *factors,
      unsigned int *q_n,
      im_s *im_s_vec,
      unsigned int *a_n,
      unsigned int *offset_vec,
      int_complex *omegas,
      FILE *out_file)
{
  unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS];
  unsigned int my_hur_ptr,my_hur_ptr1;
  int_complex omega,omega_a,z,z1,q_comb;
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
    //printf("pr of %ld is %ld\n",q,pr);
    fill_q_ns(q_n,pr,phi_q,q);
    init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
    //prep_omegas(omegas,q,phi_q,q_n,a,bs,b_star_conjs,conv_sizes[0],b_spare);
    //for(i=0;i<phi_q;i++)
    //  print_int_complex_str("prepped omegas[]= ",omegas[i]);
    s_done=0;
    my_hur_ptr1=hur_ptr;
    my_hur_ptr=hur_ptr;
    //printf("my_hur_ptr %ld\n",my_hur_ptr);
    for(i=1;i<q;i++)
      if(co_prime(i,q))
	Fft_vec[q_n[i]]=hur_vec[my_hur_ptr++];
      /*
      for(i=0;i<phi_q;i++)
	{
	  printf("Zeta(1/2,a/%ld)=",q);
	  print_int_complex_str("",Fft_vec[i]);
	}
      */

      bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare);
      // now fft_vec contains L(chi,1)*q

      int_double smallest_odd,smallest_even,largest_odd,largest_even;
      smallest_even=DBL_MAX;
      largest_even=0;
      smallest_odd=DBL_MAX;
      largest_odd=0;
      for(i=1,k=0;i<phi_q;i++)  // i=0 corresponds to principal chi
	{
	  if(prim_p(i,factors[q].primes[0])) // only do primitive charcters
	    {
	      //printf("k=%lu, conj(k)=%lu\n",k,conj_j(k,n_prims,q));
	      if(conj_j(k,n_prims,q)>=k) // pair conjugates
		{
		  x=mod(Fft_vec[i])/q;
		  //printf("%7u %7u ",q,i);if(neg_one) printf("odd  "); else printf("even ");
		  //print_int_double_str("",x);

		  //print_int_double_str("     |L(1,chi)|=",x);
		  if(i&1) // its an odd character
		    {
		      if(x.left<smallest_odd.left) // left end points held as is
			smallest_odd=x;
		      if(x.right<largest_odd.right) // right end points held negated
			largest_odd=x;
		    }
		  else
		    {
		      if(x.left<smallest_even.left)
			smallest_even=x;
		      if(x.right<largest_even.right)
			largest_even=x;
		    }
		}
	    
	      k++;
	    }
	}
      fwrite(&q,sizeof(unsigned int),1,out_file);
      fwrite(&largest_odd,sizeof(int_double),1,out_file);
      fwrite(&smallest_odd,sizeof(int_double),1,out_file);
      fwrite(&largest_even,sizeof(int_double),1,out_file);
      fwrite(&smallest_even,sizeof(int_double),1,out_file);
      /*
      printf("Q=%lu\n",q);
      print_int_double_str("Largest Odd  =",largest_odd);
      print_int_double_str("Smallest Odd =",smallest_odd);
      if(q>4)
	{
	  print_int_double_str("Largest Even =",largest_even);
	  print_int_double_str("Smallest Even =",smallest_even);
	}
      */
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
    init_fft_nd(offset_vec,q,omegas,no_dims,dims,phi_q);
    //printf("omegas prepped.\n");
    my_hur_ptr1=hur_ptr;
    my_hur_ptr=hur_ptr;
//    printf("my_hur_ptr %ld\n",my_hur_ptr);
    //while(s_done<num_s)
    //{
      //sin_cos(-log(int_double(q))*(double) (im_s_vec[s_done].im_s/2.0),&q_minus_it_2.imag,&q_minus_it_2.real);
      //q_comb=q_minus_it(q,im_s_vec[s_done].im_s/2)*im_s_vec[s_done].pi_minus_it_2*pow(int_double(q),-0.5);
      //print_int_complex_str("q^it=",q_it(q,im_s_vec[s_done].im_s/2.0));
      //q_comb=im_s_vec[s_done].pi_minus_it_2*q_minus_half;//*conj(q_it(q,im_s_vec[s_done].im_s/2.0));
      for(i=0;i<phi_q;i++)
	Fft_vec[offset_vec[i]]=hur_vec[my_hur_ptr++];
      //printf("Doing FFT.\n");
      do_nd_fft(Fft_vec,no_dims,dims,phi_q);
      //printf("FFT Done\n");
      for(i=0;i<factors[q].num_facs;i++)
        coords[i]=0;
      int_double smallest_odd,smallest_even,largest_odd,largest_even;
      smallest_even=DBL_MAX;
      largest_even=0;
      smallest_odd=DBL_MAX;
      largest_odd=0;

      for(i=0,k=0;i<phi_q;i++)
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
	      if(conj_j(k,n_prims,q)>=k)
		{
		  x=mod(Fft_vec[i])/q;
		  neg_one=neg_one_p(coords,factors[q].num_facs);
		  if(power_2&&(i<(phi_q>>1)))
		    neg_one=!neg_one;
		  //printf("%7u %7u ",q,i);if(neg_one) printf("odd  "); else printf("even ");
		  //print_int_double_str("",x);
		  if(neg_one) // its an odd character
		    {
		      if(x.left<smallest_odd.left) // left end points held as is
			smallest_odd=x;
		      if(x.right<largest_odd.right) // right end points held negated
			largest_odd=x;
		    }
		  else
		    {
		      if(x.left<smallest_even.left)
			smallest_even=x;
		      if(x.right<largest_even.right)
			largest_even=x;
		    }

		}
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
      fwrite(&q,sizeof(unsigned int),1,out_file);
      fwrite(&largest_odd,sizeof(int_double),1,out_file);
      fwrite(&smallest_odd,sizeof(int_double),1,out_file);
      fwrite(&largest_even,sizeof(int_double),1,out_file);
      fwrite(&smallest_even,sizeof(int_double),1,out_file);
      /*
      printf("Q=%lu\n",q);
      print_int_double_str("Largest Odd  =",largest_odd);
      print_int_double_str("Smallest Odd =",smallest_odd);
      if(q>4)
	{
	  print_int_double_str("Largest Even =",largest_even);
	  print_int_double_str("Smallest Even =",smallest_even);
	}
      */
  }
}

int calc_hurwitz1(unsigned int r_n_terms,unsigned int file_N,double gap,
          int_double *r_n_vals,unsigned int no_gaps,unsigned int q_start,
		  unsigned int q_end, FILE *out_file,int_double &taylor_err)
{
  int_double out_val;
  unsigned int i,j,q,num,steps,gap_ptr;
  int_double x,frac;
  int_double first_term;
  double dsteps,csteps,dfloor;
  long double lnq;

  //printf("in calc_hurwitz1\n");
  //printf("setup completed.\n");
  

  num_fracs=0;
  for(q=q_start;q<=q_end;q++)
    {
      if((q&3)==2)
	continue;      // if q = 2 mod 4 then no primitive characters

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
	    x=frac-(double)steps/no_gaps;
	    //printf("step=%d\n",steps);
	    out_val=calc_r_n(r_n_terms,x,&r_n_vals[gap_ptr]);
	    //printf("zeta_2(1,%ld/%ld)-zeta(1)=",num,q);print_int_double_str("",out_val);
	    out_val+=taylor_err+d_one/(d_one+frac)+d_one/frac;
	    // out_val now contains Hurwitz(1,num/q)-zeta(1)
	    //printf("zeta(1,%ld/%ld)-zeta(1)=",num,q);print_int_double_str("",out_val);
	    hur_vec[hur_ptr++]=out_val;
	  }
    }
  return(num_fracs);
}

int main(int argc, char **argv)
{

  _fpu_rndd();

  int_complex s,gam,*s_array,*omegas;
  im_s im_s_vec[1];
  unsigned int *a_n,*q_n,*offset_vec,num_fracs,q,max_num_prims,n_prim,i;
  int no_gaps,q_start,q_end,num_s,N,rn,file_N;
  int_double mgap,im_s_2_mod,*r_n_vals;
  double gap;
  FILE *in_file,*out_file;
  factor *factors;

  if(argc!=5)
    print_usage();
  q_start=atoi(argv[1]);
  if(q_start<3)
    q_start=3;
  q_end=atoi(argv[2]);
  if(q_end<q_start)
    q_start=q_end;
  
  if(q_end>MAX_Q)
    {
      printf("q_end (%ld) exceeds MAX_Q (%ld). Exiting.\n",q_end,MAX_Q);
      exit(0);
    }
  
  in_file=fopen(argv[3],"rb");
  if(!in_file)
    {
      perror("Error opening in_file. ");
      printf("Failed to open file |%s|\n",argv[3]);
      fatal_error("Couldn't open in_file. Exiting.\n");
    }

  out_file=fopen(argv[4],"wb");
  if(!out_file)
    {
      perror("Error opening out_file. ");
      printf("Failed to open file |%s|\n",argv[3]);
      fatal_error("Couldn't open in_file. Exiting.\n");
    }
  fwrite(&q_start,sizeof(int),1,out_file);
  fwrite(&q_end,sizeof(int),1,out_file);

  if(!(factors=(factor *) calloc(q_end+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  //printf("Making factors.\n");
  if(!make_factors(factors,q_end))
    fatal_error("Error creating factors. Exiting.\n");

  N=5; // Use 5 Taylor terms

  fread(&num_s,sizeof(int),1,in_file);
  assert(num_s==1);
  fread(&file_N,sizeof(int),1,in_file);
  if((N>file_N)||(N<=0))
    fatal_error("N<=0 or exceeds N used to save lattice file. Exiting.\n");
  fread(&rn,sizeof(int),1,in_file);
  if(rn!=2)
    fatal_error("RN must be 2. Exiting.\n");
  /*
  out_file=fopen(argv[4],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");
  */
  fread(&no_gaps,sizeof(int),1,in_file);

  gap=1.0/no_gaps;
  //printf("gap set to %10.8e\n",gap);


  //cout << "Allocating memory for s_array." << endl;
  //if(!(s_array=(int_complex *) _aligned_malloc(sizeof(int_complex)*N,16)))
  //fatal_error("Can't allocate memory for s_array. Exiting.");
  
  
  //cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=(int_double *) _aligned_malloc(sizeof(int_double)*(no_gaps+1)*N,16)))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");


  //fwrite(&q_start,sizeof(int),1,out_file);
  //fwrite(&q_end,sizeof(int),1,out_file);


  //fwrite(&max_num_prims,sizeof(int),1,out_file);
  //fwrite(&num_s,sizeof(int),1,out_file);

  s.real=int_double(1.0);

  if(!(hur_vec=(int_complex *) _aligned_malloc(q_end*sizeof(int_complex)*num_s,16)))
    fatal_error("Failed to allocate memory for hur_vec. Exiting.");

  // try and force everything resident
  /*
  for(i=0;i<num_fracs*num_s;i+=getpagesize()/sizeof(int_complex))
    hur_vec[i]=c_zero;
  */
  if(!(lattice_buff=(double *) malloc((2*(no_gaps+1)*file_N)*sizeof(double))))
    fatal_error("Failed to allocate memory for lattice_buff. Exiting.\n");
  printf("allocated memory\n");
  
  //sem_wait(lfuncsem);
  
  //printf("remember to uncomment sem_wait etc.\n");
  fread(lattice_buff,sizeof(double),(2*(no_gaps+1)*file_N),in_file);
  printf("read lattice file\n");
  fclose(in_file);
  
  //if(!(omegas=(int_complex *) _aligned_malloc((q_end-1)*sizeof(int_complex),16)))
  //fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  printf("initialising ws\n");
  init_ws();
  printf("ws initialised\n");
  /*
  for(i=3;i<100;i++)
    printf("%lu has %lu primitive characters.\n",i,num_prims(i,factors));
  exit(0);
  */
  //create_s_array(s,s_array,N,gap); // end column used for error estimation
  int_double taylor_err=int_double(-TAY_ERR,TAY_ERR);
  int j,gap_ptr=0;
  for(i=0;i<=no_gaps;i++)
  {

    for(j=0;j<N;j++)
      r_n_vals[gap_ptr++]=read_int_double();
    skip_int_double(file_N-N);
  }
  /*
  for(int i=0;i<10;i++)
    {
      printf("rn_vals[%d]",i);
      print_int_double_str("=",r_n_vals[i]);
    }
  exit(0);
  */
  printf("Only looking at primitive characters, q=0 mod 3.\n");
  for(q=q_start;q<=q_end;q++)
    if((q%3==0)&&((q&3)!=2))
      {
	num_fracs=calc_hurwitz1(N,file_N,gap,r_n_vals,no_gaps,q,q,out_file,taylor_err);
	//printf("calc_hurwitz1 completed.\n");
	//exit(0);
	if(!num_fracs)
	  fatal_error("Error running Hurwitz routine. Exiting.");
	hur_ptr=0;
	      
	//printf("Processing q=%ld\n",q);
	make_l(q,num_s,factors,q_n,im_s_vec,offset_vec,a_n,omegas,out_file);
      }

  //fclose(out_file);

  //printf("FFT Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks1)/((double) CLOCKS_PER_SEC));
  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  //printf("Lambda value straddled zero %d times\n",temp_count);
  printf("l-func successful completion on q=%d-%d file %s\n",q_start,q_end,argv[3]);


  return(0);
}
