//
// upsamdouble.cpp
//

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <complex>
#include <limits.h>
#include "../includes/int_double12.0.h"
//#include "../includes/im_s.h"
#include "../includes/fft_defs_half.h"
#include "../includes/make-factors.h"
//#include "../includes/int-fft.h"
//#include "../includes/upsamdefs.h"

#define DEBUG printf("Reached line %d\n",__LINE__)
// q is 2^a, a>=3, phi = q/4

int_complex chis[MAX_Q];

inline long int ln2(int q)
{
  long int res=0;
  while(q)
    {
      q>>=1;
      res++;
    }
  return(res);
}

inline bool power_2_p(unsigned long int i)
{
  while((i&1)==0)
    i>>=1;
  return(i==1);
}

// q is 2^a with a>=3
void make_chis1(unsigned int q, unsigned int index, unsigned int pr, unsigned int phi, int_complex *chis)
{
  unsigned int pow,chi_ptr,n,phi_by_2=phi>>1;
  int even_p=(index<phi_by_2);
  int_double theta,four_pi_over_phi=int_double(4)/phi;//d_two_pi/(phi/2);
  unsigned int w=0;
  //printf("q=%d\nln2(q)=%d\n",q,ln2(q));

  //if((ln2(q)&1)==1)
  // still not got this right..
  if (q==16) four_pi_over_phi=-four_pi_over_phi;
  chi_ptr=5;
  chis[1]=c_one;
  if(even_p)
    chis[q-1]=c_one;
  else
    chis[q-1]=-c_one;
  for(pow=1;pow<phi_by_2;pow++)
    {
      w+=index;
      while(w>=phi_by_2)
	w-=phi_by_2;
      theta=four_pi_over_phi*w;
      sin_cospi(theta,&chis[chi_ptr].imag,&chis[chi_ptr].real);
      if(even_p)
	chis[q-chi_ptr]=chis[chi_ptr];
      else
	chis[q-chi_ptr]=-chis[chi_ptr];
      chi_ptr=(chi_ptr*5)%q;
    }
  return;
}

void make_chis(unsigned int q, unsigned int index, unsigned int pr, unsigned int phi, unsigned int no_dims, int_complex *chis)
{
  unsigned long int n,a,chi_ptr,w=0;
  int even_p;
  int_double two_pi_over_phi,theta;
  //printf("%d\n",sizeof(l));
  //exit(0);
  //printf("make_chis called with q=%d index=%d pr=%d phi=%d\n",
  // q,index,pr,phi);
  if(pr==0)
    {
      make_chis1(q,index,pr,phi,chis);
      return;
    }
  chis[1]=c_one;
  even_p=((index&1)==0);
  //if (even_p) printf("even character\n"); else printf("odd character\n");
  if (even_p)
    chis[q-1]=c_one;
  else
    chis[q-1]=-c_one;
  if(phi==2)
    return;
  a=pr;
  two_pi_over_phi=-d_two_pi/phi;
  // now I cocked up my initial DFT so that all
  // n-dimensional DFT's <=50 not 2 or 4 have sum (exp(2pi i/n) not exp(-2pi i/n)
  // 1-d DFTs are OK because they just use Bluestein. Power of 2 work half the time.
  if((no_dims>1)&&(phi>4)) // single dimension and length 2,4 work fine
    {
      if(power_2_p(phi))
	{
	  //printf("phi=%d\n",phi);
	  //	  if((ln2(phi)&1)==0)
	    two_pi_over_phi=-two_pi_over_phi; // phi=16,64,256 ...
	}
      else
	{
	  if(phi<=MAX_SIMPLE_DFT)             // phi=6,10,12,14,18...MAX
	    two_pi_over_phi=-two_pi_over_phi;
	}
    }

    /*
  if((no_dims==1)||(phi>MAX_SIMPLE_DFT)||(phi==2)||(phi==4)||(power_2_p(phi)&&((ln2(phi)&1)==0)))//||(phi==8)||(phi==16)||(phi==32))
    two_pi_over_phi=-two_pi_over_phi;
    */
  for(n=1;n<phi/2;n++)
    {
      w+=index;
      if(w>=phi)
	w-=phi;
      theta=two_pi_over_phi*w;
      //print_int_double_str("theta=",theta);
      sin_cos(theta,&chis[a].imag,&chis[a].real);
      if(even_p)
	chis[q-a]=chis[a];
      else
	chis[q-a]=-chis[a];
      a=(a*pr)%q;
    }
  return;
}

inline void make_coords(long int *coords, long int q, long int index, factor *factors)
{
  long int n=index,dim;
  for(dim=factors[q].num_facs-1;dim>=0;dim--)
    {
      coords[dim]=n%factors[factors[q].facs[dim]].phi;
      n/=factors[factors[q].facs[dim]].phi;
    }
  //for(dim=0;dim<factors[q].num_facs;dim++)
  //printf("%d ",coords[dim]);
  //printf("\n");
  return;
}


int_complex L(double re_z, double im_z, long int q, long int index, factor *factors, int_complex *s_array, long int *coords)
{ 
  long int chi_ptr=0;
  long int n,a,q_phi=factors[q].phi,this_phi,this_pr,dim,this_q;
  int_complex hur,res,minus_s=int_complex(int_double(-re_z),int_double(-im_z));
  make_coords(coords,q,index,factors);
	
  for(dim=0;dim<factors[q].num_facs;dim++)
    {
      this_q=factors[q].facs[dim];
      this_pr=factors[this_q].pr;
      this_phi=factors[this_q].phi;
      make_chis(this_q,coords[dim],
		this_pr,this_phi,factors[q].num_facs,&chis[chi_ptr]);
      chi_ptr+=this_q;
    }
  //printf("chis made.\n");for(n=0,a=10;n<20;n++,a=(a*10)%q) {printf("chi[%ld]=",n);print_int_complex_str("",chis[n]);}
  //exit(0);

  res=c_zero;

  for(n=1;n<q;n++)
    {
      //if(n==35) exit(0);
      if(co_prime(n,q))
	{	
	  hur.real=hurwitz_half(int_double(n)/q,s_array);
	  hur.imag=d_zero;
	  //printf("Zeta(s,%d/%d)=\n",n,q);print_int_complex_str("",hur);
	  chi_ptr=0;
	  for(dim=0;dim<factors[q].num_facs;dim++)
	    {
	      hur*=chis[chi_ptr+(n%factors[q].facs[dim])]; // Chinese Remainder Theorem
	      //printf("chi(%d,%d)=",n,dim);print_int_complex_str("",chis[chi_ptr+(n%factors[q].facs[dim])]);
	      chi_ptr+=factors[q].facs[dim];
	    }
	  //printf("Zeta(s,%d/%d)*chi(%d)=",n,q,n);mpfi_c_print(L_hur);

	  res+=hur;
	}
    }
  res*=pow(q,minus_s);
  return(res);
}

int_complex L_same_q(double re_z, double im_z, long int q, long int index, factor *factors, int_complex *s_array, long int *coords)
{ 
  long int chi_ptr=0;
  long int n,a,q_phi=factors[q].phi,this_phi,this_pr,dim,this_q;
  int_complex hur,res,minus_s=int_complex(int_double(-re_z),int_double(-im_z));

  res=c_zero;
	
  for(n=1;n<q;n++)
    {
      //if(n==35) exit(0);
      if(co_prime(n,q))
	{	
	  hur.real=hurwitz_half(int_double(n)/q,s_array);
	  hur.imag=d_zero;
	  //printf("Zeta(s,%d/%d)=\n",n,q);print_int_complex_str("",hur);
	  chi_ptr=0;
	  for(dim=0;dim<factors[q].num_facs;dim++)
	    {
	      hur*=chis[chi_ptr+(n%factors[q].facs[dim])]; // Chinese Remainder Theorem
	      //printf("chi(%d,%d)=",n,dim);print_int_complex_str("",chis[chi_ptr+(n%factors[q].facs[dim])]);
	      chi_ptr+=factors[q].facs[dim];
	    }
	  //printf("Zeta(s,%d/%d)*chi(%d)=",n,q,n);mpfi_c_print(L_hur);

	  res+=hur;
	}
    }
  res*=pow(q,minus_s);
  return(res);
}

long int fix_index(long int index, long int phi)
{
  if(index<phi/4)
    return(index+phi/4);
  if(index<phi/2)
    return(index-phi/4);
  if(index<3*phi/4)
    return(index+phi/4);
  return(index-phi/4);
}

int main(int argc, char **argv)
{

  long int prec,q,index,i,index1;
  long int coords[MAX_DIMS];
  factor *factors;
  int_complex s_array[MAX_BERNOULLI_K],res;
  FILE *infile;

  _fpu_rndd();
  set_h_bernoulli();


  /*  check all the command line arguments are ok, if not print message
      and exit sharpish */


  if(argc!=2)
    {
      printf("usage:- l_half_double <data file>\nExiting.\n");
      exit(1);
    }

  infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Failed to open text file %s for input. Exiting.\n",argv[1]);
      exit(1);
    }


  factors=(factor *) malloc(sizeof(factor)*(MAX_Q+1));
  if(!factors)
    {
      printf("Fatal error allocating memory for factors. Exiting.\n");
      exit(1);
    }


  if(!make_factors(factors,MAX_Q))
    {
      printf("Fatal error building factors. Exiting\n");
      exit(1);
    }

   hur_init(s_array,0.5,0);
   while(true)
     {
       if(fscanf(infile,"%ld %ld\n",&q,&index)!=2)
	 return(0);
       //printf("processing q=%ld index=%ld\n",q,index);
       if((q>MAX_Q)||(index>factors[q].phi))
	 {
	   printf("q,index bad. Exiting.\n");
	   exit(1);
	 }
       /*
       if((factors[factors[q].phi].primes[0]==2)&&(factors[factors[q].phi].facs[0]>=16))
	 index1=fix_index(index,factors[q].phi);
       else
       */
	 index1=index;
       
       res=L(0.5,0.0,q,index1,factors,s_array,coords);

       if(contains_zero(res))
	 {
	   printf("L(1/2,chi_%d_%d) contains zero. Try at higher precison.",q,index);
	   print_int_complex_str("result was ",res);
	 }
       else
	 {
	   printf("L(1/2,chi_%d_%d) non-zero.",q,index);
	   print_int_complex_str("result was ",res);
	 }
     }
   //printf("Completed.\n");
  return(0);
}

