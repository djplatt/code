/*

File: upsamhigh.c

Created: 2 November 2009

Version: 0.0

Last Modified:

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 0.0 Initial implementation

Build instructions: gcc -o upsamhigh upsamhigh.c -O2 -lmpfi -lmpfr -lgmp -lm

By: DJ Platt
    Bristol University

Copyright 2008,2009.

This work is funded by the UK ESPRC. */


#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "../includes/mpfi.h"
#include "../includes/mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "../includes/fft_defs_half.h"
#include "../includes/make-factors.h"

#define SUCCESS 0
#define QUIT_CODE 0



void print_usage()
{
  printf("Usage:- ./upsamhigh <prec> <q> <index>\n");
}



mpfi_t L_a,L_two_pi_over_phi,L_theta;
mpfi_c_t L_hur,L_omega;//,L_omega1;
//mpfi_c_t L_omegas[LN2_MAX_Q]; //
mpfi_c_t L_chis[MAX_Q];

void mpfi_c_init_L()
{
	int i;
	mpfi_init(L_a);
	mpfi_init(L_two_pi_over_phi);
	mpfi_init(L_theta);
	mpfi_c_init(L_hur);
	mpfi_c_init(L_omega);
/*
	mpfi_c_init(L_omega1);
	for(i=0;i<LN2_MAX_Q;i++)
		mpfi_c_init(L_omegas[i]);
*/
	for(i=0;i<MAX_Q;i++)
		mpfi_c_init(L_chis[i]);
}

// q is 2^a, a>=3, phi = q/4
void make_chis1(unsigned int q, unsigned int index, unsigned int pr, unsigned int phi, mpfi_c_t *chis)
{
	unsigned int pow,chi_ptr,n;
	int even_p=(index<phi/2);
	mpfi_const_pi(L_two_pi_over_phi);
    mpfi_div_ui(L_two_pi_over_phi,L_two_pi_over_phi,phi/4);
    
    if(q==16) mpfi_neg(L_two_pi_over_phi,L_two_pi_over_phi);
	chi_ptr=5;
	mpfi_c_set_ui(chis[1],1,0);
	if(even_p)
		mpfi_c_set_ui(chis[q-1],1,0);
	else
		mpfi_c_set_d(chis[q-1],-1.0,0.0);

	for(pow=1;pow<q/4;pow++)
	{
		mpfi_mul_ui(L_theta,L_two_pi_over_phi,(pow*index)%(q/4));
		mpfi_cos(chis[chi_ptr]->re,L_theta);
		mpfi_sin(chis[chi_ptr]->im,L_theta);
		if(even_p)
			mpfi_c_set(chis[q-chi_ptr],chis[chi_ptr]);
		else
			mpfi_c_neg(chis[q-chi_ptr],chis[chi_ptr]);
		chi_ptr=(chi_ptr*5)%q;
	}
	return;
		for(n=0;n<q;n++)
	{
		printf("chi[%d]=",n);mpfi_c_print(chis[n]);
	}

}

inline int power_2_p(long int i)
{
  while((i&1)==0)
    i>>=1;
  return(i==1);
}

void make_chis(unsigned long int q, unsigned long int index, unsigned long int pr, unsigned long int phi, mpfi_c_t *chis, unsigned long int no_dims)
{
  unsigned long int n,a,w=0;
  long int even_p;
	//printf("make_chis called with q=%d index=%d pr=%d phi=%d\n",
	//q,index,pr,phi);
	if(pr==0)
	{
		make_chis1(q,index,pr,phi,chis);
			return;
	}
	mpfi_c_set(chis[1],c_one);
	even_p=((index&1)==0);
	//if (even_p) printf("even character\n"); else printf("odd character\n");

	if (even_p)
		mpfi_c_set(chis[q-1],c_one);
	else
		mpfi_c_set(chis[q-1],c_minus_one);
	if(phi==2)
		return;
	a=pr;
    mpfi_const_pi(L_two_pi_over_phi);
    mpfi_div_d(L_two_pi_over_phi,L_two_pi_over_phi,-(double)phi/2.0); // phi is always even
		// now I cocked up my initial DFT so that all
		// n-dimensional DFT's <=MAX_SIMPLE_DFT not a power of 2 have sum (exp(2pi i/n) not exp(-2pi i/n)
	    // 1-d DFTs are OK because they just use Bluestein

    if((no_dims>1)&&(phi>4)&&(power_2_p(phi)||(phi<=MAX_SIMPLE_DFT)))
      mpfi_neg(L_two_pi_over_phi,L_two_pi_over_phi);
    /*
    if((no_dims==1)||(phi>MAX_SIMPLE_DFT)||(phi==2)||(phi==4)||(phi==8)||(phi==16)||(phi==32))
		mpfi_neg(L_two_pi_over_phi,L_two_pi_over_phi);
    */
	for(n=1;n<phi/2;n++)
	{
	  w+=index;
	  if(w>=phi)
	    w-=phi;
	  mpfi_mul_ui(L_theta,L_two_pi_over_phi,w);
	  mpfi_cos(chis[a]->re,L_theta);
	  mpfi_sin(chis[a]->im,L_theta);
	  if(even_p)
		  mpfi_c_set(chis[q-a],chis[a]);
	  else
		  mpfi_c_neg(chis[q-a],chis[a]);
	  a=(a*pr)%q;
	}
	return;
	for(n=0;n<q;n++)
	{
		printf("chi=");mpfi_c_print(chis[n]);
	}
	
}

inline int conv_index(long int index, long int q)
{
	return(index);
}

void mpfi_c_L(mpfi_c_ptr res, double re_z, double im_z, long int q, long int index, factor *factors)
{ 
  //int even_p=((index&1)==0);
	long int chi_ptr=0;
	long int n,a,q_phi=factors[q].phi,this_phi,this_pr,dim,this_q;
	long int coords[MAX_DIMS];
	n=conv_index(index,q);
	for(dim=factors[q].num_facs-1;dim>=0;dim--)
	{
		coords[dim]=n%factors[factors[q].facs[dim]].phi;
		n/=factors[factors[q].facs[dim]].phi;
	}
	printf("coords=[");
	for(dim=0;dim<factors[q].num_facs;dim++)
		printf(" %d",coords[dim]);
	printf(" ]\n");

	for(dim=0;dim<factors[q].num_facs;dim++)
	{
	  this_q=factors[q].facs[dim];
	  this_pr=factors[this_q].pr;
	  this_phi=factors[this_q].phi;
	  make_chis(this_q,coords[dim],
		    this_pr,this_phi,&L_chis[chi_ptr],factors[q].num_facs);
	  chi_ptr+=this_q;
	}
	
	//printf("Chis calculated.\n");for(n=0;n<20;n++) mpfi_c_print_str("chi=",L_chis[n]);
	//printf("chis made.\n");for(n=0,a=10;n<20;n++,a=(a*10)%q) {printf("chi[%ld]=",n);mpfi_c_print_str("",L_chis[n]);}
	//exit(0);

	mpfi_c_set(res,c_zero);
	
	for(n=1;n<q;n++)
	{
	  if(co_prime(n,q))
	    {	
	      mpfi_set_ui(L_a,n);
	      mpfi_div_ui(L_a,L_a,q);
	      mpfi_c_hurwitz1(L_hur,re_z,im_z,L_a);
	      //printf("Zeta(s,%d/%d)=\n",n,q);mpfi_c_print(L_hur);
	      chi_ptr=0;
	      for(dim=0;dim<factors[q].num_facs;dim++)
		{
		  mpfi_c_mul(L_hur,L_hur,L_chis[chi_ptr+(n%factors[q].facs[dim])]);
		  //printf("chi(%d,%d)=",n,dim);mpfi_c_print(L_chis[chi_ptr+(n%factors[q].facs[dim])]);
		  chi_ptr+=factors[q].facs[dim];
		}
	      //printf("Zeta(s,%d/%d)*chi(%d)=",n,q,n);mpfi_c_print(L_hur);
	      
	      mpfi_c_add(res,res,L_hur);
	    }
	}
	mpfi_c_pow_double_to_doubles(L_hur,q,-re_z,-im_z); // q^(-s)
	//mpfi_c_print(L_hur);
	mpfi_c_mul(res,res,L_hur);
	//printf("L(s,chi)=");mpfi_c_print(res);
	return;
}


int main(int argc, char **argv)
{

  int prec,q,index,i;
  factor *factors;

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=4)
    {
	print_usage();
	return(QUIT_CODE);
    }

    prec=atoi(argv[1]);
    if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    {
	print_usage();
	return(QUIT_CODE);
    }
    q=atoi(argv[2]);
    index=atoi(argv[3]);
    if((q<3)||(index>=q))
      {
	print_usage();
	return(QUIT_CODE);
      }
    if(!(factors=(factor *) malloc(sizeof(factor)*(q+1))))
      {
	printf("Failed to allocate memory for factors. Exiting.\n");
	return(QUIT_CODE);
      }
    make_factors(factors,q);
    

    mpfi_c_setup(prec);
    mpfi_c_init_L();
    double im_s=0.0;
    mpfi_c_t res;
    mpfi_c_init(res);
    
    mpfi_c_hur_init(0.5,0.0);
    mpfi_c_L(res,0.5,0.0,q,index,factors);
    printf("L(1/2,chi_%ld_%ld)=",q,index);mpfi_c_print_str("",res);
    printf("Completed.\n");
    return(QUIT_CODE);
}
