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
#include "../includes/fft_defs.h"
#include "../includes/upsamdefs.h"

#define SUCCESS 0
#define QUIT_CODE 0
#define FAILURE 1

#define TRUE (1==1)
#define FALSE (0==1)

// read the factors file
// file starts at 3
// each line consists of:-
// primitive root (0 if there isn't one)
// phi(n)
// num_facs = number of factors of the form p^n with p prime
// the factors from lowest p up in the form p then p^n
int read_factors(factor *factors, unsigned int n, FILE *facs_file)
{
  unsigned int i,j,k;
	for(i=3;i<=n;i++)
	{
		k=fscanf(facs_file,"%d %d %d",&factors[i].pr,&factors[i].phi,&factors[i].num_facs);
		//printf("%d: %d %d %d\n",i,factors[i].pr,factors[i].phi,factors[i].num_facs);
		for(j=0;j<factors[i].num_facs;j++)
		{
			k=fscanf(facs_file,"%d %d",&factors[i].primes[j],&factors[i].facs[j]);
			//printf("  %d %d\n",factors[i].primes[j],factors[i].facs[j]);
		}
	}
	return(TRUE);
}

void print_usage()
{
  printf("Usage:- ./upsamhigh <prec> <facs_file>\n");
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



void make_chis(unsigned int q, unsigned int index, unsigned int pr, unsigned int phi, mpfi_c_t *chis, unsigned int no_dims)
{
	unsigned int n,a;
	int even_p;
	//printf("make_chis called with q=%d index=%d pr=%d phi=%d even_p=%d\n",
	//	q,index,pr,phi,even_p);
	if(pr==0)
	{
		make_chis1(q,index,pr,phi,chis);
			return;
	}
	mpfi_c_set(chis[1],c_one);
	even_p=((index&1)==0);
	if (even_p)
		mpfi_c_set(chis[q-1],c_one);
	else
		mpfi_c_set(chis[q-1],c_minus_one);
	if(phi==2)
		return;
	a=pr;
    mpfi_const_pi(L_two_pi_over_phi);
    mpfi_div_ui(L_two_pi_over_phi,L_two_pi_over_phi,phi/2); // phi is always even
		// now I cocked up my initial DFT so that all
		// n-dimensional DFT's <=50 not a power of 2 have sum (exp(2pi i/n) not exp(-2pi i/n)
	    // 1-d DFTs are OK because they just use Bluestein

    if((no_dims==1)||(phi>MAX_SIMPLE_DFT)||(phi==2)||(phi==4)||(phi==8)||(phi==16)||(phi==32))
		mpfi_neg(L_two_pi_over_phi,L_two_pi_over_phi);

	for(n=1;n<phi/2;n++)
	{
	  mpfi_mul_ui(L_theta,L_two_pi_over_phi,((n*index)%phi));
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

inline int conv_index(int index, int q)
{
	return(index);
}

void mpfi_c_L(mpfi_c_ptr res, double re_z, double im_z, int q, int index, factor *factors)
{ 
	int even_p=((index&1)==0);
	int chi_ptr=0;
	int n,a,q_phi=factors[q].phi,this_phi,this_pr,dim,this_q;
	int coords[MAX_DIMS];
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

mpfi_c_t lambda_tmp1;
mpfi_t lambda_tmp2;

void mpfi_c_lambda(mpfi_c_ptr res, double re_z, double im_z, int q, int index, factor *factors, mpfi_c_ptr omega, int neg_one)
{
  mpfi_c_L(res,re_z,im_z,q,index,factors); // res = L(s,chi)
  
  mpfi_c_mul(res,res,omega);               // res = omega*L(s,chi)
  
  if(neg_one)
    mpfi_c_lngamma(lambda_tmp1,(re_z+1.0)/2.0,im_z/2.0);
  else
	  mpfi_c_lngamma(lambda_tmp1,re_z/2.0,im_z/2.0);
  mpfi_mul_d(lambda_tmp2,mpfi_pi_by_4,im_z);
  
  mpfi_c_add_i(lambda_tmp1,lambda_tmp1,lambda_tmp2);
  
  mpfi_c_exp(lambda_tmp1,lambda_tmp1);
  
  mpfi_c_mul(res,res,lambda_tmp1);
  mpfi_set_ui(lambda_tmp2,q);
  mpfi_div(lambda_tmp2,lambda_tmp2,mpfi_pi);
  mpfi_c_set_d(lambda_tmp1,0.0,im_z/2.0);
  mpfi_c_pow_i_c(lambda_tmp1,lambda_tmp2,lambda_tmp1);
  mpfi_c_mul(res,res,lambda_tmp1);
}

void mpfi_c_lambda1(mpfi_c_ptr res, double re_z, double im_z, int q, int index, factor *factors, mpfi_c_ptr omega, int neg_one)
{
 mpfi_c_mul(res,res,omega);               // res = omega*L(s,chi)
  if(neg_one)
    mpfi_c_lngamma(lambda_tmp1,(re_z+1.0)/2.0,im_z/2.0);
  else
	  mpfi_c_lngamma(lambda_tmp1,re_z/2.0,im_z/2.0);
  mpfi_mul_d(lambda_tmp2,mpfi_pi_by_4,im_z);
  
  mpfi_c_add_i(lambda_tmp1,lambda_tmp1,lambda_tmp2);
  mpfi_c_exp(lambda_tmp1,lambda_tmp1);
  mpfi_c_mul(res,res,lambda_tmp1);
  mpfi_set_ui(lambda_tmp2,q);
  mpfi_div(lambda_tmp2,lambda_tmp2,mpfi_pi);
  mpfi_c_set_d(lambda_tmp1,0.0,im_z/2.0);
  mpfi_c_pow_i_c(lambda_tmp1,lambda_tmp2,lambda_tmp1);
  mpfi_c_mul(res,res,lambda_tmp1);
}

int main(int argc, char **argv)
{

  int prec,q,index,i;
  FILE *infile,*facs_file;
  factor *factors;
  q_state qs;

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=3)
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
    
	facs_file=fopen(argv[2],"r");
	if(!facs_file)
	{
	printf("Failed to open file %s for input. Exiting.\n",argv[3]);
	exit(FAILURE);
      }

	factors=(factor *) malloc(sizeof(factor)*MAX_Q);
	if(!factors)
	{
		printf("Fatal error allocating memory for factors. Exiting.\n");
		exit(FAILURE);
	}
	if(!read_factors(factors, MAX_Q, facs_file))
	{
		printf("Fatal error reading facs file. Exiting\n");
		exit(FAILURE);
	}

        mpfi_c_setup(prec);
	mpfi_c_init_L();
	double im_s=5.0/64.0;
	q=87963;
	int ind= 40093;
	int ind_bar=47007;
	int neg_one=1;
	mpfi_c_t res_0,res_1,res_0_bar,res_1_bar,omega,omega_bar;
	mpfi_t x_sqr;
	mpfi_c_init(res_0);mpfi_c_init(res_1);mpfi_c_init(res_0_bar);mpfi_c_init(res_1_bar);
	mpfi_c_init(omega);mpfi_c_init(omega_bar);
	mpfi_init(x_sqr);

	mpfi_c_hur_init(0.5,0.0);
	mpfi_c_L(res_0,0.5,0.0,q,ind,factors);
	mpfi_c_L(res_0_bar,0.5,0.0,q,ind_bar,factors);
	mpfi_sqr(x_sqr,res_0->re);
	mpfi_sqr(omega->re,res_0->im);
	mpfi_add(omega->re,omega->re,x_sqr);
	mpfi_div(omega->re,x_sqr,omega->re);
	mpfi_set_ui(omega->im,1);
	mpfi_sub(omega->im,omega->im,omega->re);
	mpfi_sqrt(omega->re,omega->re);
	mpfi_sqrt(omega->im,omega->im);
	if(mpfi_is_pos(res_0->re))
	  {
	    if(mpfi_is_pos(res_0->im))
	      mpfi_neg(omega->im,omega->im);
	  }
	else
	  {
	  if(mpfi_is_pos(res_0->im))
	    {
	      mpfi_neg(omega->re,omega->re);
	      mpfi_neg(omega->im,omega->im);
	    }
	  else
	    mpfi_neg(omega->re,omega->re);
	  }
	mpfi_c_conj(omega_bar,omega);
	//mpfi_c_print(omega_bar);
	mpfi_c_lambda1(res_0,0.5,0.0,q,ind,factors,omega,neg_one);
	printf("Lambda(0,chi)=");mpfi_c_print(res_0);printf("\n");
	mpfi_c_lambda1(res_0_bar,0.5,0.0,q,ind_bar,factors,omega_bar,neg_one);
	printf("Lambda(0,chi_bar)=");mpfi_c_print(res_0_bar);printf("\n");

        mpfi_c_hur_init(0.5,5.0/64.0);
	mpfi_c_lambda(res_1,0.5,5.0/64.0,q,ind,factors,omega,neg_one);
	printf("Lambda(5/64,chi)=");mpfi_c_print(res_1);printf("\n");
	mpfi_c_lambda(res_1_bar,0.5,5.0/64.0,q,ind_bar,factors,omega_bar,neg_one);
	printf("Lambda(5/64,chi_bar)=");mpfi_c_print(res_1_bar);printf("\n");
	fclose(facs_file);
    printf("Completed.\n");
    return(QUIT_CODE);
}
