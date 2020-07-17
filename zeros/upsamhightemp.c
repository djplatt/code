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
typedef double int_double;
#include "../includes/upsamdefs.h"

#define SUCCESS 0
#define QUIT_CODE 0
//#define FAILURE 1

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
  unsigned int i,j;
  for(i=3;i<=n;i++)
    {
      fscanf(facs_file,"%d %d %d",&factors[i].pr,&factors[i].phi,&factors[i].num_facs);
      //printf("%d: %d %d %d\n",i,factors[i].pr,factors[i].phi,factors[i].num_facs);
      for(j=0;j<factors[i].num_facs;j++)
	{
	  fscanf(facs_file,"%d %d",&factors[i].primes[j],&factors[i].facs[j]);
	  //printf("  %d %d\n",factors[i].primes[j],factors[i].facs[j]);
	}
    }
  return(TRUE);
}

void print_usage()
{
  printf("Usage:- ./upsamhigh <prec> <in_file> <facs_file>\n");
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


inline int ln2(int q)
{
  int res=0;
  while(q)
    {
      q>>=1;
      res++;
    }
  return(res);
}

inline int power_2_p(unsigned int i)
{
  while((i&1)==0)
    i>>=1;
  return(i==1);
}


// q is 2^a, a>=3, phi = q/4
void make_chis1(unsigned int q, unsigned int index, unsigned int pr, unsigned int phi, mpfi_c_t *chis)
{
  unsigned int pow,chi_ptr,n,phi_by_2=phi>>1;
  int even_p=(index<phi_by_2);
  unsigned int w=0;
  mpfi_const_pi(L_two_pi_over_phi);
  mpfi_div_ui(L_two_pi_over_phi,L_two_pi_over_phi,phi/4);
  //if((ln2(q)&1)==1)
  if(q==16) mpfi_neg(L_two_pi_over_phi,L_two_pi_over_phi);
  chi_ptr=5;
  mpfi_c_set_ui(chis[1],1,0);
  if(even_p)
    mpfi_c_set_ui(chis[q-1],1,0);
  else
    mpfi_c_set_d(chis[q-1],-1.0,0.0);

  for(pow=1;pow<q/4;pow++)
    {
      w+=index;
      while(w>=phi_by_2)
	w-=phi_by_2;
      mpfi_mul_ui(L_theta,L_two_pi_over_phi,w);
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
  unsigned int n,a,w=0;
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
  mpfi_neg(L_two_pi_over_phi,L_two_pi_over_phi);
  // now I cocked up my initial DFT so that all
  // n-dimensional DFT's <=50 not a power of 2 have sum (exp(2pi i/n) not exp(-2pi i/n)
  // 1-d DFTs are OK because they just use Bluestein

  if(((no_dims==1)&&(phi<=MAX_SIMPLE_DFT))||((no_dims>1)&&(phi>4)&&((phi<=MAX_SIMPLE_DFT)||(power_2_p(phi)))))
    mpfi_neg(L_two_pi_over_phi,L_two_pi_over_phi);

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
  //printf("coords=[");
  //for(dim=0;dim<factors[q].num_facs;dim++)
  //printf(" %d",coords[dim]);
  //printf(" ]\n");

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
      if ((n&4095)==0) printf("%d\n",n);
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
  //printf("Gamma()*exp()=");mpfi_c_print(lambda_tmp1);
  mpfi_c_mul(res,res,lambda_tmp1);
  mpfi_set_ui(lambda_tmp2,q);
  mpfi_div(lambda_tmp2,lambda_tmp2,mpfi_pi);
  mpfi_c_set_d(lambda_tmp1,0.0,im_z/2.0);
  mpfi_c_pow_i_c(lambda_tmp1,lambda_tmp2,lambda_tmp1);
  //printf("(q/pi)^it/2=");mpfi_c_print(lambda_tmp1);
  mpfi_c_mul(res,res,lambda_tmp1);
}

int my_fread(q_state *qs,int a, int b, FILE *infile)
{
  char buff[16];
  if(!fread(&qs->q,sizeof(int),1,infile))
    return(FALSE);
  fread(&qs->index,sizeof(int),1,infile);
  fread(&qs->rate,sizeof(int),1,infile);
  fread(buff,sizeof(char),4,infile);
  fread(&qs->gap,sizeof(double),1,infile);
  fread(&qs->neg_one,sizeof(char),1,infile);
  fread(&qs->type,sizeof(char),1,infile);
  fread(buff,sizeof(char),2,infile);
  fread(&qs->n0,sizeof(int),1,infile);
  fread(&qs->offset1,sizeof(int),1,infile);
  fread(&qs->offset2,sizeof(int),1,infile);
  fread(buff,sizeof(char),8,infile);
  fread(&qs->omega,sizeof(double),4,infile);
  qs->omega[1]=-qs->omega[1];
  qs->omega[3]=-qs->omega[3];
  fread(buff,sizeof(char),16,infile);
  return(TRUE);
}


inline void make_coords(int *coords, int q, int index, factor *factors)
{
  int n=index,dim;
  for(dim=factors[q].num_facs-1;dim>=0;dim--)
    {
      coords[dim]=n%factors[factors[q].facs[dim]].phi;
      n/=factors[factors[q].facs[dim]].phi;
    }
  return;
}


inline int conj_index(int index, unsigned int q, factor *factors)
{
  int coords[MAX_DIMS],res,i;

  make_coords(coords,q,index,factors);
  for(i=0;i<factors[q].num_facs;i++)
    coords[i]=factors[factors[q].facs[i]].phi-coords[i];
  res=coords[0];
  for(i=1;i<factors[q].num_facs;i++)
    res=res*factors[factors[q].facs[i]].phi+coords[i];
  return(res);
}



int main(int argc, char **argv)
{

  int prec,q,index,i;
  FILE *infile,*facs_file;
  factor *factors;
  q_state qs;
  double n0,offset;


  /*  check all the command line arguments are ok, if not print message
      and exit sharpish */


  if(argc!=6)
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
    
    

  infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Failed to open file %s for input. Exiting.\n",argv[2]);
      exit(FAILURE);
    }

  facs_file=fopen(argv[3],"r");
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
  double im_s;
  mpfi_c_t res,omega;
  mpfi_c_init(res);
  mpfi_c_init(omega);
	

  while(fread(&qs,sizeof(q_state),1,infile))
    {
      printf("index=%d\n",qs.index);
      if(qs.index<0)
	qs.index=conj_index(-qs.index,qs.q,factors); // shouldn't happen
      //
      // hardwired bit
      //
      n0=atof(argv[4]);
      offset=atof(argv[5]);
      im_s=(n0+offset/(double) qs.rate)*qs.gap;
      //im_s=(131971.0+462.75/(double) qs.rate)*qs.gap;
      printf("t=%30.28e\n",im_s);
      //
      // end hardwired bit
      //
      //      if(qs.type==CHECK_PROB)
      //	im_s=qs.n0*qs.gap+(qs.offset1+qs.offset2)/2.0*qs.gap/(double) qs.rate;
      //      else
      //	im_s=qs.n0*qs.gap+qs.offset1*qs.gap/(double) qs.rate;
      mpfi_c_hur_init(0.5,im_s);
      //for(i=0;i<4;i++)
      //printf("omega[%d]=%20.18e\n",i,qs.omega[i]);
      mpfi_interv_d(omega->re,qs.omega[0],qs.omega[1]);
      mpfi_interv_d(omega->im,qs.omega[2],qs.omega[3]);
      //printf("omega=");mpfi_c_print(omega);
      //printf("q=%d\nindex=%d\nim_s=%30.28e\n",qs.q,qs.index,im_s);
      //if(qs.neg_one)
      //printf("This character is odd.\n");
      //else
      //printf("This character is even.\n");
      mpfi_c_lambda(res,0.5, im_s, qs.q, qs.index, factors, omega, qs.neg_one);
      printf("q=%d index=%d t=%30.28e Lambda returned ",qs.q,qs.index,im_s);mpfi_c_print(res);
      if(qs.exp_sign==POS)
	printf("q=%d index=%d Wanted a +ve interval.\n",qs.q,qs.index);
      else
	printf("q=%d index=%d Wanted a -ve interval.\n",qs.q,qs.index);
    }
  fclose(infile);
  fclose(facs_file);
  printf("Completed.\n");
  return(QUIT_CODE);
}
