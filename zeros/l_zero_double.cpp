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
#include "../includes/im_s.h"
#include "../includes/fft_defs.h"
#include "../includes/int-fft.h"
#include "../includes/upsamdefs.h"
#define DEBUG printf("Reached line %d\n",__LINE__)
// q is 2^a, a>=3, phi = q/4

int_complex chis[MAX_Q];

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

inline bool power_2_p(unsigned int i)
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
  int_double four_pi_over_phi=d_pi/(phi/4),theta;
  unsigned int w=0;
  //printf("q=%d\nln2(q)=%d\n",q,ln2(q));

  //if((ln2(q)&1)==1)
  //   four_pi_over_phi=-four_pi_over_phi;
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
      sin_cos(theta,&chis[chi_ptr].imag,&chis[chi_ptr].real);
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
  unsigned int n,a,chi_ptr,w=0;
  int even_p;
  int_double two_pi_over_phi,theta;
  //printf("%d\n",sizeof(l));
  //exit(0);
  //printf("make_chis called with q=%d index=%d pr=%d phi=%d\n",
  //q,index,pr,phi);
  if(pr==0)
    {
      make_chis1(q,index,pr,phi,chis);
      return;
    }
  chis[1]=c_one;
  even_p=((index&1)==0);
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

inline void make_coords(int *coords, int q, int index, factor *factors)
{
  int n=index,dim;
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


int_complex L(double re_z, double im_z, int q, int index, factor *factors, int_complex *s_array)
{ 
  int even_p=((index&1)==0);
  int chi_ptr=0;
  int n,a,q_phi=factors[q].phi,this_phi,this_pr,dim,this_q;
  int coords[MAX_DIMS];
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
  res=c_zero;
	
  for(n=1;n<q;n++)
    {
      //if(n==35) exit(0);
      if(co_prime(n,q))
	{	
	  hur=hurwitz1(re_z,im_z,int_double(n)/q,s_array);
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

inline int_complex lambda(double re_z, double im_z, int q, int index, factor *factors, const int_complex &omega, int neg_one, int_complex *s_array)
{
  int_complex res=L(re_z,im_z,q,index,factors,s_array); // res = L(s,chi)
  res*=omega;               // res = omega*L(s,chi)
  if(neg_one)
    res*=exp(lngamma(int_complex(int_double((re_z+1.0)/2.0),int_double(im_z/2.0)))+d_pi*(im_z/4.0));
  else
    res*=exp(lngamma(int_complex(int_double(re_z/2.0),int_double(im_z/2.0)))+d_pi*(im_z/4.0));

  res*=pow(int_double(q)/d_pi,int_complex(int_double(0.0),int_double(im_z/2.0)));
  return(res);
}

int my_fread(q_state *qs,int a, int b, FILE *infile)
{
  char buff[16];
  if(!fread(&qs->q,sizeof(int),1,infile))
    return(false);
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
  fread(buff,sizeof(char),16,infile);
  return(true);
}

inline int conj_index(int index, unsigned int q, factor *factors)
{
  int coords[MAX_DIMS],res;

  make_coords(coords,q,index,factors);
  
  if((q&7)==0) // multiple of 8
    {
      res=q>>2;
      if(index<res) // top row
	coords[0]=res-coords[0];
      else // 2nd row
	coords[0]=3*res-coords[0];
    }
  else
    coords[0]=factors[factors[q].facs[0]].phi-coords[0];

  for(int i=1;i<factors[q].num_facs;i++)
    coords[i]=factors[factors[q].facs[i]].phi-coords[i];
  res=coords[0];
  for(int i=1;i<factors[q].num_facs;i++)
    res=res*factors[factors[q].facs[i]].phi+coords[i];
  return(res);
}

int main(int argc, char **argv)
{

  int prec,q,index,i;
  ifstream facs_file;
  factor *factors;
  int_complex s_array[MAX_BERNOULLI_K];

  _fpu_rndd();
  set_h_bernoulli();


  /*  check all the command line arguments are ok, if not print message
      and exit sharpish */


  if(argc!=4)
    {
      printf("usage:- l_zero_double <q> <index> <facs_file>\nExiting.\n");
      exit(1);
    }
  q=atoi(argv[1]);
  if((q<=2)||(q>MAX_Q))
    {
      printf("Bad q=%d. Exiting.\n",q);
      exit(1);
    }

  facs_file.open(argv[3]);
  if(!facs_file.is_open())
    {
      printf("Failed to open file %s for input. Exiting.\n",argv[3]);
      exit(1);
    }

  factors=(factor *) malloc(sizeof(factor)*MAX_Q);
  if(!factors)
    {
      printf("Fatal error allocating memory for factors. Exiting.\n");
      exit(1);
    }
  if(!read_factors(factors, MAX_Q, facs_file))
    {
      printf("Fatal error reading facs file. Exiting\n");
      exit(1);
    }

  index=atoi(argv[2]);
  if((index<0)||(index>=factors[q].phi))
    {
      printf("Bad index=%d. Exiting.\n",index);
      exit(1);
    }
  double im_s;
  int_complex res,omega;

   hur_init(s_array,0.5,0);
   res=L(0.5,0,q,index,factors,s_array);
   if(contains_zero(res))
     {
       printf("L(1/2,chi_%d_%d) contains zero. Try at higher precison.\n",q,index);
       print_int_complex_str("result was ",res);
     }
   else
     {
       printf("L(1/2,chi_%d_%d) non-zero.\n",q,index);
       print_int_complex_str("result was ",res);
     }
  facs_file.close();
  printf("Completed.\n");
  return(0);
}

