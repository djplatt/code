/*
File: int-l-func5.5.cpp

Created: 10th Jan 2009

Version: <v> = 5.5

Last Modified: 16th September 2009

5.0 changed summation code in calc_rn and calc_rn_neg
5.1 now multiply by (q/pi)^(it/2)
5.2 moved (q/pi)^(it/2) and q^(-s) into FFT
5.3 changed calc of q_comb for accuracy, introduced im_s_vec[].pi_minus_it_2
5.4 using semaphores to access lattice files
5.5 read lattice file in one go using read
5.6 improved calculation of hurwitz values
5.7 now uses crlibm 

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

Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "5.7"


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <semaphore.h>
#include <assert.h>
using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft-half.h"
#include "../includes/make-factors.h"

void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: int-l-func%s (q-start) (q-end) (ifname) (ofname) (N) (semaphore name)\n",VERSION);
  printf("  (q-start)   - integer >=3\n");
  printf("  (q-end)     - integer >=(q-start)\n");
  printf("  (ifname)    - file with lattice values.\n");
  printf("  (ofname)    - output file.\n");
  printf("  (N)         - number of Taylor terms to use.\n");
  printf("  (sem name)  - name of semaphore to control file access.\n");
  exit(1);
}


void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};


	    
// only a vector here. last element used for taylor error
void create_s_array(const int_complex &s,int_complex *s_array, unsigned int n, const double gap)
{
  unsigned int i;
  s_array[0].real=int_double(gap/2.0); // 0.5*gap is representable
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

//unsigned int temp_count=0;

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
  int_double q_minus_half;
  int dims[MAX_FACS];
  int no_dims;
  bool power_2,primitive,neg_one;
  unsigned int pr=factors[q].pr,n_prims;
  unsigned int i,j,offset,s_done;
  //long double lnq=ln_q(q);

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
    //printf("pr of %ld in %ld\n",q,pr);
    fill_q_ns(q_n,pr,phi_q,q);
    init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
    prep_omegas(omegas,q,phi_q,q_n,a,bs,b_star_conjs,conv_sizes[0],b_spare);
    //for(i=0;i<phi_q;i++)
    //  print_int_complex_str("prepped omegas[]= ",omegas[i]);
    s_done=0;
    my_hur_ptr1=hur_ptr;
    my_hur_ptr=hur_ptr;
    //printf("my_hur_ptr %ld\n",my_hur_ptr);
    while(s_done<num_s)
    {
      q_comb=im_s_vec[s_done].pi_minus_it_2*q_minus_half;//*conj(q_it(q,im_s_vec[s_done].im_s/2.0));
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
      // now fft_vec contains L(chi,s)*q^(s)
      /*
      for(i=0;i<phi_q;i++)
	{
	  print_int_complex_str("FFT produced ",Fft_vec[i]*q_minus_half/q_it(q,im_s_vec[s_done].im_s));
	  print_int_complex_str("L   produced ",L(int_complex(int_double(0.5),int_double(im_s_vec[s_done].im_s)),q,i,factors));
	}
      exit(0);
      */
      for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
        if(prim_p(i,factors[q].primes[0]))        
        {
          neg_one=(i&1);

          if(s_done==0) // first time through, so finish off omega values
          {
            omegas[i]=finish_omega(omegas[i],neg_one);
	    //print_int_complex_str("finished omega[]= ",omegas[i]);
            fwrite(&i,sizeof(unsigned int),1,out_file);
            fwrite(&omegas[i],sizeof(int_complex),1,out_file);
            fwrite(&neg_one,sizeof(bool),1,out_file);
          }
          z=Fft_vec[i];
	  z*=q_comb;
          if(neg_one) // chi(-1)=-1
            z*=omegas[i]*im_s_vec[s_done].lambda_s_a;
          else
            z*=omegas[i]*im_s_vec[s_done].lambda_s;

          if(!contains_zero(z.imag))
          {
            printf("Imaginary part of z does not contain zero.\n");
	    print_int_complex_str("Fft_vec[]= ",Fft_vec[i]);
	    print_int_complex_str("z= ",z);
	    printf("im_s=%20.18e\n",im_s_vec[s_done].im_s);
	    printf("q=%d\n",q);
	    print_int_complex_str("omega ",omegas[i]);
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
	  //print_int_double_str("z=",z.real);
	  if(contains_zero(z.real))
	    printf("L(1/2,chi) may contain zero for q=%ld index=%ld\n",q,i);

          fwrite(&z.real,sizeof(int_double),1,out_file);
        }
        s_done++;
	my_hur_ptr1+=num_fracs;
	my_hur_ptr=my_hur_ptr1;
    }
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
    prep_omegas_nd(offset_vec,q,omegas,no_dims,dims,phi_q);
    //printf("omegas prepped.\n");
    my_hur_ptr1=hur_ptr;
    my_hur_ptr=hur_ptr;
//    printf("my_hur_ptr %ld\n",my_hur_ptr);
    while(s_done<num_s)
    {
      //sin_cos(-log(int_double(q))*(double) (im_s_vec[s_done].im_s/2.0),&q_minus_it_2.imag,&q_minus_it_2.real);
      //q_comb=q_minus_it(q,im_s_vec[s_done].im_s/2)*im_s_vec[s_done].pi_minus_it_2*pow(int_double(q),-0.5);
      //print_int_complex_str("q^it=",q_it(q,im_s_vec[s_done].im_s/2.0));
      q_comb=im_s_vec[s_done].pi_minus_it_2*q_minus_half;//*conj(q_it(q,im_s_vec[s_done].im_s/2.0));
      for(i=0;i<phi_q;i++)
	Fft_vec[offset_vec[i]]=hur_vec[my_hur_ptr++];
      //printf("Doing FFT.\n");
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
		  //print_int_complex_str("q^(it)=",q_it(q,im_s_vec[s_done].im_s));
		  //print_int_complex_str("q^(it/2)=",q_it(q,im_s_vec[s_done].im_s/2.0));
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
	      //print_int_double_str("z=",z.real);
	      if(contains_zero(z.real))
		printf("L(1/2,chi) may contain zero for q=%ld index=%ld\n",q,i);
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
		  unsigned int q_end, FILE *out_file,int_complex &taylor_err)
{
  int_complex out_val;
  unsigned int i,j,q,num,steps,gap_ptr;
  int_double x,frac;
  int_complex q_s,first_term;
  double dsteps,csteps,dfloor;
  long double lnq;

  //printf("in calc_hurwitz1\n");
  //printf("setup completed.\n");
  

  num_fracs=0;
  __builtin_prefetch(s_array,0,3);
  for(q=q_start;q<=q_end;q++)
    {
      if((q&3)==2)
	continue;      // if q = 2 mod 4 then no primitive characters
      q_s=sqrt(int_double(q));//q_it(q,s.imag.left)*sqrt(int_double(q));//pow(int_double(q),s.real.left);

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
	    out_val=taylor_err+q_s/sqrt(int_double(num));//conj(q_it(num,s.imag.left))*q_s/sqrt(int_double(num));//pow(int_double(num),-0.5)*
	    out_val+=calc_r_n(r_n_terms,s_array,x,&r_n_vals[gap_ptr]);
	    // out_val now contains Zeta(s,num/q)
	    hur_vec[hur_ptr++]=out_val;
	  }
    }
  return(num_fracs);
}

int main(int argc, char **argv)
{

  _fpu_rndd();

  int_complex s,gam,*s_array,*r_n_vals,*omegas;
  im_s im_s_vec[1];
  unsigned int *a_n,*q_n,*offset_vec,num_fracs,q,max_num_prims,n_prim,i;
  int no_gaps,q_start,q_end,num_s,N,rn,file_N;
  int_double mgap,im_s_2_mod;
  double gap;
  FILE *in_file,*out_file;
  factor *factors;

  if(argc!=7)
    print_usage();
  q_start=atoi(argv[1]);
  if(q_start<3)
    print_usage();
  q_end=atoi(argv[2]);
  if(q_end<q_start)
    print_usage();
  in_file=fopen(argv[3],"rb");
  if(!in_file)
    {
      perror("Error opening in_file. ");
      printf("Failed to open file |%s|\n",argv[3]);
      fatal_error("Couldn't open in_file. Exiting.\n");
    }

  if(!(factors=(factor *) calloc(q_end+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  //printf("Making factors.\n");
  if(!make_factors(factors,q_end))
    fatal_error("Error creating factors. Exiting.\n");

  N=atoi(argv[5]);
  /*
  sem_t *lfuncsem=sem_open(argv[6],0);
  if(lfuncsem==SEM_FAILED)
    {
      perror("");
      fatal_error("Can't open semaphore. Exiting.\n"); 
    }
  */
  fread(&num_s,sizeof(int),1,in_file);
  assert(num_s==1);
  fread(&file_N,sizeof(int),1,in_file);
  if((N>file_N)||(N<=0))
    fatal_error("N<=0 or exceeds N used to save lattice file. Exiting.\n");
  fread(&rn,sizeof(int),1,in_file);
  if(rn!=1)
    fatal_error("RN must be 1. Exiting.\n");
  out_file=fopen(argv[4],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");
  fread(&no_gaps,sizeof(int),1,in_file);

  gap=1.0/no_gaps;
  //printf("gap set to %10.8e\n",gap);


  //cout << "Allocating memory for s_array." << endl;
  if(!(s_array=(int_complex *) _aligned_malloc(sizeof(int_complex)*N,16)))
    fatal_error("Can't allocate memory for s_array. Exiting.");
  
  
  //cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=(int_complex *) _aligned_malloc(sizeof(int_complex)*(no_gaps+1)*N,16)))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");


  //fwrite(&q_start,sizeof(int),1,out_file);
  //fwrite(&q_end,sizeof(int),1,out_file);


  //fwrite(&max_num_prims,sizeof(int),1,out_file);
  //fwrite(&num_s,sizeof(int),1,out_file);

  s.real=int_double(0.5);

  if(!(hur_vec=(int_complex *) _aligned_malloc(q_end*sizeof(int_complex)*num_s,16)))
    fatal_error("Failed to allocate memory for hur_vec. Exiting.");

  // try and force everything resident
  /*
  for(i=0;i<num_fracs*num_s;i+=getpagesize()/sizeof(int_complex))
    hur_vec[i]=c_zero;
  */
  if(!(lattice_buff=(double *) malloc((13+4*(no_gaps+1)*file_N)*sizeof(double))))
    fatal_error("Failed to allocate memory for lattice_buff. Exiting.\n");
  printf("allocated memory\n");
  
  //sem_wait(lfuncsem);
  
  //printf("remember to uncomment sem_wait etc.\n");
  fread(lattice_buff,sizeof(double),(13+4*(no_gaps+1)*file_N),in_file);
  printf("read lattice file\n");
  fclose(in_file);
  
  //sem_post(lfuncsem);
  
  im_s_vec[i].im_s=read_double();
  //printf("im_s=%10.8e\n",im_s_vec[i].im_s);
  s.imag=int_double(im_s_vec[i].im_s);
  im_s_vec[i].lambda_s=read_int_complex(); // GAMMA(s/2)*exp(Pi t/4)
  im_s_vec[i].lambda_s_a=read_int_complex(); // GAMMA((s+1)/2)*exp(Pi t/4)
  im_s_vec[i].pi_minus_it_2=read_int_complex();
  //print_int_complex_str("gamma(s/2)*exp(Pi t/4)=",im_s_vec[i].lambda_s);
  //exit(0);

  if(!(omegas=(int_complex *) _aligned_malloc((q_end-1)*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  //printf("initialising ws\n");
  init_ws();
  //printf("ws initialised\n");

  create_s_array(s,s_array,N,gap); // end column used for error estimation
  int_complex taylor_err=taylor_error(s,N,gap,s_array);
  int j,gap_ptr=0;
  for(i=0;i<=no_gaps;i++)
  {

    for(j=0;j<N;j++)
      r_n_vals[gap_ptr++]=read_int_complex();
    skip_int_complex(file_N-N);
  }


  for(q=q_start;q<=q_end;q++)
    if((q&3)!=2)
      {
	num_fracs=calc_hurwitz1(s,N,file_N,gap,s_array,r_n_vals,no_gaps,q,q,out_file,taylor_err);
	//printf("calc_hurwitz1 completed.\n");
	//exit(0);
	if(!num_fracs)
	  fatal_error("Error running Hurwitz routine. Exiting.");
	hur_ptr=0;
	      
	//printf("Processing q=%ld\n",q);
	make_l(q,num_s,factors,q_n,im_s_vec,offset_vec,a_n,omegas,out_file);
      }

  fclose(out_file);

  //printf("FFT Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks1)/((double) CLOCKS_PER_SEC));
  //printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  //printf("Lambda value straddled zero %d times\n",temp_count);
  printf("l-func successful completion on q=%d-%d file %s\n",q_start,q_end,argv[3]);
  return(0);
}
