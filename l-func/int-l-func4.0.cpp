/*
File: int-l-func4.0.cpp

Created: 10th Jan 2009

Version: <v> = 4.0

Last Modified: 13th May 2009

Dialect: C++

Requires: No Libraries required

Implementation notes: Takes Output from l-func-mpfi-1.0
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

#define VERSION "4.0"



#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
using namespace std;



#include "../includes/int_double10.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft.h"


void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: int-l-func%s (q-start) (q-end) (ifname) (ofname) (N)\n",VERSION);
  printf("  (q-start)   - integer >=3\n");
  printf("  (q-end)     - integer >=(q-start)\n");
  printf("  (ifname)    - file with lattice values.\n");
  printf("  (facs_file) - file with factors data.\n");
  printf("  (ofname)    - output file.\n");
  printf("  (N)         - number of Taylor terms to use.\n");
  exit(0);
}


void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(0);
};

// only a vector here. last element used for taylor error
void create_s_array(int_complex &s,int_complex *s_array, unsigned int n, double &gap)
{
  unsigned int i;
  s_array[0].real=int_double(gap/2); // 0.5*gap is representable
  s_array[0].imag=s.imag*gap;
  for(i=1;i<n;i++)
    s_array[i]=s_array[i-1]*(s+i)*gap/(i+1); // ditto and (i+1)
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

/*
// sum (alpha+n)^(-s) from n=R to n=RN
int_complex zeta (int_complex &s, int_double &alpha, unsigned int start, unsigned int end)
{
  int_complex tot=c_zero;
  int i,en=end,st=start;
  for(i=en;i>=st;i--)
    tot=tot+pow(alpha+(double)i,-s); // do me accurately for large Im(s)
  return(tot);
}
*/

int_complex calc_r_n(unsigned int r_n_terms, int_complex *s_array, int_double &delta,
      int_complex *r_n_vec,int_double &alpha)
{
  unsigned int i;
  int_complex res,term;
  int_double delta_n=delta;

  res=r_n_vec[0];


  for(i=0;;)
    {
      term=s_array[i]*delta_n;
    i++;
      term=term*r_n_vec[i];
      res=res+term;
    if(i==r_n_terms-1)
      break;
    else
      delta_n=times_pos(delta_n,delta);
    }
  return(res);
}

int_complex calc_r_n_neg(unsigned int r_n_terms, int_complex *s_array, int_double &delta,
      int_complex *r_n_vec,int_double &alpha)
{
  unsigned int i;
//  int_double ln_delta=log(-delta);
  int_double delta_n=delta; // will be negative
//  int sign=-1;
  int_complex res,term;

  res=r_n_vec[0];


  for(i=0;i<r_n_terms-1;)
    {
      term=s_array[i]*delta_n;
    delta_n=delta_n*delta; // can't use quick multiply
      i++;
      term=term*r_n_vec[i];
    res=res+term;
    }
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

int_double read_int_double(FILE *in_file)
{
  double x[2];
  int_double ans;
  fread(&x,sizeof(double),2,in_file);
  ans.left=x[0];
  ans.right=-x[1];
  return(ans);
}


int_complex read_int_complex(FILE *in_file)
{
  double x[4];
  int_complex ans;
  fread(x,sizeof(double),4,in_file);
  ans.real.left=x[0];
  ans.real.right=-x[1];
  ans.imag.left=x[2];
  ans.imag.right=-x[3];
  return(ans);
}

void skip_int_complex(FILE *in_file,unsigned int n)
{
  unsigned int i;
  double x[4];
  for(i=0;i<n;i++)
    fread(x,sizeof(double),4,in_file);
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
  int_complex omega,omega_a,z,z1;
  int dims[MAX_FACS];
  int no_dims;
  bool power_2,primitive,neg_one;
  unsigned int pr=factors[q].pr,n_prims;
  unsigned int i,j,offset,s_done;

//  printf("Processing Q=%d\n",q);

  fwrite(&q,sizeof(unsigned int),1,out_file);
  n_prims=num_prims(q,factors); // no of prim characters mod q
  fwrite(&n_prims,sizeof(unsigned int),1,out_file);

  j=0;
  for(i=1;i<q;i++)
    if(co_prime(i,q))
      a_n[j++]=i;
  if(pr)  // q has a primitive root so nice and easy 1-d FFT
  {
    fill_q_ns(q_n,pr,phi_q,q);
    init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
    prep_omegas(omegas,q,phi_q,q_n,a,bs,b_star_conjs,conv_sizes[0],b_spare);
    s_done=0;
    my_hur_ptr1=hur_ptr;
    my_hur_ptr=hur_ptr;
//    printf("my_hur_ptr %ld\n",my_hur_ptr);
    while(s_done<num_s)
    {
      for(i=1;i<q;i++)
        if(co_prime(i,q))
        {
		Fft_vec[q_n[i]]=hur_vec[my_hur_ptr++];
//          fread(&fft_vec[q_n[i]],sizeof(int_complex),1,hur_file);
        }
      bluestein_fft(1,Fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare);
      
      // now fft_vec contains L(chi,s)*q^(s/2)
      for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
        if(prim_p(i,factors[q].primes[0]))        
        {
          neg_one=(i&1);

          if(s_done==0) // first time through, so finish off omega values
          {
            omegas[i]=finish_omega(omegas[i],neg_one);
            fwrite(&i,sizeof(unsigned int),1,out_file);
            fwrite(&omegas[i],sizeof(int_complex),1,out_file);
            fwrite(&neg_one,sizeof(bool),1,out_file);
          }
          z=Fft_vec[i];      
          if(neg_one) // chi(-1)=-1
            z=z*omegas[i]*im_s_vec[s_done].lambda_s_a;
          else
            z=z*omegas[i]*im_s_vec[s_done].lambda_s;

          if(!contains_zero(z.imag))
          {
            printf("Imaginary part of z does not contain zero.\n");
            exit(0);
          }

          fwrite(&z.real,sizeof(int_double),1,out_file);
        }
        s_done++;
	my_hur_ptr1+=num_fracs;
	my_hur_ptr=my_hur_ptr1;
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

      for(i=0;i<phi_q;i++)
	      Fft_vec[offset_vec[i]]=hur_vec[my_hur_ptr++];
//        fread(&fft_vec[offset_vec[i]],sizeof(int_complex),1,hur_file);
/*
      if((q==12)&&(im_s_vec[s_done].im_s==150.46875))
      {
	      int_complex s=int_complex(int_double(0.5),int_double(im_s_vec[s_done].im_s));
	      for(i=0;i<4;i++)
		      print_int_complex_str("fft_vec[]=",fft_vec[offset_vec[i]]);
	      print_int_complex_str("hur(s,1/12)=",hurwitz(s,int_double(1)/12)*pow(12,-s/2));
	      print_int_complex_str("hur(s,5/12)=",hurwitz(s,int_double(5)/12)*pow(12,-s/2));
	      print_int_complex_str("hur(s,7/12)=",hurwitz(s,int_double(7)/12)*pow(12,-s/2));
	      print_int_complex_str("hur(s,11/12)=",hurwitz(s,int_double(11)/12)*pow(12,-s/2));
      }
*/      
      do_nd_fft(Fft_vec,no_dims,dims,phi_q);

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
              omegas[i]=finish_omega(omegas[i],neg_one);
              fwrite(&i,sizeof(unsigned int),1,out_file);
              fwrite(&omegas[i],sizeof(int_complex),1,out_file);
              fwrite(&neg_one,sizeof(bool),1,out_file);
            }

            if(neg_one)
              z*=omegas[i]*im_s_vec[s_done].lambda_s_a;
            else
              z*=omegas[i]*im_s_vec[s_done].lambda_s;

            if(!contains_zero(z.imag))
            {
              printf("Non Zero Imag Part %d %d %d\n",q,i,j);
			  printf("Im(s)=%10.8e\n",im_s_vec[s_done].im_s);
              exit(0);
            }

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
          unsigned int q_end, FILE *in_file, FILE *out_file)
{
  int_complex out_val;
  int_complex taylor_err;
  unsigned int i,j,gap_ptr,q,num,steps;
  int_double x,frac;
  int_complex q_minus_s_2,minus_s;
  double dsteps,csteps;

  minus_s=-s;

  create_s_array(s,s_array,r_n_terms,gap); // end column used for error estimation
  taylor_err=taylor_error(s,r_n_terms,gap,s_array);
  gap_ptr=0;
  
  for(i=0;i<=no_gaps;i++)
  {

    for(j=0;j<r_n_terms;j++)
      r_n_vals[gap_ptr++]=read_int_complex(in_file);
    skip_int_complex(in_file,file_N-r_n_terms);
  }
  num_fracs=0;
  for(q=q_start;q<=q_end;q++)
    {
    if((q&3)==2)
      continue;      // if q = 2 mod 4 then no primitive characters
    q_minus_s_2=pow1(q,-0.25,-s.imag.left/2.0); // q,Re(s),Im(s) all representable
    for(num=1;num<q;num++)
      if(co_prime(num,q))
        {
        num_fracs++;
        frac=int_double(num)/q;
        dsteps=double(num)/(gap*q);
        csteps=ceil(dsteps);
        steps=(int) csteps;
        if((csteps-dsteps)>=0.5) // use lower row of R_ns
        {
          gap_ptr=r_n_terms*(no_gaps-steps);
          x=int_double(steps)-frac/gap;
          out_val=calc_r_n(r_n_terms,s_array,x,&r_n_vals[gap_ptr],frac);
        }
        else // use upper row
        {
          gap_ptr=r_n_terms*(no_gaps-steps+1);
          x=int_double(steps)-(frac/gap+1.0);
          out_val=calc_r_n_neg(r_n_terms,s_array,x,&r_n_vals[gap_ptr],frac);
        }
        // out_val is now ~ sum_{n=RN}^\infty (n+num/q)^(-s)
        
	out_val=out_val+taylor_err;
        out_val=out_val+pow1(frac,-s.real.left,-s.imag.left);
        out_val=out_val*q_minus_s_2; // out_val now contains Zeta(s,num/q)*q^(-s/2)
		//print_int_complex_str("out_val ",out_val);
	
	hur_vec[hur_ptr++]=out_val;
//        fwrite(&out_val,sizeof(int_complex),1,out_file);
      }
  }
  return(num_fracs);
}

int main(int argc, char **argv)
{

  _fpu_rndd();

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
    fatal_error("Couldn't open in_file. Exiting.\n");

  facs_file.open(argv[4]);
    if(!facs_file.is_open())
      fatal_error("Couldnt open factors file. Exiting.\n");

  if(!(factors=(factor *) calloc(q_end+1,sizeof(factor))))
    fatal_error("Error allocating memory for factors. Exiting.\n");

  if(!read_factors(factors,q_end,facs_file))
    fatal_error("Error reading factor file. Exiting.\n");

  facs_file.close();

  N=atoi(argv[6]);
  fread(&num_s,sizeof(int),1,in_file);
  fread(&file_N,sizeof(int),1,in_file);
  if((N>file_N)||(N<=0))
    fatal_error("N<=0 or exceeds N used to save lattice file. Exiting.\n");
  fread(&rn,sizeof(int),1,in_file);
  if(rn!=1)
    fatal_error("RN must be 1. Exiting.\n");
  out_file=fopen(argv[5],"wb");
  if(!out_file)
    fatal_error("Couldn't open out_file. Exiting.\n");
  fread(&no_gaps,sizeof(int),1,in_file);
  gap=1.0/no_gaps;


//  cout << "Allocating memory for s_array." << endl;
  if(!(s_array=(int_complex *) _aligned_malloc(sizeof(int_complex)*N,16)))
    fatal_error("Can't allocate memory for s_array. Exiting.");
  
  
//  cout << "Allocating memory for R_N vectors." << endl;
  if(!(r_n_vals=(int_complex *) _aligned_malloc(sizeof(int_complex)*(no_gaps+1)*N,16)))
    fatal_error("Can't allocate memory for r_n_vals. Exiting.");

//  cout << "Allocating memory for im_s_vec." << endl;
  if(!(im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*num_s,16)))
    fatal_error("Can't allocate memory for im_s_vec. Exiting.");


  fwrite(&q_start,sizeof(int),1,out_file);
  fwrite(&q_end,sizeof(int),1,out_file);

  max_num_prims=0;
  num_fracs=0;
  for(q=q_start;q<=q_end;q++)
  {
    if((q&3)!=2)
    {
      num_fracs+=factors[q].phi;
      n_prim=num_prims(q,factors);
      if(n_prim>max_num_prims)
        max_num_prims=n_prim;
    }
  }

  fwrite(&max_num_prims,sizeof(int),1,out_file);
  fwrite(&num_s,sizeof(int),1,out_file);

  s.real=int_double(0.5);

  if(!(hur_vec=(int_complex *) _aligned_malloc(num_fracs*sizeof(int_complex)*num_s,16)))
    fatal_error("Failed to allocate memory for hur_vec. Exiting.");

  for(i=0;i<num_s;i++)
    {
    fread(&im_s_vec[i].im_s,sizeof(double),1,in_file);
    s.imag=int_double(im_s_vec[i].im_s);
    im_s_vec[i].lambda_s=read_int_complex(in_file);
    im_s_vec[i].lambda_s_a=read_int_complex(in_file);
    num_fracs=calc_hurwitz1(s,N,file_N,gap,s_array,r_n_vals,no_gaps,q_start,q_end,in_file,out_file);
    if(!num_fracs)
      fatal_error("Error running Hurwitz routine. Exiting.");
    }

  fclose(in_file); // finished with lattice values
  _aligned_free(r_n_vals);
  _aligned_free(s_array);
  
  printf("Done %ld fracs and %ld s vals from Im(s)=%10.8e, Q=%ld..%ld.\n",num_fracs,num_s,im_s_vec[0].im_s,q_start,q_end);

  printf("Hur Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  no_clicks1=clock(); // start timing again
  
  fwrite(im_s_vec,sizeof(im_s),num_s,out_file);

  if(!(omegas=(int_complex *) _aligned_malloc((q_end-1)*sizeof(int_complex),16)))
    fatal_error("Error allocating memory for omegas. Exiting.\n");

  if(!(q_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for q_n. Exiting.\n");

  if(!(a_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for a_n. Exiting.\n");

  if(!(offset_vec=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
    fatal_error("Error allocating memory for offset_vec. Exiting.\n");

  init_ws(_ws);
  hur_ptr=0;

  for(q=q_start;q<=q_end;q++)
  {
    if((q&3)!=2)
    {
      make_l(q,num_s,factors,q_n,im_s_vec,offset_vec,a_n,omegas,out_file);
      hur_ptr+=factors[q].phi;
    }
  }

  fclose(out_file);

  printf("FFT Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks1)/((double) CLOCKS_PER_SEC));
  printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
  return(0);
}
