/*

File: l-func-mpfi-1.1.c

Created: 8 November 2008

Version: 1.1

Last Modified: 3rd February 2010

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 0.0 Initial implementation
V 1.0 Cleaned up summing routines
V 1.1 Multiply Gamma factors by exp(Pi t/4)
      Fix for bluecrystalp2 (use %d not %ld to format int's)
V 1.2 Use Zeta(s,1/2)=(2^s-1)Zeta(s) trick

Build instructions: gcc l-func-mpfi-1.2.c -O2 -lmpfi -lmpfr -lgmp -lm

By: DJ Platt
    Bristol University

Copyright 2008,2009,2010.

This work is funded by the UK ESPRC. */


/*

Creates a lattice of estimates for the Hurwitz Zeta function less the first
RN values of the sum.
Lattice is N points wide and NO_GAPS+1 points deep.  */

#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "../includes/mpfi.h"
#include "../includes/mpfi_io.h"
#include "../includes/mpfi_c.h"

// define this to keep track of summation errors
#define SUM_ERROR
// define this to keep track of Taylor truncation errors
#define TAYLOR_ERROR

#define SUCCESS 0
#define QUIT_CODE 0
#define FAILURE 1
#define DEBUG printf("Reached line %d.\n",__LINE__);

#define TRUE (1==1)
#define FALSE (0==1)

#define HI_PREC (1000)

// Global Variables to save stack usage

int N,M,RN,N_SUM,NO_GAPS,N_OUT,RN_OUT,NUM_DIG; 
/*
N = width of lattice
M = number points by Taylor (N-M by series) <=N
RN = number of terms not included in Hurwitz >=2
N_SUM = upper limit for Hurwitz sum > RN
NO_GAPS+1 = number of rows in lattice > Re(s)
N_OUT = width of lattice to save to file <=N
RN_OUT = RN value for saved file <=RN
*/

mpfi_c_t *r_n_vals,*s_array,*hi_prec;

#ifdef TAYLOR_ERROR
mpfi_t *taylors;
#endif

double gap;

void print_usage()
{
  printf("Usage: l-func-mpfi <prec> <N> <M> <RN> <N_SUM> <NO_GAPS> <data file> <out file> <N_OUT> <RN_OUT> <NUM_DIG>\n");
  printf("N>=M\nRN >=2\nN_SUM > RN\n");

}

int grab_memory()
{
  int i,j,lim;
  lim=N*(NO_GAPS+1);
    r_n_vals=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*lim);
    if(!r_n_vals)
      {
	printf("Fatal error allocating memory for r_n_vals. Exiting.\n");
	return(FAILURE);
      }
    for(i=0;i<lim;i++)
      mpfi_c_init(r_n_vals[i]);

    lim=N*N;
    s_array=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*lim);
    if(!s_array)
      {
	printf("Fatal error allocating memory for s_array. Exiting.\n");
	return(FAILURE);
      }
    for(i=0;i<N;i++)
      for(j=i;j<N;j++)
	mpfi_c_init(s_array[i*N+j]);

    hi_prec=(mpfi_c_t *) malloc(sizeof(mpfi_c_t)*N);
    if(!hi_prec)
      {
	printf("Fatal error allocating memory for hi_prec. Exiting.\n");
	return(FAILURE);
      }

    for(i=0;i<N;i++)
      {
	mpfi_init2(hi_prec[i]->re,HI_PREC);
	mpfi_init2(hi_prec[i]->im,HI_PREC);
      }
#ifdef TAYLOR_ERROR
    taylors=(mpfi_t *) malloc(sizeof(mpfi_t)*M);
    if(!taylors)
      {
	printf("Fatal error allocating memory for taylors. Exiting.\n");
	return(FAILURE);
      }

    for(i=0;i<M;i++)
      mpfi_init(taylors[i]);
#endif

}

mpfr_t end_point;

// initialise variables for read_float and write_mpfi_c
void init_read_float()
{
    mpfr_init(end_point);
}

inline void read_float (FILE *infile, mpfi_ptr x)
{
    mpfi_inp_str(x,infile,10);
}

inline void read_c (FILE *infile, mpfi_c_ptr z)
{
    read_float(infile,z->re);
    read_float(infile,z->im);
}

void write_float (FILE *outfile, mpfi_ptr z)
{
  double x[2];
  mpfi_get_left(end_point,z);
  x[0]=mpfr_get_d(end_point,GMP_RNDD);
  mpfi_get_right(end_point,z);
  x[1]=mpfr_get_d(end_point,GMP_RNDU);
  fwrite(x,sizeof(double),2,outfile);
}
/*
// convert complex interval to 4 doubles and output.
void write_mpfi_c (FILE *outfile, mpfi_c_ptr z)
{
    double x[4];

    mpfi_get_left(end_point,z->re);
    x[0]=mpfr_get_d(end_point,GMP_RNDD);
    mpfi_get_right(end_point,z->re);
    x[1]=mpfr_get_d(end_point,GMP_RNDU);
    mpfi_get_left(end_point,z->im);
    x[2]=mpfr_get_d(end_point,GMP_RNDD);
    mpfi_get_right(end_point,z->im);
    x[3]=mpfr_get_d(end_point,GMP_RNDU);
    fwrite(x,sizeof(double),4,outfile);
}
*/
void print_s_array()
{
  int i,j;
  printf("s_array is \n");
  for(i=0;i<N;i++)
    for(j=i;j<N;j++)
      {
	printf("i:%d j:%d ",i,j);
	mpfi_c_print(s_array[i*N+j]);
      }
}


void build_s_array(complex s)
{
    int i,j,del=1;
    double fact=2.0;

    for(i=0;i<N;i++)
    {
	mpfi_c_set_d(s_array[i*N+i],s.re+i,s.im);
	mpfi_c_mul_d(s_array[i*N+i],s_array[i*N+i],gap);
	
    }
    while(del<N)
      {
	for(i=0;i<N-del;i++)
	  {
	    mpfi_c_div_d(s_array[i*N+i+del],s_array[i*N+i],fact);
	    mpfi_c_mul(s_array[i*N+i+del],s_array[i*N+i+del],s_array[i*N+N+i+del]);
	  }
	del++;
	fact++;
      }
}



mpfi_c_t z,zeta_sum_tmp,s_tmp,m_rn_tmp,zas_tmp,hi_prec_c,wrn_tmp;
mpfi_c_t rn_sum_tmp,two_s_1,hi_prec_c1,hi_prec_c2;
mpfi_t hi_prec_i,hi_prec_i1;

#ifdef SUM_ERROR
mpfi_c_t zae_tmp3;
mpfi_t zae_tmp,zae_tmp1,zae_tmp2,zs_err;
#endif

#ifdef TAYLOR_ERROR
mpfi_t taylor_a,taylor_err;
#endif

mpfr_t outval;

void init_process_s()
{
    mpfi_c_init(z);
    mpfi_c_init(s_tmp);
    mpfi_c_init(zeta_sum_tmp);
    mpfi_c_init(m_rn_tmp);
    mpfi_c_init(zas_tmp);
    mpfi_init2(hi_prec_i,HI_PREC);
    mpfi_init2(hi_prec_i1,HI_PREC);
    mpfi_init2(hi_prec_c->re,HI_PREC);
    mpfi_init2(hi_prec_c->im,HI_PREC);
    mpfi_init2(hi_prec_c1->re,HI_PREC);
    mpfi_init2(hi_prec_c1->im,HI_PREC);
    mpfi_init2(hi_prec_c2->re,HI_PREC);
    mpfi_init2(hi_prec_c2->im,HI_PREC);
    mpfi_init2(two_s_1->re,HI_PREC);
    mpfi_init2(two_s_1->im,HI_PREC);
    mpfi_c_init(wrn_tmp);
    mpfi_c_init(rn_sum_tmp);
#ifdef SUM_ERROR
    mpfi_init(zae_tmp);
    mpfi_init(zae_tmp1);
    mpfi_init(zae_tmp2);
    mpfi_c_init(zae_tmp3);
#endif
#ifdef TAYLOR_ERROR
    mpfi_init(taylor_err);
    mpfi_init(taylor_a);
#endif
	mpfr_init(outval);
}


#ifdef SUM_ERROR
void zeta_sum_error(mpfi_c_ptr res, complex minus_s, double alpha)
{
  mpfi_set_d(zae_tmp1,alpha+(double) N_SUM);
  mpfi_set_d(zae_tmp,minus_s.re+1.0);
  mpfi_pow2(zae_tmp1,zae_tmp);
  mpfi_div(zae_tmp1,zae_tmp1,zae_tmp);
  mpfi_neg(zae_tmp,zae_tmp1);
  mpfi_put(zae_tmp1,zae_tmp);
  mpfi_inc(res->re,zae_tmp1);
  mpfi_inc(res->im,zae_tmp1);
}
#endif

void print_r_n(mpfi_c_t *r_n_vec)
{
  int i;
  for(i=0;i<5;i++)
    {
      printf("i: %5d ",i);
      mpfi_c_print(r_n_vec[i]);
    }
}


void make_r_ns(mpfi_c_t *this_r_n_vec, mpfi_c_t *last_r_n_vec, 
	       complex minus_s, double alpha)
{
  int i,j;

  for(i=0;i<M;i++)
    {
      mpfi_c_set_ui(this_r_n_vec[i],0,0);
      for(j=N-i-2;j>=0;j--)
	{
	  mpfi_c_mul(m_rn_tmp,s_array[i*(N+1)+j],last_r_n_vec[j+i+1]);
	  mpfi_c_inc(this_r_n_vec[i],m_rn_tmp);
	}
      mpfi_c_inc(this_r_n_vec[i],last_r_n_vec[i]);
#ifdef TAYLOR_ERROR
      mpfi_add(this_r_n_vec[i]->re,this_r_n_vec[i]->re,taylors[i]);
      mpfi_add(this_r_n_vec[i]->im,this_r_n_vec[i]->im,taylors[i]);
#endif
    }
}


void make_r_ns_neg(mpfi_c_t *this_r_n_vec, mpfi_c_t *last_r_n_vec, 
		   complex minus_s, double alpha)
{
  int i,j,sign;

  for(i=0;i<M;i++)
    {
      mpfi_c_set_ui(this_r_n_vec[i],0,0);
      if((N-i)&1)
	sign=1;
      else
	sign=-1;
      for(j=N-i-2;j>=0;j--)
	{
	  mpfi_c_mul(m_rn_tmp,s_array[i*(N+1)+j],last_r_n_vec[j+i+1]);
	  if(sign==1)
	    mpfi_c_inc(this_r_n_vec[i],m_rn_tmp);
	  else
	    mpfi_c_dec(this_r_n_vec[i],m_rn_tmp);
	  sign=-sign;
	}
      mpfi_c_inc(this_r_n_vec[i],last_r_n_vec[i]);
#ifdef TAYLOR_ERROR
      mpfi_add(this_r_n_vec[i]->re,this_r_n_vec[i]->re,taylors[i]);
      mpfi_add(this_r_n_vec[i]->im,this_r_n_vec[i]->im,taylors[i]);
#endif
    }
}


// sets up top and bottom row of lattice
// called with hi_prec[j] set to zeta(1+j+it), j=0..N
// assumes rn>=2
// sets up top and bottom row of lattice
void r_n_top_bottom_rows(complex minus_s)
{
  int i,j,ptr,col,bot=N*NO_GAPS;
  complex my_s;
  my_s.re=-minus_s.re;
  my_s.im=-minus_s.im;

  mpfi_const_log2(hi_prec_i1);
  mpfi_mul_d(hi_prec_i,hi_prec_i1,my_s.im);
  mpfi_cos(hi_prec_c1->re,hi_prec_i);
  mpfi_sin(hi_prec_c1->im,hi_prec_i); // 2^(it)

  for(i=0;i<N;i++)
    {
      mpfi_mul_d(hi_prec_i,hi_prec_i1,my_s.re);
      mpfi_exp(hi_prec_i,hi_prec_i);
      mpfi_mul(two_s_1->re,hi_prec_c1->re,hi_prec_i);
      mpfi_mul(two_s_1->im,hi_prec_c1->im,hi_prec_i);
      mpfi_c_dec_ui(two_s_1,1); // 2^s -1
      //printf("2^s-1=");mpfi_c_print(two_s_1);
      mpfi_mul(hi_prec_c->re,two_s_1->re,hi_prec[i]->re);
      mpfi_mul(hi_prec_i,two_s_1->im,hi_prec[i]->im);
      mpfi_sub(hi_prec_c->re,hi_prec_c->re,hi_prec_i);
      mpfi_mul(hi_prec_c->im,two_s_1->re,hi_prec[i]->im);
      mpfi_mul(hi_prec_i,two_s_1->im,hi_prec[i]->re);
      mpfi_add(hi_prec_c->im,hi_prec_c->im,hi_prec_i); // Zeta(s)*(2^s-1)=Zeta(s,1/2)

      //mpfi_c_print(hi_prec_c);if(i==3) exit(0);
      
      for(j=0;j<RN;j++)
	{
	  mpfi_set_d(hi_prec_i,0.5+j);
	  mpfi_log(hi_prec_i,hi_prec_i);
	  mpfi_mul_d(hi_prec_c2->re,hi_prec_i,-my_s.im);
	  mpfi_mul_d(hi_prec_i,hi_prec_i,-my_s.re);
	  mpfi_exp(hi_prec_i,hi_prec_i);                         // (0.5+j)^(-x)
	  mpfi_sin(hi_prec_c2->im,hi_prec_c2->re);
	  mpfi_cos(hi_prec_c2->re,hi_prec_c2->re);
	  mpfi_mul(hi_prec_c2->re,hi_prec_c2->re,hi_prec_i);
	  mpfi_mul(hi_prec_c2->im,hi_prec_c2->im,hi_prec_i);
	  //printf("subtracting ");mpfi_c_print(hi_prec_c2);if(j==2) exit(0);
	  mpfi_c_sub(hi_prec_c,hi_prec_c,hi_prec_c2);
	}
      mpfi_c_set(r_n_vals[N*NO_GAPS/2+i],hi_prec_c);
      //mpfi_c_print(r_n_vals[N*NO_GAPS/2+i]);
      my_s.re++;
    }

 
  ptr=0;
  for(i=0;i<N;i++)
    {
      mpfi_c_dec_ui(hi_prec[i],1);
      //printf("%d ",-mpfi_rel_error(hi_prec[i]->re));
    }

  for(i=2;i<RN;i++)
    {
      mpfi_set_ui(hi_prec_i,i);
      mpfi_log(hi_prec_i,hi_prec_i);
      mpfi_set(hi_prec_i1,hi_prec_i);
      mpfi_mul_d(hi_prec_i,hi_prec_i,minus_s.re);
      mpfi_exp(hi_prec_i,hi_prec_i);
      mpfi_mul_d(hi_prec_i1,hi_prec_i1,minus_s.im);
      mpfi_sin(hi_prec_c->im,hi_prec_i1);
      mpfi_cos(hi_prec_c->re,hi_prec_i1);
      mpfi_mul(hi_prec_c->im,hi_prec_c->im,hi_prec_i);
      mpfi_mul(hi_prec_c->re,hi_prec_c->re,hi_prec_i);  // hp_c = i^(-s)
      mpfi_c_sub(hi_prec[0],hi_prec[0],hi_prec_c);
      for(j=1;j<N;j++)
	{
	  mpfi_div_ui(hi_prec_c->re,hi_prec_c->re,i);
	  mpfi_div_ui(hi_prec_c->im,hi_prec_c->im,i);
	  mpfi_c_sub(hi_prec[j],hi_prec[j],hi_prec_c);
	}
    }

  for(i=0;i<N;i++)
    mpfi_c_set(r_n_vals[NO_GAPS*N+i],hi_prec[i]); // bottom row done


  mpfi_set_ui(hi_prec_i,RN);
  mpfi_log(hi_prec_i,hi_prec_i);
  mpfi_set(hi_prec_i1,hi_prec_i);
  mpfi_mul_d(hi_prec_i,hi_prec_i,minus_s.re);
  mpfi_exp(hi_prec_i,hi_prec_i);
  mpfi_mul_d(hi_prec_i1,hi_prec_i1,minus_s.im);
  mpfi_sin(hi_prec_c->im,hi_prec_i1);
  mpfi_cos(hi_prec_c->re,hi_prec_i1);
  mpfi_mul(hi_prec_c->im,hi_prec_c->im,hi_prec_i);
  mpfi_mul(hi_prec_c->re,hi_prec_c->re,hi_prec_i);  // hp_c = RN^(-s)
  mpfi_c_sub(hi_prec[0],hi_prec[0],hi_prec_c);
  for(j=1;j<N;j++)
    {
      mpfi_div_ui(hi_prec_c->re,hi_prec_c->re,RN);
      mpfi_div_ui(hi_prec_c->im,hi_prec_c->im,RN);
      mpfi_c_sub(hi_prec[j],hi_prec[j],hi_prec_c);
    }

  for(i=0;i<N;i++)
    mpfi_c_set(r_n_vals[i],hi_prec[i]);

}

// only called if N>M
// fill r_n columns M->N-1 with sum approximations to Hrn(s,alpha)
// s=(0.5+j)+it for column j,
void do_r_n_sums(double im_s, double alpha, mpfi_c_t *rn)
{
  int i,n;
  complex minus_s;
  double beta=alpha+(double) RN;

  minus_s.re=-0.5-M;
  minus_s.im=-im_s;

  for(i=M;i<N;i++)
    mpfi_c_set_ui(rn[i],0,0);

  for(n=RN;n<=N_SUM;n++)
    {
      mpfi_c_pow_c(rn_sum_tmp,beta,minus_s);
      mpfi_c_add(rn[M],rn[M],rn_sum_tmp);
      for(i=M+1;i<N;i++)
	{
	  mpfi_c_div_d(rn_sum_tmp,rn_sum_tmp,beta);
	  mpfi_c_add(rn[i],rn[i],rn_sum_tmp);
	}
      beta++;
    }

#ifdef SUM_ERROR
  for(i=M;i<N;i++)
    {
      zeta_sum_error(rn[i],minus_s,alpha);
      minus_s.re-=1.0;
    }
#endif
}

void write_mpfi_s(FILE *outfile,mpfi_c_ptr z)
{
	
	mpfi_get_fr(outval,z->re);
	mpfr_out_str(outfile,10,NUM_DIG,outval,GMP_RNDN);
	fprintf(outfile,"\n");
	mpfi_get_fr(outval,z->im);
	mpfr_out_str(outfile,10,NUM_DIG,outval,GMP_RNDN);
	fprintf(outfile,"\n");
}

// add back sum (n+alpha)^(minus_s-j) to rn[j]
// n=RN_OUT..RN-1
void write_rn(FILE *outfile, mpfi_c_t *row, double alpha, complex minus_s)
{
  int n,i;
  double beta=(double)RN_OUT+alpha;

  for(n=RN_OUT;n<RN;n++)
    {
      mpfi_c_pow_c(wrn_tmp,beta,minus_s);
      mpfi_c_add(row[0],row[0],wrn_tmp);
      for(i=1;i<N_OUT;i++)
	{
	  mpfi_c_div_d(wrn_tmp,wrn_tmp,beta);
	  mpfi_c_add(row[i],row[i],wrn_tmp);
	}
      beta++;
    }

  for(i=0;i<N_OUT;i++)
    write_mpfi_s(outfile,row[i]);
}

char first=TRUE;
char second=FALSE;
double first_im_s;

int process_s(FILE *infile, FILE *outfile, int num_zeta)
{
  double alpha;
  int i,j,ptr;
  complex s,minus_s,minus_s_m;

  s.re=1.0;

  fscanf(infile,"%lg",&s.im);
  //fwrite(&s.im,sizeof(double),1,outfile);
  if(second)
    {
      second=FALSE;
      if((s.im-first_im_s)!=5.0/64.0)
	{
	  printf("Step size not 5/64 as expected. Exiting.\n");
	  exit(1);
	}
    }
  if(first)
  {
    printf("Processing Im(s)=%20.18e\n",s.im);
    first=FALSE;
    first_im_s=s.im;
    second=TRUE;
  }


  // two values of gamma function, just read and output
  // not needed here but used in next phase
  read_c(infile,z);  // Gamma(s/2)
  mpfi_mul_d(taylor_err,mpfi_pi,s.im/4.0);
  mpfi_exp(taylor_err,taylor_err);
  mpfi_c_mul_i(z,z,taylor_err);
  //write_mpfi_c(outfile,z);

  read_c(infile,z); // Gamma((s+1)/2)
  mpfi_mul_d(taylor_err,mpfi_pi,s.im/4.0);
  mpfi_exp(taylor_err,taylor_err);
  mpfi_c_mul_i(z,z,taylor_err);
  //write_mpfi_c(outfile,z);

  // calculate pi^(-it/2)
  mpfi_log(taylor_err,mpfi_pi); // log(pi)
  mpfi_mul_d(taylor_err,taylor_err,-s.im/2.0); // -tlog(pi)/2
  mpfi_sin(z->im,taylor_err);
  mpfi_cos(z->re,taylor_err);
  //write_mpfi_c(outfile,z);
  //mpfi_c_print(z);


//printf("%d\n",N);
//printf("here1\n");

	for(i=0;i<N;i++)
    {
      // nasty frig for debugging
	//printf("%d\n",i);
      read_c(infile,hi_prec[i]); // = zeta(s+i)
	  //mpfi_c_print(hi_prec[i]);
	  //printf("%d\n",i);
      //mpfi_c_set(r_n_vals[i],r_n_vals[NO_GAPS*N+i]);
    }

  //printf("here2\n");
    
  minus_s=c_neg(s);
  //printf("here3\n");  
  r_n_top_bottom_rows(minus_s);
  //printf("here4\n");

  //printf("num_zeta=%d\n",num_zeta);
  for(i=N;i<num_zeta;i++)                 // throw away extras
    read_c(infile,z);

  //printf("here5\n");
  build_s_array(s);
  //printf("here6\n");

#ifdef TAYLOR_ERROR
  // set taylor_err
  mpfi_set_ui(taylor_a,RN-1);
  mpfi_log(taylor_a,taylor_a);
  mpfi_mul_d(taylor_a,taylor_a,0.5-(double)N);
  mpfi_exp(taylor_a,taylor_a);
  mpfi_div_d(taylor_a,taylor_a,(double)N-0.5); // a=((RN-1)^(-N+0.5))/(N-0.5)

  mpfi_add_d(taylor_err,taylor_err,s.im+0.5+2.0+(double)N);
  mpfi_mul_d(taylor_err,taylor_err,gap);  // r>=|s+N+2|*delta
  mpfi_div_ui(taylor_err,taylor_err,N+2); // r>=|s+N+2|*delta/(N+2)
  mpfi_dec_ui(taylor_err,1);   // r-1
  mpfi_div(taylor_err,taylor_a,taylor_err); // a/(r-1)
  for(i=0;i<M;i++)
    {
      mpfi_c_mod(taylors[i],s_array[(i+1)*N-1]);
      mpfi_mul(taylors[i],taylors[i],taylor_err);
      mpfi_neg(taylor_a,taylors[i]);
      mpfi_put(taylors[i],taylor_a);
      //printf("Taylor error [%d]=",i);mpfi_print(taylors[i]);printf("\n");
    }

#endif
//  printf("here7\n");


  if(N!=M) // there are some terms to sum
    {
      for(i=1;i<NO_GAPS;i++)
	  {
	    if(i==NO_GAPS/2)
	      continue;
	    //printf("%d\n",i);
	    do_r_n_sums(s.im,1.0-gap*i,&r_n_vals[N*i]);
	  }
    }

    //printf("here8\n");

  for(i=1;i<=NO_GAPS/4;i++)
  {
	  //if((i%100)==0)
		  //printf("%d\n",i);
	    make_r_ns(&r_n_vals[i*N],&r_n_vals[(i-1)*N],minus_s,1.0-gap*i);
  }
  for(i=NO_GAPS/2+1;i<=3*NO_GAPS/4;i++)
  {
	  //if((i%100)==0)
		  //printf("%d\n",i);
	    make_r_ns(&r_n_vals[i*N],&r_n_vals[(i-1)*N],minus_s,1.0-gap*i);
  }
    
  for(i=NO_GAPS-1;i>3*NO_GAPS/4;i--)
    make_r_ns_neg(&r_n_vals[i*N],&r_n_vals[(i+1)*N],minus_s,1.0-gap*i);

  for(i=NO_GAPS/2-1;i>NO_GAPS/4;i--)
    make_r_ns_neg(&r_n_vals[i*N],&r_n_vals[(i+1)*N],minus_s,1.0-gap*i);

  /*      
  for(i=0;i<=20;i++)
  printf(" i  ");
  printf("\n");
  for(j=0;j<N;j+=3)
  printf("%4d",j);
  printf("\n");
  for(i=0;i<=NO_GAPS;i=i+NO_GAPS/16)
  {
  printf("%4d",i);
  for(j=0;j<N;j+=3)
  printf(" %3d",-mpfi_rel_error(r_n_vals[i*N+j]->re));
  printf("\n");
  }
  printf("H_rn(s,0.5)=");mpfi_c_print(r_n_vals[N*NO_GAPS/2]);    
    
  printf("saving...\n");

  */  
  
  for(i=0;i<=NO_GAPS;i++)
    write_rn(outfile,&r_n_vals[i*N],1.0-i*gap,minus_s);
  /*
    printf("rel_errors now\n");
    printf(" i  ");
    for(j=0;j<N;j+=3)
    printf("%4d",j);
    printf("\n");
    for(i=0;i<=NO_GAPS;i=i+NO_GAPS/16)
    {
    printf("%4d",i);
    for(j=0;j<N;j+=3)
    printf(" %3d",-mpfi_rel_error(r_n_vals[i*N+j]->re));
    printf("\n");
    }
  */  
  printf("H_1(s,0.25)=");mpfi_c_print(r_n_vals[3*N*NO_GAPS/4]);
  printf("H_1(s,0.5)=");mpfi_c_print(r_n_vals[N*NO_GAPS/2]);    
  printf("H_1(s,0.75)=");mpfi_c_print(r_n_vals[N*NO_GAPS/4]);
  
  for(i=0;i<N_OUT;i++)
    printf(" %d",-mpfi_rel_error(r_n_vals[N*NO_GAPS/4+i]->re));
  printf(".\n");
  //exit(0);
  return(SUCCESS);
}

int main(int argc, char **argv)
{

    int prec,num_s,num_zeta,i;
    int counter;
    FILE *infile,*outfile;

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=12)
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
    
    mpfi_c_setup(prec);

    N=atoi(argv[2]);
    if(N<=0)
      {
	print_usage();
	return(QUIT_CODE);
      }

    M=atoi(argv[3]);
    if((M>N)||(M<=0))
      {
	print_usage();
	return(QUIT_CODE);
      };

    RN=atoi(argv[4]);
    if(RN<=1)
      {
	print_usage();
	return(QUIT_CODE);
      }

    if(N!=M)
      {
	N_SUM=atoi(argv[5]);
	if(N_SUM<=RN)
	  {
	    print_usage();
	    return(QUIT_CODE);
	  }
      }
    else
      N_SUM=0;

    NO_GAPS=atoi(argv[6]);
    if(NO_GAPS<=0)
      {
	print_usage();
	return(QUIT_CODE);
      }

    infile=fopen(argv[7],"r");
    if(!infile)
      {
	printf("Failed to open file %s for input. Exiting.\n",argv[7]);
	return(QUIT_CODE);
      }

    outfile=fopen(argv[8],"w");
    if(!outfile)
      {
	printf("Failed to open file %s for output. Exiting.\n",argv[8]);
	return(QUIT_CODE);
      }

    N_OUT=atoi(argv[9]);
    if((N_OUT<=0)||(N_OUT>N))
      {
	printf("Problem with value of N_OUT. Exiting.\n");
	return(QUIT_CODE);
      }

    RN_OUT=atoi(argv[10]);
    if((RN_OUT<0)||(RN_OUT>RN))
      {
	printf("Problem with value of RN_OUT. Exiting.\n");
	return(QUIT_CODE);
      }
    NUM_DIG=atoi(argv[11]);
    if(NUM_DIG<=0)
      {
	printf("Problem with value of NUM_DIG. Exiting.\n");
	return(QUIT_CODE);
      }

    if(grab_memory()==FAILURE)
      return(QUIT_CODE);
    printf("Running with prec=%d N=%d M=%d RN=%d N_SUM=%d NO_GAPS=%d.\n",
	   prec,N,M,RN,N_SUM,NO_GAPS);

    printf("N_OUT=%d RN_OUT=%d.\n",N_OUT,RN_OUT);

    fscanf(infile,"%d",&num_s);
    fscanf(infile,"%d",&num_zeta);

    //fwrite(&num_s,sizeof(int),1,outfile);
    //fwrite(&N_OUT,sizeof(int),1,outfile);
    //fwrite(&RN_OUT,sizeof(int),1,outfile);
    //fwrite(&NO_GAPS,sizeof(int),1,outfile);

    if(num_zeta<=N)
    {
      printf("Insufficient Zeta Terms %d in infile to support width of lattice. Exiting.\n",num_zeta);
	return(QUIT_CODE);
    }

    gap=1.0/NO_GAPS;

    init_read_float();
    init_process_s();

    for(i=0;i<num_s;i++)
      if(process_s(infile,outfile,num_zeta)!=SUCCESS)
	return(QUIT_CODE);

    fclose(infile);
    fclose(outfile);
    printf("Completed.\n");
    return(QUIT_CODE);
}
