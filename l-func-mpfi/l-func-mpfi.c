/*

File: l-func-mpfi.c

Created: 8 November 2008

Version: 1.0

Last Modified: 14 November 2008

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes: Assumes bugs in trig functions of MPFI
                      are fixed in this release.

V 1.0 Initial implementation

Build instructions: gcc -ol-func-mpfi l-func-mpfi.c -O2 -lmpfi -lmpfr -lgmp -lm

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */


/*

Creates a lattice of estimates for the Hurwitz Zeta function less the first
RN values of the sum.
Lattice is N points wide and NO_GAPS+1 points deep.  */

#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
typedef struct{double re;double im;} complex;
#include "../includes/mpfi_c.h"

// define this to keep track of summation errors
// #define SUM_ERROR
// define this to keep track of Taylor truncation errors
#define TAYLOR_ERROR

#define SUCCESS 0
#define QUIT_CODE 0
#define FAILURE 1
#define DEBUG printf("Reached line %d.\n",__LINE__);

#define TRUE (1==1)
#define FALSE (0==1)


inline complex cmplx (double re, double im)
{
  complex z;
  z.re=re;
  z.im=im;
  return(z);
}

inline complex c_neg (complex z)
{
  complex z1;
  z1.re=-z.re;
  z1.im=-z.im;
  return(z1);
}

void complex_print(complex s)
{
  printf("%10.8e+i%10.8e\n",s.re,s.im);
}

// Global Variables to save stack usage

int N,M,RN,N_SUM,NO_GAPS; // N = width of lattice
                          // M = number points by Taylor (N-M by series)
                          // RN = number of terms not included in Hurwitz
                          // N_SUM = upper limit for Hurwitz sum
                          // NO_GAPS+1 = number of rows in lattice

mpfi_c_t *r_n_vals,*s_array;
double gap;

void print_usage()
{
  printf("Usage: l-func-mpfi <prec> <N> <M> <RN> <N_SUM> <NO_GAPS> <data file> <out file>\n");
  printf("RN >=2\nN_SUM > RN\n");

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


mpfi_c_t z,zeta_sum_tmp,s_tmp,m_rn_tmp,zas_tmp;
mpfi_c_t zae_tmp3;
mpfi_t zae_tmp,zae_tmp1,zae_tmp2,zs_err;

#ifdef TAYLOR_ERROR
mpfi_t max_hn,s_tmp1;
#endif

void init_process_s()
{
    mpfi_c_init(z);
    mpfi_c_init(s_tmp);
    mpfi_c_init(zeta_sum_tmp);
    mpfi_c_init(m_rn_tmp);
    mpfi_c_init(zas_tmp);
    mpfi_init(zae_tmp);
    mpfi_init(zae_tmp1);
    mpfi_init(zae_tmp2);
    mpfi_c_init(zae_tmp3);
#ifdef TAYLOR_ERROR
    mpfi_init(max_hn);
    mpfi_init(s_tmp1);
#endif
}

void zeta_sum(mpfi_c_ptr res, complex minus_s, int start, int end)
{
  int i;

  mpfi_c_set_ui(res,0,0);
  for(i=start;i<=end;i++)
    {
      mpfi_c_pow_c(zeta_sum_tmp,i,minus_s);
      mpfi_c_inc(res,zeta_sum_tmp);
    }
}

#ifdef TAYLOR_ERROR
// calculate error due to Taylor approximation
// approximating H_{RN}(s,alpha-delta) with N-j terms
void taylor_error(mpfi_c_ptr res, complex s, double delta, unsigned int j)
{
  mpfi_c_mod(zae_tmp,s_array[N*j+N-1]);
  mpfi_mul(zae_tmp,zae_tmp,max_hn);
  mpfi_mul_ui(zae_tmp1,zae_tmp,N-j+1);
  mpfi_set_d(zae_tmp,delta);
  mpfi_mul_d(zae_tmp,zae_tmp,s.im+(double)j+(double)N+0.5);
  mpfi_sub_ui(zae_tmp,zae_tmp,N-j+1);
  mpfi_div(zae_tmp1,zae_tmp1,zae_tmp);
  mpfi_neg(zae_tmp,zae_tmp1);
  mpfi_put(zae_tmp1,zae_tmp);
  //printf("taylor error for Re(s)=%3.1e n=%d is ",s.re,n);
  //mpfi_print(zae_tmp1);
  mpfi_inc(res->re,zae_tmp1);
  mpfi_inc(res->im,zae_tmp1);
  //if(n==(N-M+1))
  //  exit(0);
}
#endif

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

void zeta_alpha_sum(mpfi_c_ptr res, complex minus_s, double alpha)
{
  int i;

  mpfi_c_set_ui(res,0,0);
  for(i=RN;i<=N_SUM;i++)
    {
      mpfi_c_pow_c(zas_tmp,i+alpha,minus_s);
      mpfi_c_add(res,res,zas_tmp);
    }
}

void print_r_n(mpfi_c_t *r_n_vec)
{
  int i;
  for(i=0;i<N;i++)
    {
      printf("i: %5d ");
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
	  //printf("multiplying ");mpfi_c_print(s_array[i*(N+1)+j]);
	  //printf("by ");mpfi_c_print(last_r_n_vec[j+i+1]);
	  mpfi_c_mul(m_rn_tmp,s_array[i*(N+1)+j],last_r_n_vec[j+i+1]);
	  //printf("getting ");mpfi_c_print(m_rn_tmp);
	  mpfi_c_inc(this_r_n_vec[i],m_rn_tmp);
	  //printf("Hrn now=");mpfi_c_print(this_r_n_vec[i]);
	}
      mpfi_c_inc(this_r_n_vec[i],last_r_n_vec[i]);
#ifdef TAYLOR_ERROR
      taylor_error(this_r_n_vec[i],cmplx(0.5+(double) i,-minus_s.im) ,gap,i);
#endif
    }

  minus_s.re-=M;

  for(i=M;i<N;i++)
    {
      zeta_alpha_sum(this_r_n_vec[i],minus_s,alpha);
#ifdef SUM_ERROR
      zeta_sum_error(this_r_n_vec[i],minus_s,alpha);
#endif
      minus_s.re--;
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
	  //if(sign==1) printf("adding "); else printf("subtracting ");
	  //printf("new term ");mpfi_c_print(m_rn_tmp);
	  //printf("gives ");
	  if(sign==1)
	    mpfi_c_inc(this_r_n_vec[i],m_rn_tmp);
	  else
	    mpfi_c_dec(this_r_n_vec[i],m_rn_tmp);
	  sign=-sign;
	  //mpfi_c_print(this_r_n_vec[i]);
	}
      mpfi_c_inc(this_r_n_vec[i],last_r_n_vec[i]);
#ifdef TAYLOR_ERROR
      taylor_error(this_r_n_vec[i],cmplx(0.5+(double) i,-minus_s.im) ,gap,i);
#endif    
    }
  minus_s.re-=M;
  for(i=M;i<N;i++)
    {
      zeta_alpha_sum(this_r_n_vec[i],minus_s,alpha);
#ifdef SUM_ERROR
      zeta_sum_error(this_r_n_vec[i],minus_s,alpha);
#endif
      minus_s.re--;
    }
}

#ifdef TAYLOR_ERROR
int first=TRUE;
#endif

int process_s(FILE *infile, FILE *outfile, int num_zeta, complex s)
{
    double alpha;
    int i,j,ptr;
    complex minus_s;

    fscanf(infile,"%lg",&s.im);
    fwrite(&s.im,sizeof(double),1,outfile);

// two values of gamma function, just read and output
// not needed here but used in next phase
    read_c(infile,z);
    write_mpfi_c(outfile,z);
    read_c(infile,z);
    write_mpfi_c(outfile,z);

    printf("in process_s with s=");complex_print(s);

    minus_s=c_neg(s);
    for(i=0;i<N;i++)
      {
	read_c(infile,r_n_vals[NO_GAPS*N+i]); // bottom row = zeta(s+i)-1
	zeta_sum(s_tmp,minus_s,2,RN-1);
	mpfi_c_dec(r_n_vals[NO_GAPS*N+i],s_tmp);

	mpfi_c_pow_c(s_tmp,RN,minus_s);
	mpfi_c_sub(r_n_vals[i],r_n_vals[NO_GAPS*N+i],s_tmp);
	
	minus_s.re--;
      }
#ifdef TAYLOR_ERROR
    if(first)
      {
	first=FALSE;
	read_c(infile,z);   // z=zeta(N+0.5)
	mpfi_set(max_hn,z->re);
	for(i=2;i<=RN;i++)
	  {
	    mpfi_set_ui(s_tmp1,i);
	    mpfi_log(s_tmp1,s_tmp1);
	    mpfi_mul_d(s_tmp1,s_tmp1,-(double)N-0.5);
	    mpfi_exp(s_tmp1,s_tmp1);
	    mpfi_sub(max_hn,max_hn,s_tmp1);
	  }
	printf("max_hn set to ");mpfi_print(max_hn);

	for(i=N+1;i<num_zeta;i++)
	  read_c(infile,z);
      }
    else
#endif
      for(i=N;i<num_zeta;i++)                 // throw away extras
	read_c(infile,z);

    build_s_array(s);
    //    print_s_array();
    
    minus_s=c_neg(s);

    //for(i=1;i<=20;i++)  // debugging stuff
    for(i=1;i<=NO_GAPS/2;i++)
      {
	if((i%128)==0)
	  printf(".");
	make_r_ns(&r_n_vals[i*N],&r_n_vals[(i-1)*N],minus_s,1.0-gap*i);
      }
    
    for(i=NO_GAPS-1;i>NO_GAPS/2;i--)
      {
	if((i%128)==0)
	  printf(".");
	make_r_ns_neg(&r_n_vals[i*N],&r_n_vals[(i+1)*N],minus_s,1.0-gap*i);
      }
    
    printf("\n");

    //for(i=0;i<=20;i++)
    for(i=0;i<=NO_GAPS;i=i+NO_GAPS/16)
      {
	printf("i:%5d ",i);
	for(j=0;j<N;j++)
	  printf(" %2d",-mpfi_rel_error(r_n_vals[i*N+j]->re));
	printf("\n");
      }
    /*
      printf("i:%5d err:%4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n",
	     i,mpfi_rel_error(r_n_vals[i*N]->re),
	     mpfi_rel_error(r_n_vals[i*N+1]->re),
	     mpfi_rel_error(r_n_vals[i*N+2]->re),
	     mpfi_rel_error(r_n_vals[i*N+3]->re),
	     mpfi_rel_error(r_n_vals[i*N+4]->re),
	     mpfi_rel_error(r_n_vals[i*N+5]->re),
	     mpfi_rel_error(r_n_vals[i*N+6]->re),
	     mpfi_rel_error(r_n_vals[i*N+7]->re),
	     mpfi_rel_error(r_n_vals[i*N+8]->re),
	     mpfi_rel_error(r_n_vals[i*N+9]->re));
    */
    printf("H_rn(s,0.5)=");mpfi_c_print(r_n_vals[N*NO_GAPS/2]);

 
    /*
    printf("first row\n");
    print_r_n(r_n_vals);
    printf("second row\n");
    print_r_n(&r_n_vals[N]);
    //printf("next row\n");
    //print_r_n(&r_n_vals[N*NO_GAPS/2+N]);
    //printf("last row\n");
    //print_r_n(&r_n_vals[N*NO_GAPS]);
    */

    ptr=0;
    /*
    for(i=0;i<=NO_GAPS;i++)
      for(j=0;j<N;j++)
	write_mpfi_c(outfile,r_n_vals[ptr++]);
    */
    return(SUCCESS);
}



int main(int argc, char **argv)
{

    int prec,num_s,num_zeta,i;
    complex s;
    FILE *infile,*outfile;



/*  g_rho_sum will contain the sum of all the G(rho) terms, plus the
    error because we have ignored zeros > NUM_ZEROS, plus the error
    because we assume G(biggest rho) is zero.

    g1 will contain our estimate of G(1)

    ln2 will contain natural log(2)

    err_term is the bounds for the integral of F(s) with Im(s)=-1

    if we add all the lower bounds and all the upper bounds of these real
    intervals, we get the interval containing our result. */

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=9)
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
    if(RN<1)
      {
	print_usage();
	return(QUIT_CODE);
      }

    N_SUM=atoi(argv[5]);
    if(N_SUM<=RN)
      {
	print_usage();
	return(QUIT_CODE);
      }

    NO_GAPS=atoi(argv[6]);
    if(NO_GAPS<=0)
      {
	print_usage();
	return(QUIT_CODE);
      }

    infile=fopen(argv[7],"r");
    if(!infile)
      {
	printf("Failed to open file %s for input. Exiting.\n",argv[5]);
	return(QUIT_CODE);
      }

    outfile=fopen(argv[8],"wb");
    if(!outfile)
      {
	printf("Failed to open file %s for output. Exiting.\n",argv[6]);
	return(QUIT_CODE);
      }

    if(grab_memory()==FAILURE)
      return(QUIT_CODE);
    printf("Running with prec=%d N=%d M=%d RN=%d N_SUM=%d NO_GAPS=%d.\n",
	   prec,N,M,RN,N_SUM,NO_GAPS);


    fscanf(infile,"%ld",&num_s);
    fscanf(infile,"%ld",&num_zeta);
    fwrite(&num_s,sizeof(int),1,outfile);
    fwrite(&N,sizeof(int),1,outfile);
    fwrite(&RN,sizeof(int),1,outfile);
    fwrite(&NO_GAPS,sizeof(int),1,outfile);
    if(num_zeta<=N)
    {
	printf("Insufficient Zeta Terms in infile to support width of lattice. Exiting.\n");
	return(QUIT_CODE);
    }

    gap=1.0/NO_GAPS;

    init_read_float();
    init_process_s();
    s.re=0.5;

    for(i=0;i<num_s;i++)
	if(process_s(infile,outfile,num_zeta,s)!=SUCCESS)
	    return(QUIT_CODE);

    fclose(infile);
    fclose(outfile);



    return(QUIT_CODE);
}
