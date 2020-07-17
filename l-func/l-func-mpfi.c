/*

File: l-func-mpfi.c

Created: 6th July 2008

Version: <v> = 1.0

Last Modified: 21st August 2008

Dialect: C

Requires: GMP v. 4.2.2
          MPFR v. 2.3.1
          MPFI 1.3.4-RC3

Implementation notes: 

Build instructions: gcc -ol-func-mpfi l-func-mpfi.c -O3 -lmpfi -lmpfr -lgmp

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "1.0"

#include "stdio.h"
#include "math.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"
#include "../g/mpfi_c.c"

#define SUCCESS (1)
#define FAILURE (0)
#define TRUE (1)
#define FALSE (0)
#define OUT_VEC_SIZE (1<<20)
#define debug {printf("%d\n",__LINE__);}

#define PREC (53)

#define print_usage {\
  printf("Usage: l-func%s <q_start> <q_end> <h_terms> ",VERSION);\
  printf("<N_terms> <gap> <zeta_file> <outfile>\n");\
  printf("  (q_start)   - integer > 2\n");\
  printf("  (q_end)     - integer > q_start\n");\
  printf("  (h_terms)   - integer Hurwitz terms, 0<h_terms<=max defined in zeta_file\n");\
  printf("  (N_terms)   - integer sum terms > 0\n");\
  printf("  (gap)       - 0.0 < gap < 1.0\n");\
  printf("  (zeta_file) - file containing values of zeta(s+n).\n");\
  printf("  (outfile)   - output file.\n");\
  return(SUCCESS);\
};


inline int gcd (unsigned int a, unsigned int b)
     // Euclid algorithm gcd
{
  unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

inline int co_prime(unsigned int a, unsigned int b)
{
  return(gcd(a,b)==1);
};


mpfi_t pow_tmp1;

void mpfi_pow(mpfi_ptr res, mpfi_ptr base, mpfi_ptr power)
{
  mpfi_log(pow_tmp1,base);
  mpfi_mul(pow_tmp1,pow_tmp1,power);
  mpfi_exp(res,pow_tmp1);
};

void mpfi_c_pow (mpfi_c_ptr res, mpfi_c_ptr base, mpfi_c_ptr power)
{
};

mpfi_t pow_tlni,pow_i_to_sigma,pow_sin,pow_cos;


void mpfi_i_c_pow(mpfi_c_t res,
		   mpfi_ptr x, mpfi_c_ptr s)
     // return i^s 
{
  mpfi_log(pow_tlni,x);        // tlni<-log(x)
  mpfi_mul(pow_tlni,s->im,pow_tlni);  //tlni<-tlnx
  mpfi_pow(pow_i_to_sigma,x,s->re); // i_to_sigma<-x^sigma
  mpfi_sin(pow_sin,pow_tlni);
  mpfi_cos(pow_cos,pow_tlni);
  mpfi_mul(pow_cos,pow_cos,pow_i_to_sigma);
  mpfi_mul(pow_sin,pow_sin,pow_i_to_sigma);
  mpfi_c_set_re(res,pow_cos);
  mpfi_c_set_im(res,pow_sin);
  
};

void mpfi_ui_c_pow(mpfi_c_t res, unsigned int i,
		   mpfi_c_ptr s)
{
  mpfi_set_ui(pow_tmp1,i);
  mpfi_i_c_pow(res,pow_tmp1,s);
};



/*
dcomplex s_n(dcomplex s, double alpha, int N)
     // return sum from 0 to N of (n+alpha)^(-s) 
     // N will probably end up being 0 or 1 so might
     // revisit 
{
  int n;
  dcomplex res(0,0);

  for(n=0;n<N;n++)
    res+=pow(n+alpha,-s);
  //  cout << "s_n called with " << s << " " << alpha << " " << N << " returning " << res << endl;
  return(res);
};
*/

mpfi_c_t s_n_tmp1,s_n_res;
mpfi_t s_n_tmp2;

void s_n(mpfi_c_ptr res, mpfi_c_ptr s, mpfi_ptr alpha, unsigned int N)
{
  int n;
  mpfi_c_set_ui(s_n_res,0,0);
  mpfi_set(s_n_tmp2,alpha);
  for(n=0;n<N;n++)
    {
      mpfi_i_c_pow(s_n_tmp1,s_n_tmp2,s);
      mpfi_c_add(s_n_res,s_n_res,s_n_tmp1);
      mpfi_add_ui(s_n_tmp2,s_n_tmp2,1);
    };
  mpfi_c_set(res,s_n_res);
};
    

void create_s_array(mpfi_c_ptr s,mpfi_c_t *s_array, unsigned int n)
{
  unsigned int i,del;
  double fact;

  for(i=0;i<n;i++)
    {
      mpfi_set(s_array[i*n+i]->im,s->im);
      mpfi_add_ui(s_array[i*n+i]->re,s->re,i);
    };
  //    s_array[i*n+i]=dcomplex(real(s)+i,imag(s));

  del=1;
  fact=2;
  while(del<n)
    {
      for(i=0;i<n-del;i++)
	{
	  mpfi_c_mul(s_array[i*n+i+del],s_array[i*n+i],
		     s_array[i*n+n+i+del]);
	  mpfi_div_d(s_array[i*n+i+del]->re,s_array[i*n+i+del]->re,fact);
	  mpfi_div_d(s_array[i*n+i+del]->im,s_array[i*n+i+del]->im,fact);
	};
      //	s_array[i*n+i+del]=s_array[i*n+i]*s_array[i*n+n+i+del]/fact;
      del++;
      fact++;
    };
};

void print_s_array(mpfi_c_t *s_array, unsigned int n)
{
  unsigned int i,j,ptr;

  ptr=0;
  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      {
	if(j>=i)
	  mpfi_c_print(s_array[ptr]);
	ptr++;
      };
};



void print_r_n_vec(mpfi_c_t *r_n_vec,unsigned int n)
{
  unsigned int i;

  for(i=0;i<n;i++)
    {
      printf("r_n[%d] is ",i);
      mpfi_c_print(r_n_vec[i]);
    };
  //cout << "r_n[" << i << "] is " << r_n_vec[i] << endl;
};


mpfi_c_t make_r_ns_tmp;

void make_r_ns (mpfi_c_t *s_array,
		unsigned int n,
		mpfi_c_t *this_r_n_vec,
		mpfi_c_t *last_r_n_vec,
		mpfi_t *delta_n)
{

  int i,j;
  //print_r_n_vec(last_r_n_vec,n);
  for(i=0;i<n;i++)
    {
      mpfi_c_set(this_r_n_vec[i],last_r_n_vec[i]);
      for(j=0;j<n-i-1;j++)
	//this_r_n_vec[i]+=s_array[i*(n-1)+j+i]*delta_n[j]*last_r_n_vec[j+i+1];
	{
	  mpfi_c_mul(make_r_ns_tmp,s_array[i*(n-1)+j+i],last_r_n_vec[j+i+1]);
	  mpfi_c_mul_i(make_r_ns_tmp,make_r_ns_tmp,delta_n[j]);
	  mpfi_c_add(this_r_n_vec[i],this_r_n_vec[i],make_r_ns_tmp);
	  //printf("delta_n[%d]= ",j);
	  //mpfi_print(delta_n[j]);
	  //printf("r_n_vec[%d]= ",i);
	  //mpfi_c_print(this_r_n_vec[i]);
	};
    };
  //print_r_n_vec(this_r_n_vec,n);

};

double limit;
int counter;

mpfr_t ts_tmp;

int too_small (mpfi_ptr x)
{
  double xd;

  mpfi_get_right(ts_tmp,x);   // x is +ve 
  xd=mpfr_get_d(ts_tmp,GMP_RNDU);
  if(xd<=limit)
    {
      counter++;
      return(TRUE);
    };
  return(FALSE);
};


mpfi_c_t c_term;
mpfi_t c_delta_n;


void calc_r_n(mpfi_c_ptr res,unsigned int r_n_terms,
		  mpfi_c_t *s_array, mpfi_ptr delta,
		  mpfi_c_t *last_r_n_vec)
{
  unsigned int i;

  mpfi_set(c_delta_n,delta);
  mpfi_c_set(res,last_r_n_vec[0]);
  //res=last_r_n_vec[0];

  // could decide how many terms to take according to delta 

  for(i=0;i<r_n_terms-1;)
  {
    mpfi_c_mul_i(c_term,s_array[i],c_delta_n);
    //term=s_array[i]*delta_n;
      
    i++;
    mpfi_c_mul(c_term,c_term,last_r_n_vec[i]);
    //term*=last_r_n_vec[i];
    mpfi_c_add(res,res,c_term);
    //res+=term;
    mpfi_mul(c_delta_n,c_delta_n,delta);

    //delta_n*=delta;

    //if(too_small(c_delta_n))
    //break;

    };
  //return(res);
};



void make_out_val(mpfi_c_ptr res,mpfi_c_ptr r_n,mpfi_c_t *n_s,
		  unsigned int a, unsigned int q, unsigned int N)
{
  unsigned int n;

  mpfi_c_mul(res,r_n,n_s[q]);
  for(n=0;n<N;n++)
    mpfi_c_add(res,res,n_s[n*q+a]);
  //res+=n_s[n*q+a];


  //return(r_n*n_s[q]+res);
};

void write_mpfi_c(mpfi_c_ptr z,FILE *outfile)
{
  mpfi_out_str(outfile,10,0,z->re);
  mpfi_out_str(outfile,10,0,z->im);
  fprintf(outfile,"\n");
}; 

mpfi_t one,h_x,h_frac,h_dq;
mpfi_c_t h_out_val,h_out_val1;

void calc_hurwitz1(unsigned int q_start,
		  unsigned int q_end,
		  mpfi_c_ptr s, 
		  unsigned int h_terms,
		  unsigned int N_terms,
		  mpfi_c_t *zeta_vals,
		  mpfi_ptr gap,
		  unsigned int no_gaps,
		  mpfi_c_t *s_array,
		  mpfi_c_t *r_n_vals,
		  FILE *out_file,
		  mpfi_c_ptr s0,
		  mpfi_c_t *n_s,
		  mpfi_t *gap_n)

{
  unsigned int i,j,gap_ptr,q,num,num_fracs,cmp;

  create_s_array(s,s_array,h_terms-1);
  //print_s_array(s_array,h_terms-1);
  //  gap_ptr=0;
  //  for(i=0;i<h_terms-1;i++)
  //    for(j=0;j<h_terms-1;j++)
  //      mpfi_c_print(s_array[gap_ptr++]);
 
  //printf("n_s within calc\n");
  //for(i=0;i<N_terms*q_end;i++)
  //  mpfi_c_print(n_s[i]); 
  //mpfi_c_set_ui(n_s[1],1,0);    // done once in main
  //mpfi_c_print(n_s[1]);
  mpfi_c_neg(s0,s);
  
  for(i=2;i<N_terms*q_end;i++)
    mpfi_ui_c_pow(n_s[i],i,s0);        //n_s[i]<-i^(-s)
    
  for(i=0;i<h_terms;i++)
    {
      s_n(r_n_vals[i],s0,one,N_terms);
      mpfi_c_sub(r_n_vals[i],zeta_vals[i],r_n_vals[i]);
      //r_n_vals[i]=zeta_vals[i]-s_n(s0,1.0,N);
      //cout << s0 << " " << zeta_vals[i] << " " << s_n(s0,1.0,N) << endl;
      //cout << "r_n[" << i << "] is " << last_big_r_n_vec[i] << endl;
      mpfi_c_dec_ui(s0,1);
      //mpfi_c_print(r_n_vals[i]);
    };

  mpfi_c_set(s0,s);
  
  for(i=1;i<no_gaps;i++)
    make_r_ns(s_array,h_terms,&r_n_vals[i*h_terms],
	      &r_n_vals[(i-1)*h_terms],gap_n);
  //print_r_n_vec(r_n_vals,h_terms*no_gaps);
  
  //open_file(0,s);
  //file_no=1;
  
  //out_val=make_out_val(r_n_vals[0],n_s,1,1,N);
  //write_dcomplex(out_val);
  //cout << "1 1" << endl;
   
  for(q=q_start;q<q_end;q++)
    {
      if((q&3)==2)
	continue;      // if q = 2 mod 4 then no primitive characters

      //printf("Working with conductor = %d\n",q);
      num_fracs=0;
      gap_ptr=no_gaps*h_terms;

      mpfi_mul_ui(h_x,gap,no_gaps);
      mpfi_sub(h_x,one,h_x);
      //x=1.0-gap*no_gaps;

      mpfi_set_ui(h_dq,q);
      //dq=(double) q;
      
      //if((q%QS_PER_FILE)==1)
      //{
      //flush_buffer();
      //out_file.close();
      //open_file(file_no++,s);
      //};
      
      for(num=1;num<q;)
	{
	  if(co_prime(num,q))
	    {
	      num_fracs++;
	      mpfi_set_ui(h_frac,num);
	      mpfi_div(h_frac,h_frac,h_dq);
	      //frac=(double) num/dq;

	      //printf("Working with numerator %d\n",num);
	      for(;mpfi_cmp(h_frac,h_x)>=0;)
	      
		//while(frac>=x)
		{
		  mpfi_add(h_x,h_x,gap);
		  //x=x+gap;
		  gap_ptr-=h_terms;
		};
	      mpfi_sub(h_frac,h_x,h_frac);
	      calc_r_n(h_out_val,h_terms,s_array,
		       h_frac,&r_n_vals[gap_ptr]);
	      //mpfi_c_print(h_out_val);
	      make_out_val(h_out_val1,h_out_val,n_s,num,q,N_terms);
	      //mpfi_c_print(h_out_val1);
	      write_mpfi_c(h_out_val1,out_file);
	    };
	  num++;
	  //	  if(num==q)
	  //{
	  //  printf("%d %d\n",num,num_fracs);
	  //  //cout << num << " " << num_fracs << endl;
	  //  break;
	  //};
	};
    };
  //close_file();
	  
};

void init_temps()
{
  mpfi_init(pow_tmp1);
  mpfi_init(pow_tlni);
  mpfi_init(pow_i_to_sigma);
  mpfi_init(pow_sin);
  mpfi_init(pow_cos);
  mpfi_init(s_n_tmp2);
  mpfi_c_init(s_n_res);
  mpfi_c_init(s_n_tmp1);
  mpfi_init(one);
  mpfi_set_ui(one,1);
  mpfi_c_init(make_r_ns_tmp);
  mpfi_init(h_x);
  mpfi_init(h_frac);
  mpfi_init(h_dq);
  mpfi_c_init(c_term);
  mpfi_init(c_delta_n);
  mpfi_c_init(h_out_val);
  mpfi_c_init(h_out_val1);
  mpfr_init(ts_tmp);
};


int main(int argc, char **argv)
{
  mpfi_c_t s,s0;
  mpfi_c_t *zeta_vals,*r_n_vals,*s_array,*n_s;
  int q_start,q_end,h_terms,N_terms;
  unsigned int max_h_terms,num_s_vals,den,num,no_gaps,i,j,ptr;
  mpfr_t delta;
  mpfi_t *gap_n;
  FILE *zeta_file,*out_file;

  double dgap;
  mpfi_t gap;
  mpfr_t gap_r,h_terms_r;
  char ch;

  mpfi_c_setup(PREC);
  
  if(argc!=8)
    print_usage;
  
  q_start=atoi(argv[1]);
  if(q_start<=2)
    print_usage;
  
  q_end=atoi(argv[2]);
  if(q_end<q_start)
    print_usage;
    
  h_terms=atoi(argv[3]);
  if(h_terms<=0)
    print_usage;
    
  N_terms=atoi(argv[4]);
  if(N_terms<=0)
    print_usage;
  
  sscanf(argv[5],"%lf",&dgap);   // atof don't work!
  if((dgap<=0.0)||(dgap>=1.0))
    print_usage;
  mpfi_init(gap);
  mpfi_set_d(gap,dgap);

  /* set limit variable = gap^h_terms */
  mpfr_init(gap_r);
  mpfr_init(h_terms_r);
  mpfr_set_d(gap_r,dgap,GMP_RNDD);
  mpfr_set_ui(h_terms_r,h_terms,GMP_RNDU);
  mpfr_pow(h_terms_r,gap_r,h_terms_r,GMP_RNDD);
  limit=mpfr_get_d(h_terms_r,GMP_RNDD);

  zeta_file=fopen(argv[6],"r");
  if(zeta_file==NULL)
    {
      printf("Unable to open Zeta File %s. Exiting.\n",argv[6]);
      return(SUCCESS);
    };

  out_file=fopen(argv[7],"w");
  if(out_file==NULL)
    {
      printf("Unable to open Out file %s. Exiting.\n",argv[7]);
      return(SUCCESS);
    };

  fscanf(zeta_file,"%d",&max_h_terms);
  if(h_terms>max_h_terms)
    {
      printf("Insufficent Zeta values in Zeta file. Exiting.\n");
      return(SUCCESS);
    };
  printf("Maximum h_terms supported by file: %d\n",max_h_terms);
  fscanf(zeta_file,"%d",&num_s_vals);
  fscanf(zeta_file,"%d",&den);
  fscanf(zeta_file,"%d",&num);
  mpfr_init(delta);
  mpfr_set_ui(delta,den,GMP_RNDN);
  mpfr_div_ui(delta,delta,num,GMP_RNDN);

  mpfi_c_init(s);
  mpfi_c_set_d(s,0.5,0.0);  // s=0.5 + 0.0i

  printf("s values increasing in steps of ");
  mpfr_print(delta);
  printf("Gap size set to: ");
  mpfi_print(gap);
  gap_n=malloc(h_terms*sizeof(mpfi_t));
  if(gap_n==NULL)
    {
      printf("Unable to allocate memory for gap_n. Exiting.\n");
      return(SUCCESS);
    };
  for(i=0;i<h_terms;i++)
    mpfi_init(gap_n[i]);
  mpfi_set(gap_n[0],gap);
  for(i=1;i<h_terms;i++)
    mpfi_mul(gap_n[i],gap_n[i-1],gap);
  zeta_vals=malloc(h_terms*sizeof(mpfi_c_t));
  if(zeta_vals==NULL)
    {
      printf("Unable to allocate memory for zeta_vals. Exiting.\n");
      return(SUCCESS);
    };


  for(i=0;i<h_terms;i++)
    mpfi_c_init(zeta_vals[i]);
  
  no_gaps=1.0/dgap; //ceil(1.0/dgap);   // why no work? 
  printf("So number of gaps is: %d\n",no_gaps);
  
  r_n_vals=malloc(no_gaps*h_terms*sizeof(mpfi_c_t));
  if(r_n_vals==NULL)
    {
      printf("Unable to allocate memory for r_n_vals. Exiting.\n");
      return(SUCCESS);
    };
  ptr=0;
  for(i=0;i<no_gaps;i++)
    for(j=0;j<h_terms;j++)
      mpfi_c_init(r_n_vals[ptr++]);
  
  s_array=malloc((h_terms-1)*(h_terms-1)*sizeof(mpfi_c_t));
  if(s_array==NULL)
    {
      printf("Unable to allocate memory for s_array. Exiting.\n");
      return(SUCCESS);
    };
  ptr=0;
  for(i=1;i<h_terms;i++)
    for(j=1;j<h_terms;j++)
      mpfi_c_init(s_array[ptr++]);

  n_s=malloc(N_terms*q_end*sizeof(mpfi_c_t));
  if(n_s==NULL)
    {
      printf("Unable to allocate memory for n_s. Exiting.\n");
      return(SUCCESS);
    };
  
  for(ptr=1;ptr<N_terms*q_end;ptr++)
    {
      mpfi_c_init(n_s[ptr]);   // don't use n_s[0], meaningless
    };

  mpfi_c_set_ui(n_s[1],1,0);

  //printf("n_s in main\n");
  //for(i=0;i<N_terms*q_end;i++)
  //  mpfi_c_print(n_s[i]);

  mpfi_c_init(s0);
  init_temps();


  
  for(i=0;i<=num_s_vals;i++)  /* there are num_s_vals+1 values! */
    {
      for(j=0;j<h_terms;j++)
	{
	  mpfi_inp_str(zeta_vals[j]->re,zeta_file,10);
	  mpfi_inp_str(zeta_vals[j]->im,zeta_file,10);
	};
      for(;j<max_h_terms;j++)
	{
	    mpfi_inp_str(s0->re,zeta_file,10);
	    mpfi_inp_str(s0->im,zeta_file,10);
	};

      calc_hurwitz1(q_start,q_end,s,h_terms,N_terms,zeta_vals,gap,no_gaps,
		    s_array,r_n_vals,out_file,s0,n_s,gap_n);




      mpfi_add_fr(s->re,s->re,delta);
    };
  //for(i=0;i<h_terms;i++)
  //mpfi_c_print(zeta_vals[i]);
  //  cout << "Entering calc_hurwitz1." << endl;
  //  if(!calc_hurwitz1(s,r_n_terms,N,zeta_vals,gap))
  //    fatal_error("Error running Hurwitz routine. Exiting.");


  printf("Too Small returned TRUE %d times.\n",counter);
  return(SUCCESS);
};

