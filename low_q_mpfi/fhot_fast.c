#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"
#include "../includes/mpfi.h"
#include "../includes/mpfi_io.h"
#include "../includes/mpfi_c.h"
#include "./f_defs.h"



void print_usage(char *command)
{
  printf("Usage:- %s prec q num_files num_file outfile.\n",command);
  exit(1);
}

unsigned long int calc_N(unsigned long int q)
{
  unsigned long int res=1;
  double target=(ceil(QT/(10.0*q))*10.0+30.0)/one_over_A*7.5;
  //printf("Trying to get to N0=%ld\n",(int) target);
  while(res<target) res<<=1;
  //printf("FFT length set to %ld\n",res);
  return(res);
}

// results in attenuation of about 1/10 at height QT/q+30
double calc_eta(unsigned int q)
{
  return(1-7.03/(QT/q+30.0));
}

inline int gcd (unsigned int a, unsigned int b)
/* Euclid algorithm gcd */
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

inline unsigned int min(unsigned int x, unsigned int y)
{
  if(x<y) return(x);
  return(y);
}

// compute in double how large we need to take
// M to get a reasonable error bound
//
inline unsigned int calc_M(double x, double log2byq, double cosetaetc)
{
  double b=log2byq+x*1.5;
  //printf("x=%e\nlog(2)-log(q)*0.75=%e\ncos(eta)=%e\n",x,log2byq,cosetaetc);
  //printf("b=%e\n",b);
  //printf("lambda=%e\n",exp(2*x)*cosetaetc);
  if(b>TARGET_ACC)
    return(1);
  return((unsigned int) sqrt((TARGET_ACC-b)/exp(2*x)/cosetaetc)+1);
}

// set w_vec[i]=w^[2^i]
//
void set_w_vec(mpfi_c_t *w_vec, long unsigned int log2M, mpfi_c_ptr w)
{
  long unsigned int ptr=1;
  mpfi_c_set(w_vec[0],w);
  while(ptr<=log2M)
    {
      mpfi_c_sqr(w_vec[ptr],w_vec[ptr-1]);
      //mpfi_c_print_str("w_vec[]=",w_vec[ptr]);
      ptr++;
    }
}

// w_vec is w,w^2,w^4...
// where w=exp(- pi exp(2u(x))/q)
void calc_m_term(mpfi_c_ptr res, long unsigned int m, mpfi_c_t *w_vec, mpfi_c_ptr u)
{
  mpfi_c_set_ui(res,m,0);
  long unsigned int bit=0,m2=m*m;
  while(m2>0)
    {
      if(m2&1)
	mpfi_c_mul(res,res,w_vec[bit]);
      bit++;
      m2>>=1;
    }
  mpfi_c_mul(res,res,u);
}

int main(int argc, char **argv)
{
  FILE *outfile;

  mpfi_c_t u,v,w,*w_vec,*res_vec,u_x0,res;

  mpfi_t x,two_pi_by_B,vdel,udel,tmp;

  long unsigned int q,phi_q,M,m,n,n0,n1,N,num_files,num_file,log2M,ptr;

  double B,eta,log2byq,cosetaetc;

  printf("Command line:-");
  for(n=0;n<argc;n++)
    printf(" %s",argv[n]);
  printf("\n");
  if(argc!=6)
    print_usage(argv[0]);
  mpfi_c_setup(atol(argv[1]));

  q=atol(argv[2]);

  for(n=1,phi_q=0;n<q;n++)
    if(co_prime(n,q)) phi_q++;

  N=calc_N(q);
  eta=calc_eta(q);

  num_files=atol(argv[3]);
  num_file=atol(argv[4]);

  n0=N/2/num_files*num_file;
  if(num_file==num_files-1)
    n1=N/2+1;
  else
    n1=n0+N/2/num_files;
  //printf("N=%lu\nn0=%lu n1=%lu\n",N,n0,n1);

  B=N*one_over_A;
  outfile=fopen(argv[5],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[5]);
      exit(1);
    }
  fwrite(&q,sizeof(unsigned int),1,outfile);
  fwrite(&eta,sizeof(double),1,outfile);
  fwrite(&n0,sizeof(unsigned int),1,outfile);
  fwrite(&n1,sizeof(unsigned int),1,outfile);
  //printf("sizeof(unsigned int)=%d bytes\n",sizeof(unsigned int));
  fwrite(&N,sizeof(unsigned int),1,outfile);
  //printf("q=%ld\neta=%e\nn0=%ld\nn1=%ld\nN=%ld\n",q,eta,n0,n1,N);

  mpfi_init(two_pi_by_B);
  mpfi_div_d(two_pi_by_B,mpfi_2_pi,B);

  mpfi_init(x);
  mpfi_mul_ui(x,two_pi_by_B,n0);

  //mpfi_print_str("x0=",x);

  mpfi_c_init(u_x0);
  mpfi_set(u_x0->re,x);
  mpfi_mul_d(u_x0->im,mpfi_pi_by_4,eta);

  mpfi_c_print_str("u(x0)=",u_x0);
  // v is -Pi*exp(2u(x))/q
  mpfi_c_init(v);
  mpfi_c_mul_ui(v,u_x0,2);
  mpfi_c_exp(v,v);
  mpfi_c_div_ui(v,v,q);
  mpfi_c_mul_i(v,v,mpfi_pi);
  mpfi_c_neg(v,v);
  mpfi_c_print_str("v=",v);

  mpfi_init(vdel);
  mpfi_mul_2ui(vdel,two_pi_by_B,1);
  mpfi_exp(vdel,vdel);

  mpfi_c_init(w);
  mpfi_c_exp(w,v);
  mpfi_c_print_str("w=",w);

  mpfi_c_init(u);
  mpfi_c_mul_d(u,u_x0,1.5);
  mpfi_c_exp(u,u);
  mpfi_c_mul_ui(u,u,2);
  mpfi_init(tmp);
  mpfi_set_ui(tmp,q);
  mpfi_log(tmp,tmp);
  mpfi_mul_d(tmp,tmp,-0.75);
  mpfi_exp(tmp,tmp);
  mpfi_c_mul_i(u,u,tmp);
  // u=2*exp(1.5*u(x))*q^(-3/4)

  mpfi_init(udel);
  mpfi_mul_d(udel,two_pi_by_B,1.5);

  cosetaetc=cos(M_PI*eta/2.0)*M_PI/q;
  log2byq=log((double) 2.0)-log((double) q)*0.75; // log((2/q)^0.75)

  M=calc_M(n0*2.0*M_PI/B,log2byq,cosetaetc);
  if(M==1)
    log2M=1;
  else
    log2M=floor(log((double) M*M)/log((double) 2.0));

  fwrite(&M,sizeof(unsigned int),1,outfile);

  printf("M=%lu\nlog_2(M^2)=%lu\n",M,log2M);
  w_vec=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*(log2M+1));
  res_vec=(mpfi_c_t *)malloc(sizeof(mpfi_c_t)*phi_q);

  for(n=0;n<=log2M;n++)
    mpfi_c_init(w_vec[n]);

  for(n=0;n<min(M,phi_q);n++)
    {
      mpfi_c_init(res_vec[n]);
      mpfi_c_set_ui(res_vec[n],0,0);
    }

  set_w_vec(w_vec,log2M,w);
  //printf("w_vec initialised.\n");
  mpfi_c_init(res);

  //printf("Reached line %lu\n",__LINE__);


  for(m=1,ptr=0;m<=M;m++)
    if(co_prime(m,q))
      {
	calc_m_term(res,m,w_vec,u);
	//mpfi_c_print_str("calc_m_term returned ",res);
	//printf("m*m=%lu ptr=%lu\n",m*m,ptr);
	mpfi_c_add(res_vec[ptr],res_vec[ptr],res);
	ptr++;
	if(ptr==phi_q)
	  ptr=0;
      }

  //printf("Reached line %lu\n",__LINE__);

  for(ptr=0;ptr<min(M,phi_q);ptr++)
    {
      mpfi_write_bin(outfile,res_vec[ptr]->re);
      mpfi_write_bin(outfile,res_vec[ptr]->im);
    }

  n=n0+1;
  //printf("Reached line %lu\n",__LINE__);
  printf("going to do %lu terms\n",n1-n0);
  while(n<n1)
    {
      if((n%1000)==0)
	printf("n=%lu\n",n);
      mpfi_c_mul_i(v,v,vdel);
      mpfi_c_mul_i(u,u,udel);
      mpfi_c_exp(w,v);

      M=calc_M(n*2.0*M_PI/B,log2byq,cosetaetc);
      if(M==1)
	log2M=1;
      else
	log2M=floor(log((double) M*M)/log((double) 2.0));
      fwrite(&M,sizeof(unsigned int),1,outfile);

      set_w_vec(w_vec,log2M,w);

      for(ptr=0;ptr<min(M,phi_q);ptr++)
	mpfi_c_set_ui(res_vec[ptr],0,0);

      for(m=1,ptr=0;m<=M;m++)
	if(co_prime(m,q))
	  {
	    calc_m_term(res,m,w_vec,u);
	    mpfi_c_add(res_vec[ptr],res_vec[ptr],res);
	    ptr++;
	    if(ptr==phi_q)
	      ptr=0;
	  }

      for(ptr=0;ptr<min(M,phi_q);ptr++)
	{
	  mpfi_write_bin(outfile,res_vec[ptr]->re);
	  mpfi_write_bin(outfile,res_vec[ptr]->im);
	}
      n++;
    }
  mpfi_c_print_str("last result saved",res_vec[ptr-1]);

  return(0);
}
