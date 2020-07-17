/*
  File f_hat_split.cpp
  reads output of f_hat_even/odd
  adds F_hat_twiddle error
  splits into files ready for disk based FFT
*/

#include "../includes/int_double12.0.h"
#include "../includes/fft_defs.h"
#include "../includes/im_s.h"
#include "f_defs.h"
#include "stdio.h"
#include "assert.h"

//#define PRINT
//#define PRINT1

void print_usage()
{
  printf("Usage:- f <spec_file> <file-stub> <no_splits> <even_p>\n");
  printf("spec_file - contains number of files followed by filename one per line.\n");
  printf("file_stub - has [0..<no_splits-1>] appended for ouput\n");
  printf("no_splits - power of two, number of ways to split\n");
  printf("even_p=1 if even, 0 if odd.\n");
  printf("Exiting.\n");
  exit(1);
}



int_double A;
int_double two_pi_A;
int_double two_pi_by_B;
int_complex *WS;
int_complex *Fft_Vec;
int_double d_log_pi;
int_double d_log2;
int_double d_log_zeta_9_by_8;
int_double d_small;
int_complex c_small;
int_double exp_minus_pi_A;
int_double d_pi_by_q;
int_double pi_A;
long unsigned int N;
double B;

void setup(unsigned int q)
{
  assert(one_over_A<2*M_PI);
  d_small=exp(int_double(LN_SMALL));
  d_small.left=d_small.right;
  c_small=int_complex(d_small,d_small);
  B=one_over_A*N;
  printf("B=%f\n",B);
  two_pi_by_B=d_two_pi/B;
  A=d_one/one_over_A;
  pi_A=d_pi*A;
  two_pi_A=d_two_pi*A;
  d_log_pi=log(d_pi);
  d_log2=log(int_double(2));
  // this is log(Zeta(9/8)*(1/2)^(5/16))
  d_log_zeta_9_by_8=log(int_double(1575824,1575825)/1000000);
  exp_minus_pi_A=d_one-exp(-pi_A);
  d_pi_by_q=d_pi/q;
}

unsigned int make_N0(unsigned int q)
{
  //return(13184); // t=10060
  return((unsigned int) (ceil(QT/(q*10.0))*10.0+30)/one_over_A);
}

FILE *make_file(unsigned int n, const char *stub)
{
  char buff[1024];
  sprintf(buff,"%s%d.dat",stub,n);
  return(fopen(buff,"wb"));
}

inline bool power_2_p(unsigned int i)
{
  while((i&1)==0)
    i>>=1;
  return(i==1);
}

inline int_double X_x(const int_double &x, const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  return(exp(x*2+logdelta+d_log_pi-logq-delta));
}

inline int_double log_X_x(const int_double &x, const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  return(x*2+logdelta+d_log_pi-logq-delta);
}

int_double E_e (const int_double &t, const int_double &eta_pi_by_4, const int_double &logq)
{
  int_complex lng=lngamma(int_complex(int_double(0.25),t/2));
  int_double res=exp(lng.real+eta_pi_by_4*t+logq*0.3125-d_log_pi*0.5625+d_log_zeta_9_by_8+log(t*t+9.0/4.0)*0.15625);
  //print_int_double_str("E_o called with t=",t);
  //print_int_double_str("Returning",res);
  return(res);
}

int_double beta_e(const int_double &t)
{
  int_double res=d_pi/4-atan2(d_one,abs(t)*2)*0.5-d_one/d_pi/d_pi/abs(t*2-1.0/4.0)*4;
  return(res);
}


int_double E_o (const int_double &t, const int_double &eta_pi_by_4, const int_double &logq)
{
  int_complex lng=lngamma(int_complex(int_double(0.75),t/2));
  int_double res=exp(lng.real+eta_pi_by_4*t+logq*0.3125-d_log_pi*1.0625+d_log_zeta_9_by_8+log(t*t+9.0/4.0)*0.15625);
  //print_int_double_str("E_o called with t=",t);
  //print_int_double_str("Returning",res);
  return(res);
}

int_double beta_o(const int_double &t)
{
  int_double res=d_pi/4-atan2(d_one,abs(t)*2)*1.5-d_one/d_pi/d_pi/abs(t*2-9.0/4.0)*4;
  return(res);
}


int_complex F_hat_twiddle_e_err( const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  int_double X_Pi_A=X_x(pi_A,delta,logdelta,logq);
  assert(X_Pi_A.left>1.0);
  int_double log_res=log(int_double(8))+two_pi_A-X_Pi_A+log(int_double(0.5)/X_Pi_A+1.0)-logq*0.25-logdelta*0.5-log(exp_minus_pi_A);
  if(log_res.left<-LN_SMALL)
    return(c_small);
  log_res=exp(log_res);
  log_res.left=log_res.right;
  return(int_complex(log_res,log_res));
}


int_complex F_hat_twiddle_o_err( const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  int_double X_Pi_A=X_x(pi_A,delta,logdelta,logq);
  assert(X_Pi_A.left>1.0);
  int_double log_res=log(int_double(8))+two_pi_A*3.0-X_Pi_A+log(int_double(0.5)/X_Pi_A+1.0)*1.5-logq*0.75-logdelta*0.5-log(exp_minus_pi_A);
  if(log_res.right>-LN_SMALL)
    return(c_small);
  log_res=exp(log_res);
  log_res.left=log_res.right;
  return(int_complex(log_res,log_res));
}


inline void check_uint(FILE **infiles, unsigned int num_files, unsigned int val)
{
  unsigned int i;
  for(int n=1;n<num_files;n++)
    {
      assert(fread(&i,sizeof(unsigned int),1,infiles[n]));
      assert(i==val);
      /*
      if(!fread(&i,sizeof(unsigned int),1,infiles[n]))
	{
	  printf("Failed to read int from file %d.\nExiting.\n",n);
	  exit(0);
	}
      if(i!=val)
	{
	  printf("Was looking for %d, got %d, in file %d.\nExiting.\n",val,i,n);
	  exit(0);
	}
      */
    }
}
inline void check_luint(FILE **infiles, unsigned int num_files, unsigned long int val)
{
  unsigned long int i;
  for(int n=1;n<num_files;n++)
    {
      assert(fread(&i,sizeof(unsigned long int),1,infiles[n]));
      assert(i==val);
      /*
      if(!fread(&i,sizeof(unsigned int),1,infiles[n]))
	{
	  printf("Failed to read int from file %d.\nExiting.\n",n);
	  exit(0);
	}
      if(i!=val)
	{
	  printf("Was looking for %d, got %d, in file %d.\nExiting.\n",val,i,n);
	  exit(0);
	}
      */
    }
}

inline void check_double(FILE **infiles, unsigned int num_files, double val)
{
  double i;
  for(int n=1;n<num_files;n++)
    {
      assert(fread(&i,sizeof(double),1,infiles[n]));
      assert(i==val);
    }
}

inline void check_bool(FILE **infiles, unsigned int num_files, bool val)
{
  bool i;
  for(int n=1;n<num_files;n++)
    {
      //printf("checking bool %d\n",n);
      assert(fread(&i,sizeof(bool),1,infiles[n]));
      assert(i==val);
    }
}

inline void check_int_complex(FILE **infiles, unsigned int num_files, const int_complex &val)
{
  int_complex i;
  for(int n=1;n<num_files;n++)
    {
      assert(fread(&i,sizeof(int_complex),1,infiles[n]));
      assert(contains_zero(i-val));
    }
}

void wuint(FILE ** outfiles, unsigned int num_files, unsigned int i)
{
  for(int n=0;n<num_files;n++)
    assert(fwrite(&i,sizeof(unsigned int),1,outfiles[n]));
}
void wluint(FILE ** outfiles, unsigned int num_files, unsigned long int i)
{
  for(int n=0;n<num_files;n++)
    assert(fwrite(&i,sizeof(unsigned long int),1,outfiles[n]));
}
void wdouble(FILE ** outfiles, unsigned int num_files, double i)
{
  for(int n=0;n<num_files;n++)
    assert(fwrite(&i,sizeof(double),1,outfiles[n]));
}

void wbool(FILE ** outfiles, unsigned int num_files, bool i)
{
  for(int n=0;n<num_files;n++)
    assert(fwrite(&i,sizeof(bool),1,outfiles[n]));
}

void wint_complex(FILE ** outfiles, unsigned int num_files, int_complex &i)
{
  for(int n=0;n<num_files;n++)
    assert(fwrite(&i,sizeof(int_complex),1,outfiles[n]));
}


int main(int argc, char **argv)
{
  FILE **infiles,*infile,**outfiles;
  unsigned int num_files,*num_ss;
  char fname[1024];
  int no_splits;
  unsigned int q,N0,num_prims,index;
  double eta;
  int_double logq,delta;
  int_complex z,F_hat_t_err,omega,omega1;
  bool neg_one,real_p;
  _fpu_rndd();
  if(argc!=5)
    print_usage();
  assert(infile=fopen(argv[1],"r"));
  assert(fscanf(infile,"%d\n",&num_files));
  printf("in f_hat_split with num_files=%d\n",num_files);
  assert(infiles=(FILE **) malloc(sizeof(FILE *)*num_files));
  assert(num_ss=(unsigned int *) malloc(sizeof(unsigned int)*num_files));
  //assert(n0s=(unsigned int *) malloc(sizeof(unsigned int)*num_files));
  //printf("No. files=%d.\n",num_files);
  for(int n=0;n<num_files;n++)
    {
      assert(fscanf(infile,"%s\n",fname));
      //printf("File %d=%s\n",n,fname);
      infiles[n]=fopen(fname,"rb");
      if(!infiles[n])
	{
	  printf("Error opening file %s\nExiting.\n",fname);
	  perror("");
	  exit(0);
	}
    }
  fclose(infile);
  no_splits=atoi(argv[3]);
  if((no_splits<=0)||(!power_2_p(no_splits)))
    print_usage();
  bool even_p=(atoi(argv[4])!=0);
  assert(outfiles=(FILE **) malloc(sizeof(FILE *)*no_splits));
  for(int n=0;n<no_splits;n++)
    assert(outfiles[n]=make_file(n,argv[2]));

  assert(fread(&q,sizeof(unsigned int),1,infiles[0]));
  printf("q=%d\n",q);
  check_uint(infiles,num_files,q);
  //printf("q=%d\n");
  wuint(outfiles,no_splits,q);
  N0=make_N0(q);
  printf("N0=%d\n",N0);
  logq=log(int_double(q));
  //print_int_double_str("log(q)=",logq);
  assert(fread(&num_prims,sizeof(unsigned int),1,infiles[0]));
  check_uint(infiles,num_files,num_prims);
  wuint(outfiles,no_splits,num_prims);
  assert(fread(&N,sizeof(unsigned long int),1,infiles[0]));
  check_luint(infiles,num_files,N);
  wluint(outfiles,no_splits,N);
  printf("N=%ld\n",N);
  setup(q);
  assert(fread(&eta,sizeof(double),1,infiles[0]));
  check_double(infiles,num_files,eta);
  printf("eta=%e\n",eta);
  wdouble(outfiles,no_splits,eta);
  delta=d_pi_2*(1.0-eta);
  if(even_p)
    F_hat_t_err=F_hat_twiddle_e_err(delta,log(delta),logq);
  else
    F_hat_t_err=F_hat_twiddle_o_err(delta,log(delta),logq);
  print_int_complex_str("F_hat_twiddle error = ",F_hat_t_err);

  // file zero
  //assert(fread(&num_ss[0],sizeof(unsigned int),1,infiles[0]));
  //assert(fread(&n0s[0],sizeof(unsigned int),1,infiles[0]));
  //assert(n0s[0]==0);
  for(int f=0;f<num_files;f++)
    assert(fread(&num_ss[f],sizeof(unsigned int),1,infiles[f]));
  //assert(fread(&n0s[f],sizeof(unsigned int),1,infiles[f]));
  //printf("num_ss[%d]=%d n0s[%d]=%d\n",f,num_ss[f],f,n0s[f]);
  //assert(num_ss[f]+n0s[f-1]==n0s[f]);
  //assert(n0s[num_files-1]+num_ss[num_files-1]==N/2+1);

  for(unsigned int prim=0;prim<num_prims;prim++)
    {
      //printf("Processing character %d\n",prim);
      assert(fread(&neg_one,sizeof(bool),1,infiles[0]));
      check_bool(infiles,num_files,neg_one);
      wbool(outfiles,no_splits,neg_one);
      /*
      if(neg_one)
	printf("%d negative character.\n",prim);
      else
	printf("%d positive character.\n",prim);
      */
      assert(fread(&real_p,sizeof(bool),1,infiles[0]));
      check_bool(infiles,num_files,real_p);
      wbool(outfiles,no_splits,real_p);
      if(!real_p)
	prim++;
      if(neg_one==even_p)
	continue;
      assert(fread(&index,sizeof(unsigned int),1,infiles[0]));
      check_uint(infiles,num_files,index);
      wuint(outfiles,no_splits,index);
      //printf("%d index = %d\n",prim,index);
      assert(fread(&omega,sizeof(int_complex),1,infiles[0]));
      check_int_complex(infiles,num_files,omega);
      wint_complex(outfiles,no_splits,omega);
      int n_ptr=0;
      for(int f=0;f<num_files;f++)
	for(int n=0;n<num_ss[f];n++)
	  {
	    assert(fread(&z,sizeof(int_complex),1,infiles[f]));
	    //if((f==0)&&(n==0))
	    //print_int_complex_str("first z = ",z);
	    z+=F_hat_t_err;
	    assert(fwrite(&z,sizeof(int_complex),1,outfiles[n_ptr]));
	    n_ptr++;
	    if(n_ptr==no_splits)
	      n_ptr=0;
	  }
      if(real_p)
	continue;
      assert(fread(&index,sizeof(unsigned int),1,infiles[0]));
      check_uint(infiles,num_files,index);
      wuint(outfiles,no_splits,index);
      //printf("index = %d\n",index);
      assert(fread(&omega1,sizeof(int_complex),1,infiles[0]));
      check_int_complex(infiles,num_files,omega1);
      wint_complex(outfiles,no_splits,omega1);
      assert(contains_zero(omega-conj(omega1)));
      n_ptr=0;
      for(int f=0;f<num_files;f++)
	for(int n=0;n<num_ss[f];n++)
	  {
	    assert(fread(&z,sizeof(int_complex),1,infiles[f]));
	    z+=F_hat_t_err;
	    assert(fwrite(&z,sizeof(int_complex),1,outfiles[n_ptr]));
	    n_ptr++;
	    if(n_ptr==no_splits)
	      n_ptr=0;
	  }
    }
  return(0);
}

