/*
  File f.cpp
*/

#include "../includes/int_double12.0.h"
#include "../includes/fft_defs.h"
#include "../includes/im_s.h"
#include "f_defs.h"
#include "stdio.h"
#include "assert.h"

//#define PRINT
//#define PRINT1

#define INDEX_PER_FILE (100)

void print_usage()
{
  printf("Usage:- f (spec_file) (file_stub)\n");
  printf("spec_file - contains number of files followed by filename one per line.\n");
  printf("file_stub - has _q_index.dat appended for ouput\n");
  printf("Exiting.\n");
  exit(1);
}

// perform an in place FFT ,n a power of 2
// F[k]=sum_{j=0..N-1} F[j]e(-kj/N)
void fft(int_complex *x,unsigned int n,int_complex *w)
{
	unsigned int i,j,k,l;
	int_complex *p,*xend=x+n,ctemp;

	/* swap each element with one with bit-reversed index */
	for (i=0,l=n>>1;i<l;++i) {
		/* j = bit reversal of i */
		for (k=1,j=0;k<n;k<<=1) {
			j <<= 1;
			if (i & k)
				j |= 1;
		}
		if (i < j)
		{
			ctemp=x[i];
			x[i]=x[j];
			x[j]=ctemp;
		}
		else if (i > j)
		{
			ctemp=x[n-1-i];
			x[n-1-i]=x[n-1-j];
			x[n-1-j]=ctemp;
		}
		++i, j |= l;
		ctemp=x[i];
		x[i]=x[j];
		x[j]=ctemp;
	}

	for (k=1,l=n/2;k<n;k<<=1,l>>=1)
		for (p=x;p<xend;p+=k)
			for (j=0;j<n/2;j+=l,p++) {
				ctemp=p[k]*w[j];
				p[k]=p[0]-ctemp;
				p[0]+=ctemp;
			}
}

/* perform an in place inverse FFT */
void ifft(int_complex *x,unsigned int n,int_complex *w) {
	unsigned int i,l=n>>1;
	int_complex ctemp;

	fft(x,n,w);
	x[0]=x[0]/n;
	x[l]=x[l]/n;

	for (i=1;i<l;i++) {
		x[i]=x[i]/n;
		x[n-i]=x[n-i]/n;
		ctemp=x[i];
		x[i]=x[n-i];
		x[n-i]=ctemp;
	}
}

#define SMALL ((double) 1e-307)
#define ln_SMALL ((double) 706.89362355)  // exp(-ln_SMALL)<=SMALL


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
unsigned int N;
double B;

void setup(unsigned int q)
{
  assert(one_over_A<2*M_PI);
  c_small=int_complex(int_double(-SMALL,SMALL),int_double(-SMALL,SMALL));
  B=one_over_A*N;
  two_pi_by_B=d_two_pi/B;
  //printf("B=%f\n",B);
  A=d_one/one_over_A;
  pi_A=d_pi*A;
  two_pi_A=d_two_pi*A;
  WS[0]=c_one;
  int_double omega=d_two_pi/N;
  for(int i=1;i<N/2;i++)
    {
      sin_cos(omega*i,&WS[i].imag,&WS[i].real);
      //if((i%1000)==0)
      //printf("rel error WS[%d]=%d\n",i,-rel_error(WS[i].real));
      //print_int_complex_str("WS[]=",WS[i]);
    }
  d_log_pi=log(d_pi);
  d_log2=log(int_double(2));
  // this is log(Zeta(9/8)*(1/2)^(5/16))
  d_log_zeta_9_by_8=log(int_double(1575824,1575825)/1000000);
  exp_minus_pi_A=d_one-exp(-pi_A);
  d_pi_by_q=d_pi/q;
}

unsigned int make_N0(unsigned int q)
{
  return(1280384);
  return(51584); // t = 4030
  //return((unsigned int) ceil(10000000.0/q/one_over_A)*10.0+30);
}

FILE *make_file(unsigned int q, unsigned int index, const char *stub)
{
  char buff[1024];
  sprintf(buff,"%s_%d_%d.dat",stub,q,index);
  return(fopen(buff,"wb"));
} 

void save_data(unsigned int q,unsigned int index,const int_complex &omega,bool neg_one,bool real_p,unsigned int N0,FILE *outfile)
{
  fwrite(&q, sizeof(unsigned int),1,outfile);
  fwrite(&index, sizeof(unsigned int),1,outfile);
  fwrite(&omega, sizeof(int_complex),1,outfile);
  fwrite(&neg_one, sizeof(bool),1,outfile);
  fwrite(&real_p, sizeof(bool),1,outfile);
  fwrite(&N0, sizeof(unsigned int),1,outfile);
  for(int n=0;n<N0;n++)
    fwrite(&Fft_Vec[n].real, sizeof(int_double),1,outfile);
}

void save_data1(unsigned int index,const int_complex &omega,unsigned int N0,FILE *outfile)
{
  fwrite(&index, sizeof(unsigned int),1,outfile);
  fwrite(&omega, sizeof(int_complex),1,outfile);
  for(int n=0;n<N0;n++)
    fwrite(&Fft_Vec[n].real, sizeof(int_double),1,outfile);
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

int_double F_twiddle_e_err(const unsigned int m, const int_double &eta_pi_by_4, const int_double &logq)
{
  int_double t1=int_double(m)/A+B;
  int_double t2=int_double(m)/A-B;
  int_double beta1=beta_e(t1);
  assert(beta1.left>-eta_pi_by_4.right);
  int_double beta2=beta_e(t2);
  assert(beta2.left>-eta_pi_by_4.right);
  int_double res=E_e(t1,eta_pi_by_4,logq)/(d_one-exp(-(beta1-eta_pi_by_4)*B));
  res+=E_e(t2,eta_pi_by_4,logq)/(d_one-exp(-(beta2+eta_pi_by_4)*B));
  res.left=res.right;
  return(res);
}


int_double F_twiddle_o_err(const unsigned int m, const int_double &eta_pi_by_4, const int_double &logq)
{
  int_double t1=int_double(m)/A+B;
  int_double t2=int_double(m)/A-B;
  int_double beta1=beta_o(t1);
  assert(beta1.left>-eta_pi_by_4.right);
  int_double beta2=beta_o(t2);
  assert(beta2.left>-eta_pi_by_4.right);
  int_double res=E_o(t1,eta_pi_by_4,logq)/(d_one-exp(-(beta1-eta_pi_by_4)*B));
  res+=E_o(t2,eta_pi_by_4,logq)/(d_one-exp(-(beta2+eta_pi_by_4)*B));
  res.left=res.right;
  return(res);
}

/*
int_double F_hat_twiddle_o_err(const int_double &x, const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  //print_int_double_str("x=",x);
  int_double Xx=X_x(x,delta,logdelta,logq);
  if(Xx.left<=1.0)
    {
      print_int_double_str("X(x)=",Xx);
      print_int_double_str("x=",x);
      assert(Xx.left>1.0);
    }
  int_double res=exp(x*1.5-Xx-logq*0.75-logdelta*0.5+log(d_one+d_one/Xx/2.0)*1.5)/exp_minus_pi_A;
  res.right*=4.0;
  res.left=res.right;
  return(res);
}
*/

int_complex F_hat_twiddle_e_err( const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  int_double X_Pi_A=X_x(pi_A,delta,logdelta,logq);
  assert(X_Pi_A.left>1.0);
  int_double log_res=log(int_double(8))+two_pi_A-X_Pi_A+log(int_double(0.5)/X_Pi_A+1.0)-logq*0.25-logdelta*0.5-log(exp_minus_pi_A);
  if(log_res.left<-ln_SMALL)
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
  if(log_res.left<-ln_SMALL)
    return(c_small);
  log_res=exp(log_res);
  log_res.left=log_res.right;
  return(int_complex(log_res,log_res));
}
/*
int_double F_hat_twiddle_e_err(const int_double &x, const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  //print_int_double_str("x=",x);
  int_double Xx=X_x(x,delta,logdelta,logq);
  if(Xx.left<=1.0)
    {
      print_int_double_str("X(x)=",Xx);
      print_int_double_str("x=",x);
      assert(Xx.left>1.0);
    }
  int_double res=exp(x*0.5-Xx-logq*0.25-logdelta*0.5)*(d_one+d_one/Xx/2.0)/exp_minus_pi_A;
  res.right*=4.0;
  res.left=res.right;
  return(res);
}
*/
inline void do_F(const int_double &delta, const int_double &logdelta, const int_double &logq, unsigned int N0, const int_double &eta_pi_by_4,
		      const int_complex &F_hat_t_err, int_double *F_twiddle_errs)
{
  int i;
  int_double n_two_pi_by_B;

  Fft_Vec[0]+=F_hat_t_err;
  Fft_Vec[N/2]+=F_hat_t_err;

  for(int n=1;n<N/2;n++)
    {
      Fft_Vec[N-n]=conj(Fft_Vec[n]);
    }

  //for(int n=0;n<N0;n+=1000)
  //  print_int_complex_str("Pre-FFT",Fft_Vec[n]);

  //printf("Running fft.\n");      
  fft(Fft_Vec,N,WS);
  for(i=0;i<N0;i++)
    Fft_Vec[i]*=two_pi_by_B;
  //for(int n=0;n<N0;n+=1000)
  //  print_int_complex_str("Post-FFT",Fft_Vec[n]);
  

  //printf("fft finished.\n");

  double x;
  for(x=0.0,i=0;i<N0;i++,x+=one_over_A)
    {
      Fft_Vec[i].real+=F_twiddle_errs[i];
      Fft_Vec[i].imag+=F_twiddle_errs[i];
      Fft_Vec[i].real*=exp(x*(d_pi/4.0-eta_pi_by_4));
      if(!contains_zero(Fft_Vec[i].imag))
	{
	  printf("Fft_Vec[%d] not real.\n",i);
	  print_int_complex_str("=",Fft_Vec[i]);
	  printf("Exiting./n");
	  //print_int_complex_str("F(t)=",L(q,int_double(i)*one_over_A,eta,chis,epsilon));
	  exit(1);
	}
    }
  //for(int n=0;n<N0;n+=1000)
  //  print_int_complex_str("Post-errs",Fft_Vec[n]);
  //exit(0);
}


inline void check_uint(FILE **infiles, unsigned int num_files, unsigned int val)
{
  unsigned int i;
  for(int n=1;n<num_files;n++)
    {
      assert(fread(&i,sizeof(unsigned int),1,infiles[n]));
      assert(i==val);
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


inline void F(FILE **infiles, unsigned int num_files, const char *file_stub)
{
  unsigned int q,num_prims,*num_ss,*n0s,N0,index,N_ptr,num_indices=0,file_count=0;
  double eta_o,eta_e;
  bool neg_one,real_p;
  FILE* outfile;
  int_complex omega,omega1;
  int_double delta_o=d_pi_2*(1.0-eta_o);
  int_double logdelta_o=log(delta_o);
  int_double delta_e=d_pi_2*(1.0-eta_e);
  int_double logdelta_e=log(delta_e);
  int_double logq;
  int_complex F_hat_t_o_err,F_hat_t_e_err;
  int_double *F_twiddle_o_errs,*F_twiddle_e_errs;
  int_double n_two_pi_by_B;

  printf("In F.\n");
  assert(num_ss=(unsigned int *) malloc(sizeof(unsigned int)*num_files));
  assert(n0s=(unsigned int *) malloc(sizeof(unsigned int)*num_files));
  assert(fread(&q,sizeof(unsigned int),1,infiles[0]));
  check_uint(infiles,num_files,q);
  printf("q=%d\n",q);
  N0=make_N0(q);
  logq=log(int_double(q));
  assert(fread(&num_prims,sizeof(unsigned int),1,infiles[0]));
  check_uint(infiles,num_files,num_prims);
  printf("num_prims=%d\n",num_prims);
  assert(fread(&N,sizeof(unsigned int),1,infiles[0]));
  check_uint(infiles,num_files,N);
  printf("N=%d N0=%d\n",N,N0);
  assert(N>8*N0);
  assert(Fft_Vec=(int_complex *) _aligned_malloc(sizeof(int_complex)*N,16));
  assert(WS=(int_complex *) _aligned_malloc(sizeof(int_complex)*N/2,16));
  assert(F_twiddle_o_errs=(int_double *) _aligned_malloc(sizeof(int_double)*N0,16));
  assert(F_twiddle_e_errs=(int_double *) _aligned_malloc(sizeof(int_double)*N0,16));

  assert(fread(&eta_o,sizeof(double),1,infiles[0]));
  printf("eta_o=%f\n",eta_o);
  check_double(infiles,num_files,eta_o);
  assert(fread(&eta_e,sizeof(double),1,infiles[0]));
  printf("eta_e=%f\n",eta_e);
  check_double(infiles,num_files,eta_e);
  assert(fread(&num_ss[0],sizeof(unsigned int),1,infiles[0]));
  assert(fread(&n0s[0],sizeof(unsigned int),1,infiles[0]));
  assert(n0s[0]==0);
  for(unsigned int n=1;n<num_files;n++)
    {
      assert(fread(&num_ss[n],sizeof(unsigned int),1,infiles[n]));
      assert(fread(&n0s[n],sizeof(unsigned int),1,infiles[n]));
      assert(n0s[n]==(n0s[n-1]+num_ss[n-1]));
    }
  assert((n0s[num_files-1]+num_ss[num_files-1])==N/2+1);
  outfile=make_file(q,file_count++,file_stub);
  setup(q);
  printf("Making F_hat_twiddle errors..\n");
  F_hat_t_e_err=F_hat_twiddle_e_err(delta_e,logdelta_e,logq);
  F_hat_t_o_err=F_hat_twiddle_o_err(delta_o,logdelta_o,logq);
  print_int_complex_str("F_hat_twiddle_o_err=",F_hat_t_o_err);
  print_int_complex_str("F_hat_twiddle_e_err=",F_hat_t_e_err);
  for(unsigned int prim=0;prim<num_prims;prim++)
    {
      assert(fread(&real_p,sizeof(bool),1,infiles[0]));
      check_bool(infiles,num_files,real_p);
      assert(fread(&index,sizeof(unsigned int),1,infiles[0]));
      printf("index=%d\n",index);
      check_uint(infiles,num_files,index);
      if(num_indices==INDEX_PER_FILE)
	{
	  fclose(outfile);
	  outfile=make_file(q,file_count++,file_stub);
	  num_indices=0;
	}
      else
	num_indices++;
      assert(fread(&omega,sizeof(int_complex),1,infiles[0]));
      //print_int_complex_str("omega=",omega);
      check_int_complex(infiles,num_files,omega);
      assert(fread(&neg_one,sizeof(bool),1,infiles[0]));
      check_bool(infiles,num_files,neg_one);
      N_ptr=0;
      for(int n=0;n<num_files;n++)
	{
	  assert(fread(&Fft_Vec[N_ptr],sizeof(int_complex),num_ss[n],infiles[n]));
	  N_ptr+=num_ss[n];
	}
      for (int n=1;n<N/2;n<<=1)
	{
	  printf("q:%d index :%d f_hat(%d)=",q,index,n);
	  print_int_complex_str("",Fft_Vec[n]);
	}
      if(neg_one)
	do_F(delta_o,logdelta_o,logq,N0,d_pi*(eta_o/4.0),F_hat_t_o_err,F_twiddle_o_errs);
      else
	do_F(delta_e,logdelta_e,logq,N0,d_pi*(eta_e/4.0),F_hat_t_e_err,F_twiddle_e_errs);
      //printf("Saving index = %d.\n",index);
      save_data(q,index,omega,neg_one,real_p,N0,outfile);
      if(real_p)
	  continue;
      prim++;
      assert(fread(&index,sizeof(unsigned int),1,infiles[0]));
      //printf("index1=%d\n",index);
      check_uint(infiles,num_files,index);
      assert(fread(&omega1,sizeof(int_complex),1,infiles[0]));
      //print_int_complex_str("omega1=",omega1);
      check_int_complex(infiles,num_files,omega1);
      assert(contains_zero(omega-conj(omega1)));
      N_ptr=0;
      for(int n=0;n<num_files;n++)
	{
	  assert(fread(&Fft_Vec[N_ptr],sizeof(int_complex),num_ss[n],infiles[n]));
	  N_ptr+=num_ss[n];
	}
      for (int n=1;n<N/2;n<<=1)
	{
	  printf("q:%d index :%d f_hat(%d)=",q,index,n);
	  print_int_complex_str("",Fft_Vec[n]);
	}

      if(neg_one)
	do_F(delta_o,logdelta_o,logq,N0,d_pi*(eta_o/4.0),F_hat_t_o_err,F_twiddle_o_errs);
      else
	do_F(delta_e,logdelta_e,logq,N0,d_pi*(eta_e/4.0),F_hat_t_e_err,F_twiddle_e_errs);
      //printf("Saving composite index %d\n",index);
      save_data1(index,omega1,N0,outfile);
    }
}

int main(int argc, char **argv)
{
  FILE **infiles,*infile;
  unsigned int num_files;
  char fname[1024];
  _fpu_rndd();
  //printf("argc=%d\n",argc);
  if(argc!=3)
    print_usage();
  assert(infile=fopen(argv[1],"r"));
  assert(fscanf(infile,"%d\n",&num_files));
  assert(infiles=(FILE **) malloc(sizeof(FILE *)*num_files));
  //printf("No. files=%d.\n",num_files);
  for(int n=0;n<num_files;n++)
    {
      assert(fscanf(infile,"%s\n",fname));
      //printf("File %d=%s\n",n,fname);
      assert(infiles[n]=fopen(fname,"rb"));
    }
  fclose(infile);
  F(infiles,num_files,argv[2]);
  


  return(0);
}

