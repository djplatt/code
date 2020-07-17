/*
  File f_even1.3.cpp
  read output from split1.2
  do FFT's and merge
*/

#include "../includes/int_double12.0.h"
#include "../includes/fft_defs.h"
#include "../includes/im_s.h"
#include "f_defs.h"
#include "stdio.h"
#include "assert.h"
#include "qt.h"

//#define PRINT
//#define PRINT1

int NUM_FILES;

void print_usage(const char *cmd)
{
  printf("Usage:- %s (# in files) (in file stub) (out file stub) (even_p) \n",cmd);
  printf("        even_p = 0 => odd characters\n");
  printf("         otherwise => even characters\n"); 
  printf("Exiting.\n");
  exit(1);
}

// perform an in place FFT ,n a power of 2
// F[k]=sum_{j=0..N-1} F[j]e(kj/N) (note not -kj/N)
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
  two_pi_by_B=d_two_pi/B;
  //printf("B=%f\n",B);
  A=d_one/one_over_A;
  pi_A=d_pi*A;
  two_pi_A=d_two_pi*A;
  WS[0]=c_one;
  int_double theta=d_two_pi/(N/NUM_FILES); // make negative to swap FFT/IFFT
  for(int i=1;i<N/(NUM_FILES<<1);i++)
    sin_cos(theta*i,&WS[i].imag,&WS[i].real);
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
  return(fopen(buff,"r+b"));
}

FILE *make_out_file(unsigned int n, const char *stub)
{
  char buff[1024];
  sprintf(buff,"%s%d.dat",stub,n);
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

void save_data2(unsigned int q,unsigned int index,const int_complex &omega,bool neg_one,bool real_p,unsigned int N0,FILE *outfile)
{
  fwrite(&q, sizeof(unsigned int),1,outfile);
  fwrite(&index, sizeof(unsigned int),1,outfile);
  fwrite(&omega, sizeof(int_complex),1,outfile);
  fwrite(&neg_one, sizeof(bool),1,outfile);
  fwrite(&real_p, sizeof(bool),1,outfile);
  fwrite(&N0, sizeof(unsigned int),1,outfile);
}

void save_data3(unsigned int index,const int_complex &omega,FILE *outfile)
{
  fwrite(&index, sizeof(unsigned int),1,outfile);
  fwrite(&omega, sizeof(int_complex),1,outfile);
}


inline bool power_2_p(unsigned int i)
{
  while((i&1)==0)
    i>>=1;
  return(i==1);
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

inline void check_uint(FILE **infiles, unsigned int num_files, unsigned int val)
{
  unsigned int i;
  for(int n=1;n<num_files;n++)
    {
      assert(fread(&i,sizeof(unsigned int),1,infiles[n]));
      assert(i==val);
    }
}
inline void check_luint(FILE **infiles, unsigned int num_files, unsigned long int val)
{
  unsigned long int i;
  for(int n=1;n<num_files;n++)
    {
      assert(fread(&i,sizeof(unsigned long int),1,infiles[n]));
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
 
inline int_complex calc_z(FILE **infiles, unsigned int m, int_double theta, unsigned int no_files)
{
  int_complex e1;
  sin_cos(theta*m,&e1.imag,&e1.real);

  if(no_files==2)
    {
      int_complex z,res;
      assert(fread(&res,sizeof(int_complex),1,infiles[0])); // n=0,2,4...
      assert(fread(&z,sizeof(int_complex),1,infiles[1]));   // n=1,3,5...
      return(res+z*e1);
    }
  if(no_files==4)
    {
      int_complex z,res,e2,e3;
      e2=e1*e1;
      e3=e2*e1;
      assert(fread(&res,sizeof(int_complex),1,infiles[0])); // n=0,4,8...
      assert(fread(&z,sizeof(int_complex),1,infiles[1]));   // n=1,5,9...
      res+=z*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[2]));   // n=2,6,10...
      res+=z*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[3]));   // n=3,7,11...
      return(res+z*e3);
    }
  if(no_files==8)
    {
      int_complex z,res,e2,e4;
      e2=e1*e1;
      e4=e2*e2;
      assert(fread(&res,sizeof(int_complex),1,infiles[0])); // n=0,8...
      assert(fread(&z,sizeof(int_complex),1,infiles[1]));   // n=1,9...
      res+=z*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[2]));   // n=2,10...
      res+=z*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[3]));   // n=3,11...
      res+=z*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[4]));   // n=4,12...
      res+=z*e4;
      assert(fread(&z,sizeof(int_complex),1,infiles[5]));   // n=5,13...
      res+=z*e4*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[6]));   // n=6,14...
      res+=z*e4*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[7]));   // n=7,15...
      return(res+z*e4*e2*e1);
    }
  if(no_files==16)
    {
      int_complex z,res,e2,e4,e8;
      e2=e1*e1;
      e4=e2*e2;
      e8=e4*e4;
      assert(fread(&res,sizeof(int_complex),1,infiles[0])); // n=0,16...
      assert(fread(&z,sizeof(int_complex),1,infiles[1]));   // n=1,17...
      res+=z*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[2]));   // n=2,18...
      res+=z*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[3]));   // n=3,19...
      res+=z*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[4]));   // n=4,20...
      res+=z*e4;
      assert(fread(&z,sizeof(int_complex),1,infiles[5]));   // n=5,21...
      res+=z*e4*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[6]));   // n=6,22...
      res+=z*e4*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[7]));   // n=7,23...
      res+=z*e4*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[8]));   // n=8,24...
      res+=z*e8;
      assert(fread(&z,sizeof(int_complex),1,infiles[9]));   // n=9,25...
      res+=z*e8*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[10]));   // n=10,26...
      res+=z*e8*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[11]));   // n=11,27...
      res+=z*e8*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[12]));   // n=12,28...
      res+=z*e8*e4;
      assert(fread(&z,sizeof(int_complex),1,infiles[13]));   // n=13,29...
      res+=z*e8*e4*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[14]));   // n=14,30...
      res+=z*e8*e4*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[15]));   // n=15,31...
      return(res+z*e8*e4*e2*e1);
    }
  if(no_files==32)
    {
      int_complex z,res,e2,e4,e8,e16;
      e2=e1*e1;
      e4=e2*e2;
      e8=e4*e4;
      e16=e8*e8;
      assert(fread(&res,sizeof(int_complex),1,infiles[0])); // n=0,16...
      assert(fread(&z,sizeof(int_complex),1,infiles[1]));   // n=1,17...
      res+=z*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[2]));   // n=2,18...
      res+=z*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[3]));   // n=3,19...
      res+=z*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[4]));   // n=4,20...
      res+=z*e4;
      assert(fread(&z,sizeof(int_complex),1,infiles[5]));   // n=5,21...
      res+=z*e4*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[6]));   // n=6,22...
      res+=z*e4*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[7]));   // n=7,23...
      res+=z*e4*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[8]));   // n=8,24...
      res+=z*e8;
      assert(fread(&z,sizeof(int_complex),1,infiles[9]));   // n=9,25...
      res+=z*e8*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[10]));   // n=10,26...
      res+=z*e8*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[11]));   // n=11,27...
      res+=z*e8*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[12]));   // n=12,28...
      res+=z*e8*e4;
      assert(fread(&z,sizeof(int_complex),1,infiles[13]));   // n=13,29...
      res+=z*e8*e4*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[14]));   // n=14,30...
      res+=z*e8*e4*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[15]));   // n=15,31...
      res+=z*e8*e4*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[16])); // n=0,16...
      res+=z*e16;
      assert(fread(&z,sizeof(int_complex),1,infiles[17]));   // n=1,17...
      res+=z*e16*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[18]));   // n=2,18...
      res+=z*e16*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[19]));   // n=3,19...
      res+=z*e16*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[20]));   // n=4,20...
      res+=z*e16*e4;
      assert(fread(&z,sizeof(int_complex),1,infiles[21]));   // n=5,21...
      res+=z*e16*e4*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[22]));   // n=6,22...
      res+=z*e16*e4*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[23]));   // n=7,23...
      res+=z*e16*e4*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[24]));   // n=8,24...
      res+=z*e16*e8;
      assert(fread(&z,sizeof(int_complex),1,infiles[25]));   // n=9,25...
      res+=z*e16*e8*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[26]));   // n=10,26...
      res+=z*e16*e8*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[27]));   // n=11,27...
      res+=z*e16*e8*e2*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[28]));   // n=12,28...
      res+=z*e16*e8*e4;
      assert(fread(&z,sizeof(int_complex),1,infiles[29]));   // n=13,29...
      res+=z*e16*e8*e4*e1;
      assert(fread(&z,sizeof(int_complex),1,infiles[30]));   // n=14,30...
      res+=z*e16*e8*e4*e2;
      assert(fread(&z,sizeof(int_complex),1,infiles[31]));   // n=15,31...
      return(res+z*e16*e8*e4*e2*e1);
    }

  printf("Can't handle a %d-way FFT merge yet. Exiting.\n",no_files);
  exit(1);
  return(c_zero);
}

int main(int argc, char **argv)
{
  FILE **infiles,*outfile;
  fpos_t file_pos,*fp;
  unsigned int q,num_prims,N0,index;
  bool even_p,neg_one,real_p;
  int_double err,eta_pi_by_4,theta;
  int_complex omega,z;
  double eta;
  _fpu_rndd();
  //printf("argc=%d\n",argc);
  if(argc!=5)
    print_usage(argv[0]);
  NUM_FILES=atoi(argv[1]);
  printf("Running %s with %d files\n",argv[0],NUM_FILES);
  assert(infiles=(FILE **) malloc(sizeof(FILE *)*NUM_FILES));
  assert(fp=(fpos_t *) malloc(sizeof(fpos_t)*NUM_FILES));
  for(int i=0;i<NUM_FILES;i++)
    assert(infiles[i]=make_file(i,argv[2]));
  even_p=(atoi(argv[4])!=0);
  fread(&q,sizeof(unsigned int),1,infiles[0]);
  printf("q=%d\n",q);
  check_uint(infiles,NUM_FILES,q);
  fread(&num_prims,sizeof(unsigned int),1,infiles[0]);
  printf("# prims=%d\n",num_prims);
  check_uint(infiles,NUM_FILES,num_prims);
  fread(&N,sizeof(unsigned long int),1,infiles[0]);
  printf("N=%lu\n",N);
  check_luint(infiles,NUM_FILES,N);
  assert(Fft_Vec=(int_complex *) _aligned_malloc(sizeof(int_complex)*N/NUM_FILES,16));
  assert(WS=(int_complex *) _aligned_malloc(sizeof(int_complex)*N/2/NUM_FILES,16));
  set_QT_even(q);
  N0=make_N0(q);
  printf("N0=%d\n",N0);
  int_double logq=log(int_double(q));
  setup(q);
  fread(&eta,sizeof(double),1,infiles[0]);
  printf("eta=%10.8e\n",eta);
  check_double(infiles,NUM_FILES,eta);
  eta_pi_by_4=d_pi*eta/4.0;
  for(int f=0;f<NUM_FILES;f++)
    assert(fgetpos(infiles[f],&fp[f])==0);
  //for(int f=0;f<NUM_FILES;f++)
  //printf("file_pos[%d]=%ld\n",f,fp[f]);
  for(int p=0;p<num_prims;p++)
    {
      //printf("reading neg_one for p=%d\n",p);
      fread(&neg_one,sizeof(bool),1,infiles[0]);
      check_bool(infiles,NUM_FILES,neg_one);
      //printf("reading real_p for p=%d\n",p);
      fread(&real_p,sizeof(bool),1,infiles[0]);
      check_bool(infiles,NUM_FILES,real_p);
      if(!real_p)
	p++;
      if(neg_one==even_p)
	continue;
      fread(&index,sizeof(unsigned int),1,infiles[0]);
      //printf("index read =%u\n",index);
      check_uint(infiles,NUM_FILES,index);
      fread(&omega,sizeof(int_complex),1,infiles[0]);
      //print_int_complex_str("omega read =",omega);
      check_int_complex(infiles,NUM_FILES,omega);
      for(int f=0;f<NUM_FILES;f++)
	{
	  assert(fgetpos(infiles[f],&file_pos)==0);
	  fread(Fft_Vec,sizeof(int_complex),N/NUM_FILES,infiles[f]);
	  //printf("Running FFT on index %d with n=%d\n",index,N/NUM_FILES);
	  //print_int_complex_str("before FFT FFT[0]=",Fft_Vec[0]);
	  fft(Fft_Vec,N/NUM_FILES,WS);
	  for(int m=0;m<10;m++) print_int_complex_str("after FFT FFT[]=",Fft_Vec[m]);
	  //printf("not writing fft_vec\n");	  
	  assert(fsetpos(infiles[f],&file_pos)==0);
	  fwrite(Fft_Vec,sizeof(int_complex),N/NUM_FILES,infiles[f]);
	}
      if(real_p)
	continue;
      //printf("Character is complex, so reading its conjugate.\n");
      fread(&index,sizeof(unsigned int),1,infiles[0]);
      //printf("Index read=%u\n",index);
      check_uint(infiles,NUM_FILES,index);
      fread(&omega,sizeof(int_complex),1,infiles[0]);
      //print_int_complex_str("omega read=",omega);
      check_int_complex(infiles,NUM_FILES,omega);
      for(int f=0;f<NUM_FILES;f++)
	{
	  assert(fgetpos(infiles[f],&file_pos)==0);
	  fread(Fft_Vec,sizeof(int_complex),N/NUM_FILES,infiles[f]);
	  //printf("Running FFT on index %d with n=%d\n",index,N/NUM_FILES);
	  //print_int_complex_str("before FFT FFT[0]=",Fft_Vec[0]);
	  fft(Fft_Vec,N/NUM_FILES,WS);
	  //print_int_complex_str("after FFT FFT[0]=",Fft_Vec[0]);
	  assert(fsetpos(infiles[f],&file_pos)==0);
	  fwrite(Fft_Vec,sizeof(int_complex),N/NUM_FILES,infiles[f]);
	}
    }

  // now merge the files together via twiddle factors and
  // multiply by two_Pi_by_B
  // add f_hat_error to each for m=0..N0
  // multiply by exp(x*(d_pi/4.0-eta_pi_by_4))
  //for(int f=0;f<NUM_FILES;f++)
  //printf("file_pos[%d]=%ld\n",f,fp[f]);
  for(int f=0;f<NUM_FILES;f++)
    assert(fsetpos(infiles[f],&fp[f])==0);
  theta=d_two_pi/N;
  for(int p=0;p<num_prims;p++)
    {
      //printf("reading neg_one\n");
      fread(&neg_one,sizeof(bool),1,infiles[0]);
      check_bool(infiles,NUM_FILES,neg_one);
      //printf("reading real_p\n");
      fread(&real_p,sizeof(bool),1,infiles[0]);
      check_bool(infiles,NUM_FILES,real_p);
      if(!real_p)
	p++;
      if(neg_one==even_p)
	continue;
      //printf("reading index\n");

      fread(&index,sizeof(unsigned int),1,infiles[0]);
      check_uint(infiles,NUM_FILES,index);
      //printf("reading omega\n");
      fread(&omega,sizeof(int_complex),1,infiles[0]);
      //print_int_complex_str("Omega=",omega);
      check_int_complex(infiles,NUM_FILES,omega);
      outfile=make_out_file(p,argv[3]);
      save_data2(q,index,omega,neg_one,real_p,N0,outfile);
      for(int f=0;f<NUM_FILES;f++)
	assert(fgetpos(infiles[f],&fp[f])==0);
      
      double x=0.0;
      for(int m=0,n=0;m<N0;m++,x+=one_over_A)
	{
	  z=calc_z(infiles,m,theta,NUM_FILES);
	  z*=two_pi_by_B;
	  err=F_twiddle_e_err(m,eta_pi_by_4,logq);
	  z.real+=err;
	  z.imag+=err;
	  if(!contains_zero(z.imag))
	    {
	      printf("prim # %d, index=%d, m=%d\n",p,index,m);
	      print_int_complex_str("z must be real ",z);
	      printf("Exiting.\n");
	      exit(1);
	    }
	  z.real*=exp(x*(d_pi/4.0-eta_pi_by_4));
	  fwrite(&z.real,sizeof(int_double),1,outfile);
	  n++;
	  if(n==N/NUM_FILES)
	    {
	      for(int f=0;f<NUM_FILES;f++)
		assert(fsetpos(infiles[f],&fp[f])==0);
	      n=0;
	    }
	}
      for(int f=0;f<NUM_FILES;f++)
	{
	  assert(fsetpos(infiles[f],&fp[f])==0);
	  assert(fseek(infiles[f],(N/NUM_FILES)*sizeof(int_complex),SEEK_CUR)==0);
	}
      if(real_p)
	{
	  fclose(outfile);
	  continue;
	}
      fread(&index,sizeof(unsigned int),1,infiles[0]);
      check_uint(infiles,NUM_FILES,index);
      fread(&omega,sizeof(int_complex),1,infiles[0]);
      check_int_complex(infiles,NUM_FILES,omega);
      save_data3(index,omega,outfile);
      for(int f=0;f<NUM_FILES;f++)
	assert(fgetpos(infiles[f],&fp[f])==0);
      x=0.0;
      for(int m=0,n=0;m<N0;m++,x+=one_over_A)
	{
	  z=calc_z(infiles,m,theta,NUM_FILES);
	  z*=two_pi_by_B;
	  err=F_twiddle_e_err(m,eta_pi_by_4,logq);
	  z.real+=err;
	  z.imag+=err;
	  if(!contains_zero(z.imag))
	    {
	      printf("prim # %d, index=%d, m=%d\n",p,index,m);
	      print_int_complex_str("z must be real ",z);
	      printf("Exiting.\n");
	      exit(1);
	    }
	  z.real*=exp(x*(d_pi/4.0-eta_pi_by_4));
	  fwrite(&z.real,sizeof(int_double),1,outfile);
	  n++;
	  if(n==N/NUM_FILES)
	    {
	      for(int f=0;f<NUM_FILES;f++)
		assert(fsetpos(infiles[f],&fp[f])==0);
	      n=0;
	    }
	}
      for(int f=0;f<NUM_FILES;f++)
	{
	  assert(fsetpos(infiles[f],&fp[f])==0);
	  assert(fseek(infiles[f],(N/NUM_FILES)*sizeof(int_complex),SEEK_CUR)==0);
	}
      fclose(outfile);
    }

  for(int f=0;f<NUM_FILES;f++)
    fclose(infiles[f]);

  printf("%s Successful completion.\n",argv[0]);
  


  return(0);
}

