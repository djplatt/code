#include "../includes/int_double11.0.h"
#include "../includes/fft_defs.h"
#include "../includes/im_s.h"
#include "stdio.h"
#include "assert.h"

//#define PRINT
//#define PRINT1

#define N (1<<20)
#define one_over_A ((double) 5.0/64.0)
#define B ((double) ((double) N)*one_over_A)

#define POS (0)
#define NEG (1)
#define UNK (2)
typedef char sign_t;

sign_t sign(int_double x)
{
  if(x.left>0.0)
    return(POS);
  if(x.right>0.0)
    return(NEG);
  return(UNK);
}

void print_sign(int_double x)
{
  if(x.left>0.0)
    {
      printf("+");
    }
  else
    {
      if(x.right>0.0)
	printf("-");
      else
	printf("?");
    }
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


int_double A;
int_double two_pi_A;
int_double four_pi_A;
int_complex *WS;
int_complex chis[MAX_Q];
int_complex *Fft_Vec;
int_double d_log_pi;
int_double d_log2;
int_double d_log_zeta_9_by_8;
int_double d_small;
int_complex c_small;
int_double exp_minus_pi_A;
int_double d_pi_by_q;

void setup(unsigned int q)
{
  assert(one_over_A<2*M_PI);
  A=d_one/one_over_A;
  two_pi_A=d_two_pi*A;
  four_pi_A=two_pi_A*2;
  WS[0]=c_one;
  int_double omega=d_two_pi/N;
  for(int i=1;i<N/2;i++)
    {
      sin_cos(omega*i,&WS[i].imag,&WS[i].real);
      //if((i%1000)==0)
      //printf("rel error WS[%d]=%d\n",i,-rel_error(WS[i].real));
      //print_int_complex_str("WS[]=",WS[i]);
    }
  chis[0]=c_zero;
  chis[1]=c_one;
  chis[2]=-c_one;
  d_log_pi=log(d_pi);
  d_log2=log(int_double(2));
  // this is log(Zeta(9/8)*(1/2)^(5/16))
  d_log_zeta_9_by_8=log(int_double(1575824,1575825)/1000000);
  double small=nextafter(0.0);
  d_small=int_double(-small,small);
  c_small=int_complex(d_small,d_small);
  exp_minus_pi_A=d_one-exp(-d_pi*A);
  d_pi_by_q=d_pi/q;
}

int_complex L(unsigned int q, const int_double &t, const double eta, int_complex *chis, int_complex epsilon)
{
  int_complex s=int_complex(int_double(0.5),t);
  int_complex res=c_zero;
  for(int i=1;i<q;i++)
    if(co_prime(i,q))
      res+=chis[i]*hurwitz(s,int_double(i)/q);
  res*=pow(int_double(q),int_complex(int_double(-0.5),-t/2));
  res*=pow(d_pi,(-s-1)/2);
  res*=exp(lngamma((s+1)/2)+d_pi*t*eta/4);
  res*=epsilon;
  return(res);
}

inline bool power_2_p(unsigned int i)
{
  while((i&1)==0)
    i>>=1;
  return(i==1);
}

int_complex make_epsilon(unsigned int q, unsigned int index, unsigned int pr, unsigned int phi, int_complex *chis, bool even_p)
{
  int_double omega=d_two_pi/q;
  int_complex w;
  int_complex res=c_zero;
  for(int n=1;n<q;n++)
    if(co_prime(n,q))
      {
	sin_cos(omega*n,&w.imag,&w.real);
	res+=chis[n]*w;
      }
  res/=sqrt(int_double(q));
  if(even_p)
    return(conj(sqrt(res)));
  // multiply by -i
  w.real=res.imag;
  w.imag=-res.real;
  return(conj(sqrt(w)));
}


// q is 2^a with a>=3
void make_chis1(unsigned int q, unsigned int index, unsigned int pr, unsigned int phi, int_complex *chis)
{
  unsigned int pow,chi_ptr,n,phi_by_2=phi>>1;
  int even_p=(index<phi_by_2);
  int_double four_pi_over_phi=d_two_pi/(phi/2),theta;
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

inline int_double X_x(const int_double &x, const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  return(exp(x*2+logdelta+d_log_pi-logq-delta));
}

inline int_double log_X_x(const int_double &x, const int_double &delta, const int_double &logdelta, const int_double &logq)
{
  return(x*2+logdelta+d_log_pi-logq-delta);
}

int_double E_o (const int_double &t, const int_double &eta_pi_by_4, const int_double &logq)
{
  int_complex lng=lngamma(int_complex(int_double(0.75),t/2));
  int_double res=exp(lng.real+eta_pi_by_4*t+logq*0.3125-d_log_pi*0.5625+d_log_zeta_9_by_8+log(t*t+9.0/4.0)*0.15625);
  //print_int_double_str("E_o called with t=",t);
  //print_int_double_str("Returning",res);
  return(res);
}

int_double beta_o(const int_double &t)
{
  int_double res=d_pi/4-atan2(d_one,abs(t)*2)*1.5-d_one/d_pi/d_pi/abs(t*2-9.0/4.0)*4;
  return(res);
}

int_complex F_twiddle_o_err(const unsigned int m, const int_double &eta_pi_by_4, const int_double &logq, bool *res_ok)
{
  int_double t1=int_double(m)/A+B;
  int_double t2=int_double(m)/A-B;
  int_double beta1=beta_o(t1);
  if(beta1.left<=-eta_pi_by_4.right)
    {
      printf("Beta1 fails test on m=%d\n",m);
      print_int_double_str("Beta()=",beta1);
      res_ok[0]=false;
      return(c_zero);
    } 
  int_double beta2=beta_o(t2);
  if(beta2.left<=eta_pi_by_4.right)
    {
      printf("Beta2 fails test on m=%d\n",m);
      res_ok[0]=false;
      return(c_zero);
    }
  int_double res=E_o(t1,eta_pi_by_4,logq)/(d_one-exp(-(beta1-eta_pi_by_4)*B));
  res+=E_o(t2,eta_pi_by_4,logq)/(d_one-exp(-(beta2+eta_pi_by_4)*B));
  res.left=res.right;
  return(int_complex(res,res));
  //return(c_zero);
}


int_complex F_hat_twiddle_o_err(const int_double &x, const int_double &delta, const int_double &logdelta, const int_double &logq)
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
  //print_int_double_str("F_hat_twiddle_o_err returning",res);
  return(int_complex(res,res));
  //return(c_zero);
}

#ifdef PRINT1
int F_hat_counter=0;
#endif

// error from truncating sum of F_hat to M terms
inline int_complex F_hat_o_err(const int_double &x,const unsigned int q, const int_double &logq, unsigned int *M, const int_double &delta, const int_double &logdelta)
{
  int_double logXx=log_X_x(x,delta,logdelta,logq);
  int_double Xx=exp(logXx);
  int_double M_min=exp(logXx*(-0.5))*6.0; // sqrt(36/Xx)
  M[0]=ceil(-M_min.right);
  int_double XxMM=Xx*M[0]*M[0];
  int_double res=exp(d_log2+x*1.5-XxMM-logq*0.75-logdelta*0.5-logXx+log(d_one+d_one/(XxMM*2.0))*0.5);
  res.left=res.right;
#ifdef PRINT1
  if((F_hat_counter%1000)==0)
    {
      print_int_double_str("Log(X(x))=",logXx);
      print_int_double_str("X(x)=",Xx);
      print_int_double_str("X(x)*M^2=",XxMM);
      printf("M=%d\n",M[0]);
      print_int_double_str("x=",x);
      print_int_double_str("F_hat_error=",res);
    }
  F_hat_counter++;
#endif
  return(int_complex(res,res));
  //return(c_zero);
}


inline int_complex F_hat_o(const int_double &x, const unsigned int q, const int_double &logq, const int_complex &epsilon, int_complex *chis,
		    const int_double &eta_pi_by_4, const int_double &delta, const int_double &logdelta)
{
  unsigned int M;
  int_complex res=F_hat_o_err(x,q,logq,&M,delta,logdelta);
  //print_int_complex_str("F_hat_o error=",res);
  //printf("M=%d\n",M);
  int_complex res1;
  int_complex u_x=int_complex(int_double(x),eta_pi_by_4);
  int_complex exp_2_u_x=exp(u_x*2)*d_pi_by_q; // Pi/q*exp(2*U(x))
  int_complex inner_exp;

  int_complex outer_exp=u_x*1.5+d_log2-logq*0.75;

  // make this sum the other way?
  for(int n=1,r=1,n_sqr;n<=M;n++,r++)
    {
      if(r==q)
	{
	  r=0;
	  continue;
	}
      if(!co_prime(r,q))
	continue;
      n_sqr=n*n;
      inner_exp=outer_exp-exp_2_u_x*n_sqr+log(int_double(n));
                
      if(inner_exp.real.right>745) // exp(-745)<d_small
	{
	  //printf("saving %d\n",M-n+1);
	  res+=c_small*(M-n+1);
	  break;
	}
	res1=exp(inner_exp)*chis[r];
	res+=res1;
#ifdef PRINT1
	if(M==4)
	  {
	    print_int_complex_str("inner_exp=",inner_exp);
	    print_int_complex_str("chi(r)=",chis[r]);
	    print_int_complex_str("This term=",res1);
	    print_int_complex_str("New total=",res);
	  }
#endif
    }
  //printf("after summing:\nrel_error(res)=%d\n",-rel_error(res.real));
  //print_int_complex_str("",res)
  //exit(0);
  res*=epsilon;
#ifdef PRINT1
	if(M==4)
	  {
	    print_int_complex_str("Final res=",res);
	    print_int_complex_str("outer_exp=",outer_exp);
	    exit(0);
	  }
#endif

  //printf("after epsilon:\nrel_error(res)=%d\n",-rel_error(res.real));
  //print_int_double_str("x=",x);
  //print_int_complex_str("F_hat(x)=",res);
  //res+=err;
  //printf("Fo error: rel_error(res)=%d\n",-rel_error(res.real));
  //print_int_complex_str("F_hat_o returning",res);
  return(res);
}

int main()
{
  _fpu_rndd();

  FILE *ofile;


  unsigned int q=9973;
  unsigned int phi=9972;
  unsigned int pr=11;
  unsigned int N0;
  unsigned int two=2;
  im_s this_im;
  assert(Fft_Vec=(int_complex *) _aligned_malloc(sizeof(int_complex)*N,16));
  assert(WS=(int_complex *) _aligned_malloc(sizeof(int_complex)*N/2,16));


  ofile=fopen("foo.dat","wb");
  setup(q);
  int_double logq=log(int_double(q));
  double eta=16770505/(double) (1<<24);
  double t0=ceil(10000000.0/q)*10.0;
  N0=(unsigned int) t0/one_over_A;
  fwrite(&N0,sizeof(unsigned int),1,ofile);
  for(int i=0;i<N0;i++)
    {
      this_im.im_s=i*one_over_A;
      fwrite(&this_im,sizeof(im_s),1,ofile);
    }
  fwrite(&q,sizeof(unsigned int),1,ofile);
  fwrite(&two,sizeof(unsigned int),1,ofile);
  printf("Looking for zeros to height %8.7e\n",t0);
  printf("Trying to use %f of the vector.\n",(double) N0/ (double) N);
  int_double delta=d_pi_2*(1.0-fabs(eta));
  int_double logdelta=log(delta);
  int_double eta_pi_by_4=d_pi*eta/4;
  int_double two_pi_by_B=d_two_pi/B;
  int_double n_two_pi_by_B;
  bool res_ok=true;
  bool even_p,odd_p;
  int index=1;
  while(true)
    {
      fwrite(&index,sizeof(unsigned int),1,ofile);
      printf("Running on q=%d index = %d.\n",q,index);
      make_chis(q,index,pr,phi,1,chis);
      //for(int i=1;i<phi;i++)
      //  printf("rel error of Re chi[%d]=%d.\n",i,rel_error(chis[i].real));
      //exit(0);
      even_p=(contains_zero(chis[q-1]-c_one));
      odd_p=!even_p;
      int_complex epsilon=make_epsilon(q,index,pr,phi,chis,even_p);
      fwrite(&epsilon,sizeof(int_complex),1,ofile);
      fwrite(&odd_p,sizeof(bool),1,ofile);
      //printf("rel error Re epsilon=%d\n",-rel_error(epsilon.real));
      //exit(0);
      //print_int_complex_str("epsilon=",epsilon);
      if(even_p)
	exit(0);
      
      printf("Running with:-\n\tN\t\t=\t%d\n\tq\t\t=\t%d\n\teta\t\t=\t%f\n\t1/A\t\t=\t%f\n\tB\t\t=\t%f\n",N,q,eta,one_over_A,B);
      printf("\tdelta\t\t=\t");print_int_double(delta);
      printf("\n\tlog(delta)\t=\t");print_int_double(logdelta);
      printf("\n");
      
      Fft_Vec[0]=F_hat_o(d_zero,q,logq,epsilon,chis,eta_pi_by_4,delta,logdelta)+
	F_hat_twiddle_o_err(two_pi_A,delta,logdelta,logq)*2;

      //printf("rel_error(F_o(%f))=%d\t",(double) 0 * one_over_A,-rel_error(Fft_Vec[0].real));
      //print_int_complex_str("",Fft_Vec[0]);
      
      for(int n=1;n<N/2;++n)
	{
	  n_two_pi_by_B=two_pi_by_B*n;
	  Fft_Vec[n]=F_hat_o(n_two_pi_by_B,q,logq,epsilon,chis,eta_pi_by_4,delta,logdelta);
	  Fft_Vec[n]+=F_hat_twiddle_o_err(n_two_pi_by_B+two_pi_A, delta, logdelta, logq);
	  Fft_Vec[n]+=F_hat_twiddle_o_err(two_pi_A-n_two_pi_by_B, delta, logdelta, logq);
	  //if(fabs(Fft_Vec[n].real.left)<1e-300)
	  //  Fft_Vec[n]=c_zero;
	  Fft_Vec[N-n]=conj(Fft_Vec[n]);
#ifdef PRINT
	  if((n%1000)==0)
	    {
	      printf("rel_error(F_o(%f))=%d\t",(double) n * one_over_A,-rel_error(Fft_Vec[n].real));
	      print_int_complex_str("",Fft_Vec[n]);
	    }
#endif
	}
      n_two_pi_by_B=two_pi_by_B*(N/2);
      Fft_Vec[N/2]=F_hat_o(n_two_pi_by_B,q,logq,epsilon,chis,eta_pi_by_4,delta,logdelta)+
	F_hat_twiddle_o_err(n_two_pi_by_B+two_pi_A, delta, logdelta, logq)+
	F_hat_twiddle_o_err(two_pi_A-n_two_pi_by_B,delta,logdelta,logq);

      printf("Running fft.\n");
      
      fft(Fft_Vec,N,WS);

      for(int i=0;i<N0;i++)
	Fft_Vec[i]*=two_pi_by_B;

      printf("fft finished.\n");

      //Fft_Vec[0]+=F_twiddle_o_err(0,eta_pi_by_4,logq,&res_ok);
      //print_int_complex_str("FFT[0]=",Fft_Vec[0]);
      //print_int_complex_str("F(0)  =",L(q,d_zero,eta,chis,epsilon));
      res_ok=true;
      for(int i=0;i<N0;i++)
	{
	  Fft_Vec[i]+=F_twiddle_o_err(i,eta_pi_by_4,logq,&res_ok);
	  fwrite(&Fft_Vec[i].real,sizeof(int_double),1,ofile);
	  if(res_ok)
	    {
	      if(!contains_zero(Fft_Vec[i].imag))
		{
		  printf("Fft_Vec[%d] not real.\n",i);
		  print_int_complex_str("=",Fft_Vec[i]);
		  print_int_complex_str("F(t)=",L(q,int_double(i)*one_over_A,eta,chis,epsilon));
		  break;
		}
#ifdef PRINT
	      if((i%1000)==0)
		{
		  printf("rel_error(F_o(%f))=%d\t",(double) i * one_over_A,-rel_error(Fft_Vec[i].real));
		  print_int_complex_str("",Fft_Vec[i]);
		}
#endif
	    }
	  else
	    {
	      res_ok=true;
	      Fft_Vec[i]=c_zero;
	    }
	}
      if(index==9971)
	break;
      index=9971;
    }

    

  /*  
  sign_t last_sign,this_sign;
  unsigned int zero_count;
  unsigned int k=0;
  zero_count=0;
  while(true)
    {
      this_sign=sign(Fft_Vec[k++].real);
      if(this_sign!=UNK)
	{
	  last_sign=this_sign;
	  break;
	}
    }
  printf("N0=%d\n",N0);
  while(k<=N0)
    {
      this_sign=sign(Fft_Vec[k++].real);
      if(this_sign==UNK)
	continue;
      if(this_sign==last_sign)
	continue;
      //printf("Zero Found at %d\n",k);
      last_sign=this_sign;
      zero_count++;
    }
  printf("I found %d zeros.\n",zero_count);

  //      exit(0);
  */
  fclose(ofile);
  return(0);
}
