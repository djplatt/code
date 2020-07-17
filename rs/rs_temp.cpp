#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <qd/qd_real.h>
#include <qd/fpu.h>

#define A_LIM (200000)
#define C_LIM (100)

#define K (10) // number of RS error terms
#define T (100) // largest power in Taylor Exp of RS
#define STEP (1250) // test every 625'th on average

//#define DEBUG

typedef struct{
  qd_real _log;
  qd_real _sqrt;
} tabs;

tabs table[A_LIM];

qd_real _pi8;

void init_tabs()
{
  for(long unsigned int i=0;i<A_LIM;i++)
    {
      qd_real qi=(double) (A_LIM-i);
      table[i]._log=log(qi);
      table[i]._sqrt=1.0/sqrt(qi);
    }
  _pi8=qd_real::_pi4/2.0;
}

inline qd_real calc_theta(const qd_real &t, const qd_real &t2pi)
{
  qd_real t2=t*t,t5=t2*t2*t,theta=((t2*1680+98.0)*t2+31.0)/80640.0/t5;
  return(theta+t/2.0*(log(t2pi)-1)-_pi8);
}
 
qd_real rs(qd_real t, qd_real coeffs[K+1][T+1])
{
  qd_real t2pi=t/qd_real::_2pi,a=sqrt(t2pi);
  long unsigned int V=floor(a.x[0]);
  //std::cout << "In rs with t=" << t << " a=" << a << " V=" << V << std::endl;
  if(a>=A_LIM)
    {
      printf("rs: t too large. Exiting.\n");
      exit(0);
    }
  qd_real theta=calc_theta(t,t2pi);
  //std::cout << "theta=" << theta << std::endl;
  qd_real res=0.0;
  for(long int n=V;n>0;n--)
    res+=table[A_LIM-n]._sqrt*cos(theta-t*table[A_LIM-n]._log);
  res+=res;

#ifdef DEBUG
  // just do the C0 correction
  qd_real p=a-V;
  qd_real C0=cos(qd_real::_2pi*(p*p-p-1.0/16.0))/cos(qd_real::_2pi*p)/sqrt(a);
  //std::cout << "Main sum=" << res << std::endl;
  //std::cout << "C0 term =" << C0 << std::endl;
  if(V&1)
    return(res+C0);
  else
    return(res-C0);
#endif

  qd_real r[T+1];
  r[0]=1.0;
  r[1]=a-V-0.5; // power series is round r=0.5
  //std::cout << "r=" << r[1] << std::endl;
  for(long unsigned int n=3;n<=T;n++)
    r[n]=1.0;
  for(long unsigned int n=2;n<=T;n<<=1)
    {
      r[n]=r[n>>1]*r[n>>1];
      for(long unsigned int n1=n+1,b=1,m=n1;(n1<n<<1)&&(n1<=T);n1++,m=n1,b=1)
	while(m)
	  {
	    if(m&1)
	      r[n1]*=r[b];
	    m>>=1;b<<=1;
	  }
    }
  //for(long unsigned int n=0;n<=T;n++) std::cout << "p^" << n << "=" << p[n] << std::endl;
      

  //std::cout << "Result of main sum =" << res << std::endl;

  qd_real as[K+1];
  as[0]=1.0/sqrt(a); // (t/2pi)^(-1/4)
  qd_real a_min=1.0/a;
  for(long unsigned int k=1;k<=K;k++)
    as[k]=as[k-1]*a_min;

  //for(long int t1=0;t1<10;t1++) std::cout << "Coefficient of t^" << t1 <<" in C3 is "<< coeffs[3][t1] << std::endl;

  qd_real ksum=0.0;
  for(long int k=K;k>=0;k--)
    {
      qd_real r1=0.0;
      for(long int t1=T;t1>=0;t1--)
	r1+=r[t1]*coeffs[k][t1];
      r1*=as[k];
      ksum+=r1;
      //std::cout << "C" << k << "=" << r1 << " total so far=" << ksum << std::endl; 
    }

  if(V&1)
    res+=ksum;
  else
    res-=ksum;

  return(res);
}

void read_coeffs(qd_real coeffs[K+1][T+1],FILE *cfile)
{
#ifdef DEBUG
  return;
#endif
  long int k,t;
  for(k=0;k<=K;k++)
    for(t=0;t<=T;t++)
      coeffs[k][t]=0.0;

  double d[4];
  fread(&k,sizeof(long int),1,cfile);
  if(k!=K)
    {
      printf("cfile K=%ld, was expecting %ld. Exiting.\n",k,K);
      exit(0);
    }
  fread(&t,sizeof(long int),1,cfile);
  if(t!=T)
    {
      printf("cfile T=%ld, was expecting %ld. Exiting.\n",t,T);
      exit(0);
    }
  for(k=0;k<=K;k++)
    {
      fread(&t,sizeof(long int),1,cfile);
      if(t!=k)
	{
	  printf("Data error in cfile. Exiting.\n");
	  exit(0);
	}
      while(true)
	{
	  fread(&t,sizeof(long int),1,cfile);
	  if(t==-1)
	    break;
	  fread(d,sizeof(double),4,cfile);
	  coeffs[k][t]=qd_real(d[0],d[1],d[2],d[3]);
	}
    }
}

// read our 13 byte structure representing a zero gap
// into a quad
qd_real in_bytes(FILE *infile)
{
  uint64_t a;uint32_t b;uint8_t c;
  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  double r1,r2;
  //printf("Read a=%lX b=%x c=%x\n",a,b,c);
  r1=(double) c/32.0 + (double) b/(32.0*256.0*256.0*256.0*256.0)+(double) (a&0xfff8000000000000)/(32.0*65536.0*65536.0*65536.0*65536.0*65536.0*65536.0); // /2^3, /2^37, /2^101
  r2=(double) (a&0x0007ffffffffffff)/(32.0*65536.0*65536.0*65536.0*65536.0*65536.0*65536.0); // /2^101
  //printf("r1=%50.48e\nr2=%50.48e\n",r1,r2);

  return(qd_real(r1,r2,0,0));
}


int main(int argc, char ** argv)
{
  std::cout << "Command line was:- ";
  for(unsigned long int i=0;i<argc;i++)
    std::cout << argv[i] << " ";
  std::cout << std::endl;

  if(argc!=2)
    {
      printf("Usage:- rs <coeff file>\n");
      exit(0);
    }
  FILE *cfile=fopen(argv[1],"rb");
  if(!cfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }


  unsigned int old_cw;
  fpu_fix_start(&old_cw);
  std::cout.precision(60);

  qd_real coeffs[K+1][T+1];
  read_coeffs(coeffs,cfile);
  
  init_tabs();
  qd_real rl,t=3.0610046000e10;
  rl=rs(t,coeffs);
  std::cout << "RS(" << t <<") returned " << rl << std::endl;
  t=10000.0;
  rl=rs(t,coeffs);
  std::cout << "RS(" << t <<") returned " << rl << std::endl;
  return(0);
}
