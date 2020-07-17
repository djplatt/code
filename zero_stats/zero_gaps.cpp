#include "../includes/int_double12.0.h"
#include "inttypes.h"

#define two_64 ((double) 18446744073709551616.0)
#define two_32 ((double) 4294967296.0)

#define two_m_101 ((double) 3.9443045261050590270586428264139311483660321755451150238513946533203125e-31)

#define two_m_102 ((double) two_m_101/2.0)

class zero_t{
  // a delta is (a+b<<64+c<<96)>>102
public:
  uint64_t a;
  uint32_t b;
  uint64_t c;

  zero_t ()
  {
  }

  zero_t(const uint64_t a1, const uint32_t b1, const uint64_t c1)
  {
    a=a1;b=b1;c=c1;
  }
  /*
  operator double() const
  {
    double ret=this->c*two_32;
    ret+=this->b;
    ret*=two_64;
    ret+=this->a;
    ret*=two_m_101;
    return(ret);
  }
  */
  operator int_double() const
  {
    int_double ret=this->c;
    ret*=two_32;
    ret+=this->b;
    ret*=two_64;
    ret+=this->a;
    ret*=two_m_101;
    return(ret);
  }

  friend zero_t operator + (const zero_t &lhs, const zero_t &rhs);
};

inline zero_t operator + (const zero_t &lhs, const zero_t &rhs)
{
  zero_t res;
__asm__("mov %3,%0\n\t"
        "add %4,%0\n\t"
	"mov %5,%1\n\t"
        "adc %6,%1\n\t"
        "mov %7,%2\n\t"
	"adc %8,%2\n\t"
	:"=r" (res.a), "=r" (res.b), "=r" (res.c)
	: "m" (lhs.a), "m" (rhs.a), "m" (lhs.b), "m" (rhs.b), "m" (lhs.c), "m" (rhs.c)
	  :);      
 return(res);
}

zero_t get_zero(FILE *zfile)
{
  zero_t res;
  uint8_t c;
  if(fread(&res.a,sizeof(uint64_t),1,zfile)!=1)
    {
      printf("Error reading 64 bit part of zero from file. Exiting.\n");
      exit(0);
    }
  if(fread(&res.b,sizeof(uint32_t),1,zfile)!=1)
    {
      printf("Error reading 32 bit part of zero from file. Exiting.\n");
      exit(0);
    }
  if(fread(&c,sizeof(uint8_t),1,zfile)!=1)
    {
      printf("Error reading 8 bit part of zero from file. Exiting.\n");
      exit(0);
    }
  res.c=c;
  return(res);
}

int main(int argc, char **argv)
{

  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>.\n",argv[0]);
      exit(0);
    }

  _fpu_rndd();

  int_double zero_error=int_double(-two_m_102,two_m_102);
  int_double diff_error=int_double(-two_m_101,two_m_101);

  FILE* zfile=fopen(argv[1],"rb");
  if(!zfile)
    {
      printf("Error opening file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  uint64_t num_its;
  if(fread(&num_its,sizeof(uint64_t),1,zfile)!=1)
    {
      printf("Error reading number of iterations from zeros file. Exiting.\n");
      exit(0);
    }
  double st[2],t0;
  int_double t,last_t;
  uint64_t zs[2],z,it;
  zero_t sum_zeros,del_t;
  for(it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,zfile);
      fread(zs,sizeof(uint64_t),1,zfile);
      if(st[0]==0.0)
	{
	  printf("st[0] was 0.0. Exiting.\n");
	  exit(0);
	}
      fread(zs+1,sizeof(uint64_t),1,zfile);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  if((it==0)&&(z==zs[0]+1))
	    {
	      sum_zeros=get_zero(zfile);
	      if((sum_zeros.a==0)&&(sum_zeros.b==0)&&(sum_zeros.c==0))
		{
		  printf("get_zero returned 0. Exiting.\n");
		  exit(0);
		}
	      last_t=st[0]+int_double(sum_zeros)+zero_error;
	      continue;
	    }
	  del_t=get_zero(zfile); // exact
	  if((del_t.a==0)&&(del_t.b==0)&&(del_t.c==0))
	    {
	      printf("get_zero returned 0. Exiting.\n");
	      exit(0);
	    }
	  sum_zeros=sum_zeros+del_t; // exact
	  t=st[0]+(int_double) sum_zeros+zero_error;
	  int_double mid_pt=(t+last_t)/2.0;
	  int_double norm_gap=(int_double)del_t*log(mid_pt)/d_two_pi;
	  printf("%10.8ee\n",(norm_gap.left-norm_gap.right)/2.0);
	  last_t=t;
	}
    }
  fclose(zfile);
  return(0);
}
