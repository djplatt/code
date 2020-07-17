#include "../includes/int_double12.0.h"
#include "../includes/hurwitz1.0.h"
#include "inttypes.h"

#define two_64 ((double) 18446744073709551616.0)
#define two_32 ((double) 4294967296.0)

#define two_m_101 ((double) 3.9443045261050590270586428264139311483660321755451150238513946533203125e-31)

#define two_m_102 ((double) two_m_101/2.0)

void print_binary1 (int_double x)
{
  uint64_t *i;
  i=(uint64_t *)&x;
  printf(" %lu %lu",i[0],i[1]);
}
void print_binary (int_double x)
{
  uint64_t *i;
  i=(uint64_t *)&x;
  printf(" %lu %lu\n",i[0],i[1]);
}

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
  friend zero_t operator - (const zero_t &lhs, const zero_t &rhs);
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
inline zero_t operator - (const zero_t &lhs, const zero_t &rhs)
{
  zero_t res;
__asm__("mov %7,%2\n\t"
        "sub %8,%2\n\t"
	"mov %5,%1\n\t"
        "sbb %6,%1\n\t"
        "mov %3,%0\n\t"
	"sbb %4,%0\n\t"
	:"=r" (res.a), "=r" (res.b), "=r" (res.c)
	: "m" (lhs.a), "m" (rhs.a), "m" (lhs.b), "m" (rhs.b), "m" (lhs.c), "m" (rhs.c)
	  :);      
 return(res);
}

zero_t uint64_t_to_zero_t (const uint64_t i)
{
  zero_t res;
  res.a=0;
  res.b=0;
  res.c=i<<5;
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

// N_t is number of zeros not including this one
// N(t)=arg gamma_r(s)+S(t)+1
// gamma_r(s)=pi^{-s/2} gamma(s/2)
int_double S(int_double t, uint64_t N_t)
{
  int_complex z=int_complex(0.25,t*0.5);
  int_complex lng=lngamma(z);
  int_double argg=(lng.imag-d_ln_pi*z.imag)/d_pi;
  return(N_t-argg);
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
  int_double sum_rho_re=0.0,sum_rho_im=0.0,sum_gam_2=0.0,sum_gam_3=0.0,sum_rho_abs=0.0,sum_l2_g2=0.0,sum_l_g3=0.0;
  for(uint64_t it=0;it<num_its;it++)
    {
      double st[2];
      uint64_t zs[2];
      // carry_over = start of block - last zero in prev block, 0 if first block
      fread(st,sizeof(double),2,zfile);
      fread(zs,sizeof(uint64_t),1,zfile);
      if(st[0]==0.0)
	{
	  printf("st[0] was 0.0. Exiting.\n");
	  exit(0);
	}
      fread(zs+1,sizeof(uint64_t),1,zfile);

      zero_t sum_zeros=zero_t(0,0,0);
      for(uint64_t z=zs[0]+1;z<=zs[1];z++)
	{
	  zero_t del_t=get_zero(zfile); // exact
	  if((del_t.a==0)&&(del_t.b==0)&&(del_t.c==0))
	    {
	      printf("get_zero returned 0. Exiting.\n");
	      exit(0);
	    }
	  sum_zeros=sum_zeros+del_t; // exact
	  int_double t=st[0]+(int_double) sum_zeros+zero_error;
	  int_double t2=sqr(t);
	  sum_gam_2+=1.0/t2;
	  sum_gam_3+=1.0/(t2*t);
	  sum_l2_g2+=sqr(log(t)/t);
	  sum_l_g3+=log(t)/t/sqr(t);

	  t2+=0.25;
	  sum_rho_re+=0.5/t2;
	  sum_rho_im-=t/t2;
	  sum_rho_abs+=1.0/sqrt(t2);
	}
    }
  printf("%s Sum 1/gamma^2         = ",argv[0]);print_int_double(sum_gam_2);
  print_binary(sum_gam_2);
  printf("%s Sum 1/gamma^3         = ",argv[0]);print_int_double(sum_gam_3);
  print_binary(sum_gam_3);
  printf("%s Sum 1/rho= ",argv[0]);print_int_double(sum_rho_re);printf(" + I* ");print_int_double(sum_rho_im);
  print_binary1(sum_rho_re);
  print_binary(sum_rho_im);
  printf("%s Sum 1/|rho|           = ",argv[0]);print_int_double(sum_rho_abs);
  print_binary(sum_rho_abs);
  printf("%s Sum log^2/gamma^2/    = ",argv[0]);print_int_double(sum_l2_g2);
  print_binary(sum_l2_g2);
  printf("%s Sum log/gamma^3       = ",argv[0]);print_int_double(sum_l_g3);
  print_binary(sum_l_g3);



  fclose(zfile);
  return(0);
}
