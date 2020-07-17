#include "../includes/int_double12.0.h"
#include "../includes/hurwitz1.0.h"
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
  double st[2],t0;
  uint64_t zs[2],z,it;
  int calced=0;
  zero_t sum_zeros;
  int_double t,last_t;
  int_double sum_rho=0.0; // sum 1/rho+1/(1-rho)
  int_double sum_gam=0.0; // sum 2/gamma
  double min_abs_smallest_gap=DBL_MAX,max_abs_smallest_gap=DBL_MAX;
  double min_abs_largest_gap=0.0,max_abs_largest_gap=0.0;
  double min_rel_smallest_gap=DBL_MAX,max_rel_smallest_gap=DBL_MAX;
  double min_rel_largest_gap=0.0,max_rel_largest_gap=0.0;
  double min_smallest_st=DBL_MAX,max_smallest_st=DBL_MAX;
  double min_largest_st=0.0,max_largest_st=0.0;
  uint64_t z_rec[12];
  for(it=0;it<num_its;it++)
    {
      sum_zeros=zero_t(0,0,0);
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
	  zero_t del_t=get_zero(zfile); // exact
	  if((del_t.a==0)&&(del_t.b==0)&&(del_t.c==0))
	    {
	      printf("get_zero returned 0. Exiting.\n");
	      exit(0);
	    }
	  sum_zeros=sum_zeros+del_t; // exact
	  t=st[0]+(int_double) sum_zeros+zero_error;
	  sum_gam+=1.0/t;
	  sum_rho+=1.0/(0.25+sqr(t));
	  int_double S_t=S(t,z-1);

	  //print_int_double_str("S(t)=",S_t);
	  if(S_t.left>min_largest_st)
	    {z_rec[0]=z;min_largest_st=S_t.left;}
	  if(-S_t.right>max_largest_st)
	    {z_rec[1]=z;max_largest_st=-S_t.right;}
	  S_t-=1;
	  //print_int_double_str("S(t)=",S_t);
	  if(S_t.left<min_smallest_st)
	    {z_rec[2]=z;min_smallest_st=S_t.left;}
	  if(-S_t.right<max_smallest_st)
	    {z_rec[3]=z;max_smallest_st=-S_t.right;}
	  if(!calced)
	    {
	      printf("0: First Zero %lu is at ",z);print_int_double_str("",t);
	      calced=1;
	    }
	  else
	    {
	      int_double diff=(int_double) del_t+diff_error;
	      if(diff.left>min_abs_largest_gap)
		{z_rec[4]=z;min_abs_largest_gap=diff.left;}
	      if(-diff.right>max_abs_largest_gap)
		{z_rec[5]=z;max_abs_largest_gap=-diff.right;}
	      if(diff.left<min_abs_smallest_gap)
		{z_rec[6]=z;min_abs_smallest_gap=diff.left;}
	      if(-diff.right<max_abs_smallest_gap)
		{z_rec[7]=z;max_abs_smallest_gap=-diff.right;}
	      int_double rel_diff=diff*log(last_t)/d_two_pi;
	      if(rel_diff.left>min_rel_largest_gap)
		{z_rec[8]=z;min_rel_largest_gap=rel_diff.left;}
	      if(-rel_diff.right>max_rel_largest_gap)
		{z_rec[9]=z;max_rel_largest_gap=-rel_diff.right;}
	      if(rel_diff.left<min_rel_smallest_gap)
		{z_rec[10]=z;min_rel_smallest_gap=rel_diff.left;}
	      if(-rel_diff.right<max_rel_smallest_gap)
		{z_rec[11]=z;max_rel_smallest_gap=-rel_diff.right;}
	    }	      
	  last_t=t;
	}
    }
  printf("1: Last Zero %lu is at ",z-1);print_int_double_str("",t);

  printf("2: There is no absolute gap larger than            %20.18e %lu\n",max_abs_largest_gap,z_rec[5]);
  printf("3: There is an absolute as large as                %20.18e %lu\n",min_abs_largest_gap,z_rec[4]);
  printf("4: There is no absolute gap smaller than           %20.18e %lu\n",min_abs_smallest_gap,z_rec[6]);
  printf("5: There is an absolute gap as small as            %20.18e %lu\n",max_abs_smallest_gap,z_rec[7]);


  printf("6: There is no relative gap larger than            %20.18e %lu\n",max_rel_largest_gap,z_rec[9]);
  printf("7: There is a relative as large as                 %20.18e %lu\n",min_rel_largest_gap,z_rec[8]);
  printf("8: There is no relative gap smaller than           %20.18e %lu\n",min_rel_smallest_gap,z_rec[10]);
  printf("9: There is a relative gap as small as             %20.18e %lu\n",max_rel_smallest_gap,z_rec[11]);

  printf("10: There is an S(t) larger than                   %20.18e %lu\n",min_largest_st,z_rec[0]);
  printf("11: There is no S(t) larger than                   %20.18e %lu\n",max_largest_st,z_rec[1]);
  printf("12: There is no S(t) smaller than                  %20.18e %lu\n",min_smallest_st,z_rec[2]);
  printf("13: There is an S(t) smaller than                  %20.18e %lu\n",max_smallest_st,z_rec[3]);
 
  print_int_double_str("Sum 1/gamma         = ",sum_gam);
  print_int_double_str("Sum 1/rho+1/(1-rho) = ",sum_rho);

  fclose(zfile);
  return(0);
}
