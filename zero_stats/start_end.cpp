#include "../includes/int_double12.0.h"
#include "../includes/hurwitz1.0.h"
#include "inttypes.h"

#define two_64 ((double) 18446744073709551616.0)
#define two_32 ((double) 4294967296.0)

#define two_m_101 ((double) 3.9443045261050590270586428264139311483660321755451150238513946533203125e-31)

#define two_m_102 ((double) two_m_101/2.0)

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

void print_zero(const zero_t &z)
{
  printf("%lu %u %lu\n",z.a,z.b,z.c);
}

int main(int argc, char **argv)
{

  if(argc!=2)
    {
      printf("Usage:- %s <sorted start end file>.\n",argv[0]);
      exit(0);
    }

  _fpu_rndd();

  int_double zero_error=int_double(-two_m_102,two_m_102);
  int_double diff_error=int_double(-two_m_101,two_m_101);

  FILE* zfile=fopen(argv[1],"r");
  if(!zfile)
    {
      printf("Error opening file %s for character input. Exiting.\n",argv[1]);
      exit(0);
    }

  double min_abs_smallest_gap=DBL_MAX,max_abs_smallest_gap=DBL_MAX;
  double min_abs_largest_gap=0.0,max_abs_largest_gap=0.0;
  double min_rel_smallest_gap=DBL_MAX,max_rel_smallest_gap=DBL_MAX;
  double min_rel_largest_gap=0.0,max_rel_largest_gap=0.0;
  double min_smallest_st=DBL_MAX,max_smallest_st=DBL_MAX;
  double min_largest_st=0.0,max_largest_st=0.0;
  uint64_t z_rec[12];


  uint64_t n1,n2,line=1;
  double d1,d2;
  zero_t z1,z2;
  if(fscanf(zfile,"%lu %lf %lu %u %lu\n",&n1,&d1,&z1.a,&z1.b,&z1.c)!=5)
    {
      printf("Error reading first line of %s. Exiting.\n",argv[0]);
      exit(0);
    }
  while(true)
    {
      line++;
      if(fscanf(zfile,"%lu %lf %lu %u %lu\n",&n1,&d1,&z1.a,&z1.b,&z1.c)!=5)
	{
	  printf("Error reading %lu'th line of %s. Exiting.\n",line,argv[0]);
	  exit(0);
	}
      if(fscanf(zfile,"%lu %lf %lu %u %lu\n",&n2,&d2,&z2.a,&z2.b,&z2.c)!=5)
	{
	  printf("Last zero read was %lu on line %lu. Finished.\n",n1,line);
	  break;
	}
      line++;
      z1.c+=d1*32; // normalise the zeros by including the double prec offset
      z2.c+=d2*32; // into the top word
      zero_t zdiff=z2-z1;
      //printf("z1=");print_zero(z1);
      //printf("z2=");print_zero(z2);
      //printf("z2-z1=");print_zero(zdiff);

      int_double diff=zdiff;
      if(diff.left>min_abs_largest_gap)
	{z_rec[4]=n1;min_abs_largest_gap=diff.left;}
      if(-diff.right>max_abs_largest_gap)
	{z_rec[5]=n1;max_abs_largest_gap=-diff.right;}
      if(diff.left<min_abs_smallest_gap)
	{z_rec[6]=n1;min_abs_smallest_gap=diff.left;}
      if(-diff.right<max_abs_smallest_gap)
	{z_rec[7]=n1;max_abs_smallest_gap=-diff.right;}
      int_double last_t=z1;
      int_double rel_diff=diff*log(last_t)/d_two_pi;
      if(rel_diff.left>min_rel_largest_gap)
	{z_rec[8]=n1;min_rel_largest_gap=rel_diff.left;}
      if(-rel_diff.right>max_rel_largest_gap)
	{z_rec[9]=n1;max_rel_largest_gap=-rel_diff.right;}
      if(rel_diff.left<min_rel_smallest_gap)
	{z_rec[10]=n1;min_rel_smallest_gap=rel_diff.left;}
      if(-rel_diff.right<max_rel_smallest_gap)
	{z_rec[11]=n1;max_rel_smallest_gap=-rel_diff.right;}

    }
  /*
  for(it=0;it<num_its;it++)
    {
      // carry_over = start of block - last zero in prev block, 0 if first block
      fread(st,sizeof(double),2,zfile);
      fread(zs,sizeof(uint64_t),1,zfile);
      if(st[0]==0.0)
	{
	  printf("st[0] was 0.0. Exiting.\n");
	  exit(0);
	}
      fread(zs+1,sizeof(uint64_t),1,zfile);

      sum_zeros=zero_t(0,0,0);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  zero_t del_t=get_zero(zfile); // exact
	  if((del_t.a==0)&&(del_t.b==0)&&(del_t.c==0))
	    {
	      printf("get_zero returned 0. Exiting.\n");
	      exit(0);
	    }
	  sum_zeros=sum_zeros+del_t; // exact
	  if(z==zs[0]+1) // first zero in block, so add delta from last block
	    del_t=del_t+carry_over;
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
	      printf("0: First Zero %lu is at %10.1f + %lu %lu %lu = ",z,st[0],sum_zeros.a,sum_zeros.b,sum_zeros.c);print_int_double_str("",t);
	      calced=1;
	    }
	  else
	    {
	      int_double diff=(int_double) del_t+diff_error;
	    }	      
	  last_t=t;
	}
      double width=st[1]-st[0];
      uint64_t iwidth=width;
      if(width-iwidth!=0.0)
	{
	  printf("File contained a region of non-integral width. Exiting.\n");
	  exit(0);
	}
      zero_t zwidth=uint64_t_to_zero_t(iwidth);
      //printf("width=%10.1f iwidth=%lu zwidth= %lu %lu %lu sum_zeros= %lu %lu %lu\n",width,iwidth,zwidth.a,zwidth.b,zwidth.c,sum_zeros.a,sum_zeros.b,sum_zeros.c);
      carry_over=zwidth-sum_zeros;
      //printf("Carry over= %lu %lu %lu\n",carry_over.a,carry_over.b,carry_over.c);    
    }
  */
  //printf("1: Last Zero %lu is at %10.1f + %lu %lu %lu = ",z-1,st[0],sum_zeros.a,sum_zeros.b,sum_zeros.c);print_int_double_str("",t);

  printf("2: There is no absolute gap larger than            %20.18e %lu\n",max_abs_largest_gap,z_rec[5]);
  printf("3: There is an absolute as large as                %20.18e %lu\n",min_abs_largest_gap,z_rec[4]);
  printf("4: There is no absolute gap smaller than           %20.18e %lu\n",min_abs_smallest_gap,z_rec[6]);
  printf("5: There is an absolute gap as small as            %20.18e %lu\n",max_abs_smallest_gap,z_rec[7]);


  printf("6: There is no relative gap larger than            %20.18e %lu\n",max_rel_largest_gap,z_rec[9]);
  printf("7: There is a relative as large as                 %20.18e %lu\n",min_rel_largest_gap,z_rec[8]);
  printf("8: There is no relative gap smaller than           %20.18e %lu\n",min_rel_smallest_gap,z_rec[10]);
  printf("9: There is a relative gap as small as             %20.18e %lu\n",max_rel_smallest_gap,z_rec[11]); 
  fclose(zfile);
  return(0);
}
