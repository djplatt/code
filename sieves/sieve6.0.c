#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"
#include "../includes/pi_x.h"
#include "gmp.h"

// Version 6.0
// Uses Fermat pseudoprime tests

#define debug printf("Reached line number %d.\n",__LINE__)



inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}

void print_usage()
{
  printf("usage:- sieve6.0 <sieve_num>  <num_its> <outfile>\n");
  printf("<num_its> is number of 2^34 sieves, <sieve_num>=0 starts at x-(<num_sieves>-<sieve_num>)/*2^33\n");
  exit(0);

}

inline long unsigned int gcd (long unsigned int a, long unsigned int b)
/* Euclid algorithm gcd */
{
	long unsigned int c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

inline int co_prime(long unsigned int a, long unsigned int b)
{
	return(gcd(a,b)==1);
};

inline __uint64_t rem_80_48(__uint128_t x, __uint64_t y)
{
  
  __uint64_t w=x&0xFFFF,z=x>>16;
  return((((z%y)<<16)+w)%y);
  
  /*
  // inline assembler is actually slower!
  __uint64_t temp,*_x,x1,x2;
  _x=(__uint64_t *) &x;
  //printf("rem_80_40 called with x=");print_bigint(x);printf(" and y=%lu\n",y);
  x1=_x[0];
  x2=_x[1];
  //printf("x1=%lu x2=%lu y=%lu\n",x1,x2,y);
  __asm__("movq %2,%%rax\n\t"
	  "xor %%rdx,%%rdx\n\t"
	  "divq %3\n\t"
	  "movq %1,%%rax\n\t"
	  "divq %3\n\t"
	  "movq %%rdx,%0\n\t"
	  :"=r" (temp)
	  :"r" (x1), "r" (x2), "r" (y)
	  :"rax", "rdx");
  //printf("rem_80_40 returning %lu\n",temp);
  return(temp);
  */
}


mpz_t mpz_2,res,mp,mp1,mp_1;



int pseudoprimep(bigint p)
{
  unsigned long int *ps;
  //printf("In pseudoprime with p=");print_bigint(p);printf("\n");
  ps=(unsigned long int *) &p;
  mpz_set_ui(mp,ps[1]);
  mpz_mul_2exp(mp1,mp,64);
  mpz_add_ui(mp,mp1,ps[0]);
  //mpz_out_str(NULL,10,mp);printf("\n");
  mpz_sub_ui(mp_1,mp,1);
  mpz_powm(res,mpz_2,mp_1,mp);
  return(mpz_cmp_ui(res,1)==0);
}

mpz_t mpow,mm,mm1,mp;

int pseudoprimep1(bigint p)
{
  bigint p1,ind;
  unsigned long int *ps;
  //printf("In pseudoprime with p=");print_bigint(p);printf("\n");
  ps=(unsigned long int *) &p;
  mpz_set_ui(mp,ps[1]);
  mpz_mul_2exp(mp1,mp,64);
  mpz_add_ui(mp,mp1,ps[0]);
  //mpz_out_str(NULL,10,mp);printf("\n");
  p1=p-1;
  mpz_set_ui(mm,1);
  mpz_set_ui(mpow,4);
  for(ind=2;ind<p;ind+=ind)
    {
      if(p1&ind)
	{
	  mpz_mul(mm1,mm,mpow);
	  mpz_fdiv_r(mm,mm1,mp);
	}
      mpz_mul(mm1,mpow,mpow);
      mpz_fdiv_r(mpow,mm1,mp);
    }
  return(mpz_cmp_ui(mm,1)==0);
}

int main(int argc, char **argv)
{
  long int sieve_num,num_its;
  ptype i,j,it,seg;
  bigint target_start,del2_sum,X0,p,mid;
  ptype count,segs_per_sieve=SEGS_PER_SIEVE;
  __int64_t del,del2,del_sum;
  FILE *outfile;

  mpz_init(mpz_2);
  mpz_init(res);
  mpz_init(mp);
  mpz_init(mp1);
  mpz_init(mp_1);
  mpz_set_ui(mpz_2,2);
  mpz_init(mm);
  mpz_init(mm1);
  mpz_init(mpow);

  if(argc!=4)
    {
      printf("argc was %ld\n",argc);
      print_usage();
    }
  sieve_num=atol(argv[1]);
  
  if(sieve_num<0)
    {
      printf("sieve_num was %ld\n",sieve_num);
      print_usage();
    }
  
  num_its=atol(argv[2]);
  if(num_its<1)
    {
      printf("num_its was %ld\n",num_its);
      print_usage();
    }
  if(!(outfile=fopen(argv[3],"wb")))
    fatal_error("Failed to open outfile for binary write.\n");

  fwrite(&num_its,sizeof(long int),1,outfile);
  fwrite(&segs_per_sieve,sizeof(ptype),1,outfile);

  for(i=0,X0=1;i<LOG_10_X0;i++,X0*=10);

  p=X0-XI*SEGS_PER_SIEVE*NUM_SIEVES/2;
  p+=XI*SEGS_PER_SIEVE*sieve_num+1;

  for(i=0,count=0;i<100000000;i+=30,p+=30)
    if(pseudoprimep(p))
      count++;
  printf("Count=%lu\n",count);
  exit(0);

  mid=p+XI/2-1;
  for(it=0;it<num_its;it++,target_start+=XI*SEGS_PER_SIEVE)
    {
      printf("Running sieve %lu starting at ",it+sieve_num);print_bigint(p);printf("\n");
      fwrite(&it,sizeof(ptype),1,outfile);
      fwrite(&p,sizeof(bigint),1,outfile);

      for(seg=0;seg<SEGS_PER_SIEVE;seg++,mid+=XI)
	{
	  printf("Segment %lu\n",seg);
	  printf("mid=");print_bigint(mid);printf("\n");
	  printf("p=");print_bigint(p);printf("\n");
	  for(count=0,del_sum=0,del2_sum=0,i=0;i<XI/2;i++,p+=2)
	    {
	      if(rem_80_48(p,3)==0)
		continue;
	      if(rem_80_48(p,5)==0)
		continue;
	      if(rem_80_48(p,7)==0)
		continue;
	      if(i%1000000==0)
		{
		  printf("Reached ");print_bigint(p);printf("\n");
		}
	      if(pseudoprimep(p))
		{
		  //print_bigint(p);printf(" is pseudoprime\n");
		  count++;
		  del=p-mid;
		  del_sum+=del;
		  del2=del*del;
		  del2_sum+=del2;
		}
	    }
	  fwrite(&count,sizeof(ptype),1,outfile);
	  fwrite(&del_sum,sizeof(__int64_t),1,outfile);
	  fwrite(&del2_sum,sizeof(bigint),1,outfile);
	  printf("we have %lu primes\nsigma (p-t0)=%ld\nsigma (p-t0)^2=",count,del_sum);
	  print_bigint(del2_sum);
	  printf("\n");
	}
    }
  return(0);
}
