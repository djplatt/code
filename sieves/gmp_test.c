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


mpz_t mod_tmp,mod_tmp1,mod_bit,mod_pow,mpz_2,mod_mp,mp_1;

void mpz_myinc(mpz_ptr x)
{
  mpz_add_ui(x,x,1);
}

inline int pseudoprimep1(mpz_ptr mp)
{
  mpz_sub_ui(mp_1,mp,1);
  mpz_set_ui(mod_tmp,1); // the result
  mpz_set_ui(mod_pow,4); // 2^bit
  mpz_fdiv_q_2exp(mod_bit,mp_1,1);
  while(mpz_cmp_ui(mod_bit,0)!=0)
    {
      if(mpz_odd_p(mod_bit))
	{
	  mpz_mul(mod_tmp1,mod_tmp,mod_pow);
	  mpz_fdiv_r(mod_tmp,mod_tmp1,mp);
	}
      mpz_fdiv_q_2exp(mod_bit,mod_bit,1); // shift right 1
      mpz_mul(mod_tmp1,mod_pow,mod_pow); // square
      mpz_fdiv_r(mod_pow,mod_tmp1,mp); // remainder
    }
  //printf("pseudoprimep computed ");mpz_out_str(NULL,10,mod_tmp);printf("\n");
  //exit(0);
  return(mpz_cmp_ui(mod_tmp,1)==0);
}

inline int pseudoprimep2(bigint p_mp)
{
  //printf("In pseudoprimep2 with ");print_bigint(p_mp);printf("\n");
  bigint mod_bit=p_mp>>1;
  //printf("In pseudoprimep2 with ");print_bigint(mod_bit);printf("\n");

  mpz_set_ui(mod_tmp,1); // the result
  mpz_set_ui(mod_pow,4); // 2^bit

  unsigned long int *ps;
  //printf("In pseudoprime with p=");print_bigint(p);printf("\n");
  ps=(unsigned long int *) &p_mp;
  mpz_set_ui(mod_mp,ps[1]);
  mpz_mul_2exp(mod_tmp1,mod_mp,64);
  mpz_add_ui(mod_mp,mod_tmp1,ps[0]);
  //printf("mod_mp set to ");mpz_out_str(NULL,10,mod_mp);printf("\n");
  while(mod_bit!=0)
    {
      //printf("mod_bit= ");print_bigint(mod_bit);printf("\n");
      //printf("res= ");mpz_out_str(NULL,10,mod_tmp);printf("\n");
      //printf("pow= ");mpz_out_str(NULL,10,mod_pow);printf("\n");
      if(mod_bit&1)
	{
	  mpz_mul(mod_tmp1,mod_tmp,mod_pow);   // *2^bit
	  mpz_fdiv_r(mod_tmp,mod_tmp1,mod_mp); // remainder
	}
      mod_bit>>=1;
      mpz_mul(mod_tmp1,mod_pow,mod_pow); // square
      mpz_fdiv_r(mod_pow,mod_tmp1,mod_mp); // remainder
    }
  //printf("pseudoprimep computed ");mpz_out_str(NULL,10,mod_tmp);printf("\n");
  //exit(0);
  return(mpz_cmp_ui(mod_tmp,1)==0);
}

inline int pseudoprimep(mpz_ptr mp)
{
  mpz_sub_ui(mp_1,mp,1);
  mpz_powm(mod_tmp,mpz_2,mp_1,mp);
  return(mpz_cmp_ui(mod_tmp,1)==0);
}

int main(int argc, char**argv)
{
  mpz_init(mod_tmp);
  mpz_init(mod_tmp1);
  mpz_init(mod_bit);
  mpz_init(mod_pow);
  mpz_init(mpz_2);
  mpz_set_ui(mpz_2,2);
  mpz_init(mp_1);

  long unsigned int i,count=0;
  
  mpz_t mp;
  mpz_init(mp);
  mpz_init(mp_1);

  mpz_set_ui(mp,1);
  for(i=0;i<atoi(argv[1]);i++)
    mpz_mul_ui(mp,mp,10);
  mpz_set(mp_1,mp);
  mpz_myinc(mp);
  printf("Running 10^6 pseudoprimep tests from ");mpz_out_str(NULL,10,mp);printf("\n");
  for(i=0;i<1000000;i++)
    {
      if(pseudoprimep(mp))
	count++;
      mpz_add_ui(mp,mp,30);
    }
    /*
  bigint mp=1;
  mpz_init(mod_mp);
  for(i=0;i<24;i++) mp*=10;
  mp++;
  printf("Running 10^6 pseudoprimep tests from ");print_bigint(mp);printf("\n");
  for(i=0;i<1000000;i++)
    {
      if(pseudoprimep2(mp))
	count++;
      //exit(0);
      mp+=30;
    }
    */
  printf("count=%lu\n",count);
  return(0);
}
