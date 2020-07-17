/*

File: sieve1.0.c

Created: 8 April 2008

Version: 1.0

Last Modified: 8 April 2008

Dialect: C

Requires: GMP  v. 4.2.2
          MPFR v. 2.3.1
          MPFI v. 1.3.4-RC3

Implementation notes:

V 1.0 Initial implementation

Build instructions: gcc -osieve1.0 sieve1.0.c -O2 -lmpfi -lmpfr -lgmp

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */


/*

sieves [floor(x/width),ceil(x*width)] for primes and prime powers
for u=p^n, sums 1/n*phi(u,x,lambda) where 
phi(u,x,lambda)=1/2(1+erf(ln(u/x)/(lambda*sqrt(2)))) (u<x)
               =1/2(-1+erf(ln(u/x)/(lambda*sqrt(2)))) otherwise

*/

#include "stdio.h"
#include "math.h"
#include "gmp.h"
#include "mpfr.h"
#include "mpfi.h"
#include "mpfi_io.h"

#define SUCCESS 0
#define QUIT_CODE 0
#define FAILURE 1
#define TRUE (1==1)
#define FALSE (0==1)

void mpfi_print(mpfi_ptr w)
{
    mpfi_out_str(stdout,10,0,w);
    printf("\n");
};
void mpfr_print(mpfr_ptr w)
{
    mpfr_out_str(stdout,10,0,w,GMP_RNDN);
    printf("\n");
};

void mpz_print(mpz_ptr w)
{
    mpz_out_str(stdout,10,w);
    printf("\n");
};

mpz_t x,sieve_start,sieve_end;
mpfr_t width,erf_tmp1,erf_tmp2;
unsigned long int prime_len, char_len,sieve_char_len,sieve_len;



void sieve_setup(int prec)
{
   mpfr_set_default_prec(prec);
   mpfr_init(erf_tmp1);
   mpfr_init(erf_tmp2);
};

unsigned char *prime_array;
unsigned char *sieve_array;

unsigned char pattern(unsigned char pos)
{
   switch(pos)
   {
      case 0: return(254);
      case 1: return(253);
      case 2: return(251);
      case 3: return(247);
      case 4: return(239);
      case 5: return(223);
      case 6: return(191);
      case 7: return(127);
      default: printf("Fatal error in pattern.\n");
               exit(SUCCESS);
   };
};

int low_prime(unsigned long int i)
{
   unsigned long int c;
   unsigned char pos;
   c=i>>3;
   pos=i&7;
   return(prime_array[c]&(~pattern(pos)));
};


unsigned long int next_low_prime(unsigned long int i)
{
   for(i++;i<=prime_len;i++)
      if(low_prime(i))
         return(i);
   return(0);
};

void clear_low_prime(unsigned long int i)
{
   unsigned long int c;
   unsigned char pos;
   c=i>>3;
   pos=i&7;
   prime_array[c]=prime_array[c]&pattern(pos);
};

void clear_sieve_prime(unsigned long int char_pos, unsigned char bit_pos)
{
   sieve_array[char_pos]=sieve_array[char_pos]&pattern(bit_pos);
};

int sieve_prime(unsigned long int char_pos, unsigned char bit_pos)
{
   return(sieve_array[char_pos]&(~pattern(bit_pos)));
};


int build_sieves()
{
   unsigned long int i,j,sqrt_len,char_offset;
   unsigned long int char_del;
   unsigned char bit_offset,bit_del;
   mpz_t tmp;
   
   printf("Allocating Prime Array for 1..%Lu.\n",prime_len); 
   char_len=(prime_len>>3)+1;
   printf("Using %Lu bytes of storage.\n",char_len);
   prime_array=malloc(char_len);
   if(prime_array==NULL)
   {
      printf("Error allocating memory for prime array.\n");
      return(FAILURE);
   };

   printf("Allocating Sieve Array of width %Lu.\n",sieve_len); 
   sieve_char_len=(sieve_len>>3)+1;
   printf("Using %Lu bytes of storage.\n",sieve_char_len);
   sieve_array=malloc(sieve_char_len);
   if(sieve_array==NULL)
   {
      printf("Error allocating memory for sieve array.\n");
      return(FAILURE);
   };
   




   for(i=1;i<char_len;i++)
      prime_array[i]=255;
   prime_array[0]=254;
   sqrt_len=ceil(sqrt(prime_len));
   printf("Sqrt(len) set to %Lu.\n",sqrt_len);
   printf("Primes <= %Lu.\n",prime_len);
   for(i=2;i<=sqrt_len;i++)
      if(low_prime(i))
         for(j=i*i;j<=prime_len;j=j+i)
            clear_low_prime(j);
//   for(i=2;i!=0;i=next_low_prime(i))
//      printf("%Lu.\n",i);

   for(i=0;i<sieve_char_len;i++)
      sieve_array[i]=255;

   for(i=2;i!=0;i=next_low_prime(i))
   {
      char_offset=mpz_fdiv_ui(sieve_start,i);
      if(char_offset!=0)
         char_offset=i-char_offset;
//      printf("Remainder %Lu ",char_offset);
      bit_offset=char_offset&7;
      char_offset=char_offset>>3;
      char_del=i>>3;
      bit_del=i&7;
//      printf("For prime %Lu, Char Offset %Lu, Bit Offset %Lu, Char Del %Lu, Bit Del %Lu.\n",
//             i,char_offset,bit_offset,char_del,bit_del);
      while(char_offset<sieve_char_len)
      {
         clear_sieve_prime(char_offset,bit_offset);
         char_offset=char_offset+char_del;
         bit_offset=bit_offset+bit_del;
         if (bit_offset>7)
         {
            bit_offset=bit_offset&7;
            char_offset++;
         };
      };
   };
/*
   printf("Prime in sieve range.\n");
   bit_offset=0;
   mpz_init(tmp);
   mpz_set(tmp,sieve_start);
   for(char_offset=0;char_offset<sieve_char_len;)
   {
      if(sieve_prime(char_offset,bit_offset))
         mpz_print(tmp);
      bit_offset++;
      if(bit_offset==8)
      {
         bit_offset=0;
         char_offset++;
      };
      mpz_add_ui(tmp,tmp,1);
   };     
*/
   return(SUCCESS);
};


void mpfi_erf(mpfi_ptr res, mpfi_ptr op1)
{

// a nice strictly increasing function so nothing clever required
   
      mpfi_get_left(erf_tmp1,op1);
      mpfr_erf(erf_tmp1,erf_tmp1,GMP_RNDD);
      mpfi_get_right(erf_tmp2,op1);
      mpfr_erf(erf_tmp2,erf_tmp2,GMP_RNDU);
      mpfi_interv_fr(res,erf_tmp1,erf_tmp2);   
};
  

void phi(mpfi_ptr res,mpz_ptr u,mpz_ptr x, mpfi_ptr root_2_lambda)
{
      mpfi_set_z(res,u);
      mpfi_div_z(res,res,x);
      mpfi_log(res,res);
      mpfi_div(res,res,root_2_lambda);
      mpfi_erf(res,res);
      if(mpz_cmp(u,x)>0)    // u>x
      {
         mpfi_sub_ui(res,res,1);
         mpfi_div_ui(res,res,2);
      }
      else
      {
         mpfi_add_ui(res,res,1);
         mpfi_div_ui(res,res,2);
      };
//      printf("Phi of ");mpz_print(u);printf("   returning ");mpfi_print(res);
};

void do_sieve (mpfi_ptr res, mpz_ptr x, mpfi_ptr lambda)
{
   mpz_t p;
   mpfi_t p_phi,root_2_lambda;
   double lnx;
   unsigned char bit_offset;
   unsigned long int char_offset,n,pr;

   mpz_init(p);
   mpfi_init(p_phi);
   mpfi_init(root_2_lambda);
   mpfi_set_ui(res,0);
   mpfi_set_ui(root_2_lambda,2);
   mpfi_sqrt(root_2_lambda,root_2_lambda);
   mpfi_mul(root_2_lambda,root_2_lambda,lambda);
   char_offset=0;
   bit_offset=0;
   mpz_set(p,sieve_start);
   while(mpz_cmp(sieve_end,p)>=0)
   {
      if(sieve_prime(char_offset,bit_offset))
      {
         phi(p_phi,p,x,root_2_lambda);
         mpfi_add(res,res,p_phi);
      };
      bit_offset++;
      if(bit_offset==8)
      {
         bit_offset=0;
         char_offset++;
      };
      mpz_add_ui(p,p,1);
   };

   lnx=log(mpz_get_d(sieve_start));
   pr=2;
   while(pr!=0)
   {
      n=ceil(lnx/log(pr));
      mpz_set_ui(p,pr);
      mpz_pow_ui(p,p,n);
      while(mpz_cmp(sieve_end,p)>=0)
      {
         phi(p_phi,p,x,root_2_lambda);
//         printf("found prime power at ");mpz_print(p);
//         printf("phi returned ");mpfi_print(p_phi);
         mpfi_div_ui(p_phi,p_phi,n);
         mpfi_add(res,res,p_phi);
         mpz_mul_ui(p,p,pr);
         n++;
      };
   pr=next_low_prime(pr);
   };
};

void print_usage()
{
   printf("Usage: sieve1.0 <prec> <x> <lam_power> <width>\n");
   printf("<prec> an unsigned integer in [%d,%d], the number of bits of precsion to use.\n",MPFR_PREC_MIN,MPFR_PREC_MAX);
   printf("<x> an arbitary precision unsigned integer > 2.\n");
   printf("<lam_power> an unsigned long integer, lambda<-2^(-lam_power).\n");
   printf("<width> an arbitary precision float > 1.0.\n");
};


int main(int argc, char **argv)
{

   unsigned long int prec,lam_power,i;
   mpfr_t tmp1,tmp2;
   mpz_t tmp3;
   mpfi_t lambda,res;

/*  check all the command line arguments are ok, if not print message
    and exit sharpish */


    if(argc!=5)
    {
	print_usage();
	return(QUIT_CODE);
    };

    prec=atoi(argv[1]);
    if((prec<MPFR_PREC_MIN)||(prec>MPFR_PREC_MAX))
    {
	print_usage();
	return(QUIT_CODE);
    };
    
    sieve_setup(prec);

    mpz_init(x);
    if(mpz_set_str(x,argv[2],10)!=SUCCESS)
    {
	print_usage();
	return(QUIT_CODE);
    };
    if(mpz_cmp_ui(x,2)<0)
    {
	print_usage();
	return(QUIT_CODE);
    };


    mpfi_init(lambda);
    mpfi_set_ui(lambda,1);
    lam_power=atoi(argv[3]);
    for(i=0;i<lam_power;i++)
        mpfi_div_ui(lambda,lambda,2);


    mpfr_init(width);
    if(mpfr_set_str(width,argv[4],10,GMP_RNDN)!=SUCCESS)
      {
	print_usage();
	return(QUIT_CODE);
      };
    if((mpfr_sgn(width)<=0)||(mpfr_cmp_d(width,1.0)<=0))
      {
	print_usage();
	return(QUIT_CODE);
      };
    printf("x set to ");
    mpz_print(x);
    printf("lambda set to\n");
    mpfi_print(lambda);
    printf("width set to ");
    mpfr_print(width);
/* call standard mpfi routine to assign ln2:=ln(2) */


    mpfr_init(tmp1);
    mpfr_set_z(tmp1,x,GMP_RNDU);
    mpfr_init(tmp2);
    mpfr_set(tmp2,tmp1,GMP_RNDN);
    mpfr_mul(tmp1,tmp1,width,GMP_RNDU);
    mpz_init(sieve_end);
    mpz_init(sieve_start);
    
    mpfr_get_z(sieve_end,tmp1,GMP_RNDU);
    mpfr_div(tmp2,tmp2,width,GMP_RNDD);
    mpfr_get_z(sieve_start,tmp2,GMP_RNDD);
    mpfr_sqrt(tmp1,tmp1,GMP_RNDU);
    prime_len=mpfr_get_ui(tmp1,GMP_RNDU);
    mpz_init(tmp3);
    mpz_sub(tmp3,sieve_start,sieve_end);
    sieve_len=mpz_get_ui(tmp3)+1;
    if(build_sieves()!=SUCCESS)
        return(SUCCESS);
    mpfi_init(res);
    do_sieve(res,x,lambda);
    printf("Result of sieve ");mpfi_print(res);

    return(SUCCESS);
}


