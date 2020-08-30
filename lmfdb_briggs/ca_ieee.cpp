/*

Compute Robin's inequality sigma(n)/n-exp(gamma)log log n

n is specified by its prime factoriastion in the file given in command line
n must contain all primes < its largest prime factor
each line 
p1 p2 e
denote the product p^e where p is each prime in [p1,p2]

Note that all CA number are of this form

*/
  
#include "stdio.h"
#include "stdlib.h"
#include "primesieve.h" // use Kim Walisch's primesieve
#include "../includes/int_double14.0.h"

int main(int argc, char **argv)
{
  // echo the command line
  printf("Command line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  // check command line parameters
  if(argc!=2)
    {
      printf("Usage:- %s <filename>\n",argv[0]);
      return 0;
    }
  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Failed to open file %s for input. Exiting.\n",argv[1]);
      return 0;
    }

  uint64_t p1,p2,e;
  primesieve_iterator it;
  primesieve_init(&it);
  int_double pe,pe1,lhs,tmp1,tmp2,tmp3,rhs;
  _fpu_rndd();
  lhs=1.0;rhs=0.0;
  uint64_t tp=primesieve_next_prime(&it); // 2

  while(fscanf(infile,"%lu %lu %lu\n",&p1,&p2,&e)==3)
    {
      if(tp!=p1) // check the prime read is as expected
	{
	  printf("Error, we have missed a prime somewhere. Exiting.\n");
	  return 0;
	}
      printf("Processing record %lu %lu %lu\n",p1,p2,e);
      while(tp<=p2)
	{
	  pe=pow(tp,e);
	  pe1=pe*tp;
	  lhs*=(pe1-1)/(pe1-pe);
	  rhs+=e*log(int_double(tp));
	  //printf("p= %lu lhs=",tp);arb_printd(lhs,30);printf("\n");
	  tp=primesieve_next_prime(&it);
	}
    }
  rhs=exp(d_gamma)*log(rhs);
  print_int_double_str("lhs = ",lhs);
  print_int_double_str("rhs = ",rhs);
  print_int_double_str("lhs-rhs = ",lhs-rhs);
  return 0;
}

