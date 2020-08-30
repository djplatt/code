/*

Compute Robin's inequality sigma(n)/n-exp(gamma)log log n

n is specified by its prime factorisation in the file given in command line
n must contain all primes < its largest prime factor
each line 
p1 p2 e
denotes the product p^e where p is each prime in [p1,p2]

Note that all CA numbers are of this form

*/
  
#include "stdio.h"
#include "stdlib.h"
#include "primesieve.h" // use Kim Walisch's primesieve
#include "arb.h" // and Fredrik's arb

int main(int argc, char **argv)
{
  // echo the command line
  printf("Command line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  // check command line parameters
  if(argc!=3)
    {
      printf("Usage:- %s <filename> <prec>\n",argv[0]);
      return 0;
    }
  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Failed to open file %s for input. Exiting.\n",argv[1]);
      return 0;
    }
  uint64_t prec=atol(argv[2]); // working prec for arb.

  uint64_t p1,p2,e;
  primesieve_iterator it;
  primesieve_init(&it);
  arb_t pe,pe1,lhs,tmp1,tmp2,tmp3,rhs;
  arb_init(pe);arb_init(pe1);arb_init(lhs);arb_init(rhs);
  arb_set_ui(lhs,1); // will hold prod sigma(p^e)/p^e
  arb_set_ui(rhs,0); // a nop, but makes me happy. will hold log n
  arb_init(tmp1);arb_init(tmp2);arb_init(tmp3);
  uint64_t tp=primesieve_next_prime(&it); // 2

  while(fscanf(infile,"%lu %lu %lu\n",&p1,&p2,&e)==3)
    {
      if(tp!=p1) // check the prime read is as expected
	{
	  printf("Error, we have missed a prime somewhere. Exiting.\n");
	  return 0;
	}
      printf("Processing record %lu %lu %lu\n",p1,p2,e);
      while(tp<=p2) // run from p1 to p2 inclusive
	{
	  arb_ui_pow_ui(pe,tp,e,prec);
	  arb_mul_ui(pe1,pe,tp,prec);
	  arb_sub(tmp1,pe1,pe,prec); // p^(e+1)-p^e
	  arb_sub_ui(tmp2,pe1,1,prec); // p^(e+1)-1
	  arb_div(tmp3,tmp2,tmp1,prec); // [p^(e+1)-1]/[p^(e+1)-p^e]
	  arb_mul(lhs,lhs,tmp3,prec);
	  arb_set_ui(tmp2,tp);
	  arb_log(tmp1,tmp2,prec); // log(p)
	  arb_mul_ui(tmp2,tmp1,e,prec); // log(p^e)
	  arb_add(rhs,rhs,tmp2,prec);
	  //printf("p= %lu lhs=",tp);arb_printd(lhs,30);printf("\n");
	  tp=primesieve_next_prime(&it);
	}
    }

  arb_exp(tmp1,rhs,prec); // rhs contains log n, print n
  printf("n = ");arb_printd(tmp1,30);printf("\n");

  arb_log(tmp1,rhs,prec); // log log n
  arb_const_euler(tmp2,prec); // gamma
  arb_exp(tmp3,tmp2,prec); // exp(gamma)
  arb_mul(rhs,tmp3,tmp1,prec); // exp(gamma) log log n
  printf("lhs=");arb_printd(lhs,30);printf("\n");
  printf("rhs=");arb_printd(rhs,30);printf("\n");
  arb_sub(tmp1,lhs,rhs,prec);
  printf("lhs-rhs=");arb_printd(tmp1,30);printf("\n");

  return 0;
}

