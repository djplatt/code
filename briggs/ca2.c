/*

Compute successive Collossally Abundant Numbers
and check Robin's inequality 
sigma(n)/n-exp(gamma)log log n

n is specified by its prime factorisation in the file given in command line
n must contain all primes < its largest prime factor
each line 
p1 p2 e
denotes the product p^e where p is each prime in [p1,p2]

Note that all CA numbers are of this form

See Briggs - Abundant Numbers and the RH - Exp Math 15(2) (2006)
D.J.Platt 2018
*/
  
#include "stdio.h"
#include "stdlib.h"
#include "primesieve.h" // use Kim Walisch's primesieve
#include "arb.h" // and Fredrik's arb
#include "pari.h"

#define false (1==0)
#define true (0==0)

// data structure for n
// p1 p2 e represents prod p in [p1,p2] p^e
// log_p1 is log(p1) so we can compute log_p easily
// eps is the epsilon in Briggs
typedef struct
{
  uint64_t p1;
  uint64_t p2;
  uint64_t e;
  arb_t log_p1;
  arb_t eps;
} entry_t;

// print one
void print_entry(entry_t *ent)
{
  printf("%lu %lu %lu ",ent->p1,ent->p2,ent->e);
  arb_printd(ent->log_p1,20);
  printf(" ");
  arb_printd(ent->eps,20);
  printf("\n");
}

#define MAX_ENTRIES (128)

typedef struct
{
  uint64_t n_entries;
  uint64_t next_entry;
  entry_t *entries;
} ca_t;


// compute a new epsilon for a new prime = log_p(1+1/p)
// Note Briggs has log_p(1+p) but this looks like a typo
void eps_ext(entry_t *ent, int64_t prec) // new prime (exponent == 0)
{
  static int init=false;
  static arb_t tmp,tmp1;
  if(!init)
    {
      init =true;
      arb_init(tmp);
      arb_init(tmp1);
    }
  /*
  arb_set_ui(tmp,ent->p1+1);
  arb_log(tmp1,tmp,prec);
  arb_div(ent->eps,tmp,ent->log_p1,prec);
  return;
  */


  arb_set_ui(tmp,ent->p1);
  arb_inv(tmp1,tmp,prec);
  arb_log1p(tmp,tmp1,prec);
  arb_div(ent->eps,tmp,ent->log_p1,prec);
}

void eps_inc(entry_t *ent, int64_t prec) // exponent == 1
{
  static int init=false;
  static arb_t tmp,tmp1;
  if(!init)
    {
      init =true;
      arb_init(tmp);
      arb_init(tmp1);
    }
  arb_set_ui(tmp,ent->p1+1);
  arb_mul_ui(tmp1,tmp,ent->p1,prec); // p^2+p
  arb_inv(tmp,tmp1,prec);
  arb_log1p(tmp1,tmp,prec);
  arb_div(ent->eps,tmp1,ent->log_p1,prec); //log_p((p+1+1/p)/(p+1))
}

void eps_max(entry_t *ent, int64_t prec) // exponent > 1
{
  static int init=false;
  static arb_t tmp,tmp1,tmp2;
  if(!init)
    {
      init =true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
    }
  //printf("in eps_max with ");print_entry(ent);
  arb_set_ui(tmp,ent->p1-1);
  arb_set_ui(tmp2,ent->p1);
  arb_pow_ui(tmp1,tmp2,ent->e+2,prec); // p^e+2
  arb_sub_ui(tmp2,tmp1,ent->p1,prec); // p^e+2-p
  arb_div(tmp1,tmp,tmp2,prec); // p-1/p^e+2-p
  arb_log1p(tmp2,tmp1,prec);
  arb_div(ent->eps,tmp2,ent->log_p1,prec);
}

// rhs contains log n
// lhs contains sigma(n)/n
void check_Robin(arb_t lhs, arb_t rhs, int64_t prec, int verbosep)
{
  static int init=false;
  static arb_t tmp1,tmp2,tmp3,egamma;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(egamma);
      arb_const_euler(tmp1,prec); // gamma
      arb_exp(egamma,tmp1,prec); // exp(gamma)
    }

  arb_log(tmp1,rhs,prec); // log log n
  arb_mul(tmp2,egamma,tmp1,prec); // exp(gamma) log log n
  arb_sub(tmp3,lhs,tmp2,prec);
  if(verbosep)
    {
      printf("log log n = ");arb_printd(tmp1,20);printf(" lhs-rhs=");arb_printd(tmp3,20);printf("\n");fflush(stdout);
    }
  if(!arb_is_negative(tmp3))
    {
      printf("Robin check failed.\n");
      arb_exp(tmp1,rhs,prec); // rhs contains log n, print n
      printf("n = ");arb_printd(tmp1,30);printf("\n");
      printf("sigma(n)/n= ");arb_printd(lhs,30);printf("\n");
      printf("exp(gamma) log log n= ");arb_printd(tmp2,30);printf("\n");
      printf("lhs-rhs=");
      arb_printd(tmp3,30);
      printf("\n");
      exit(0);
    }
}

// find the largest epsilon in the database
// simple linear search. Very few entries in practice
// but even so a heap might be quicker.
int64_t find_biggest(entry_t *entries, uint64_t next_entry, int64_t prec)
{
  static int init=false;
  static arb_t tmp,biggest;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(biggest);
    }
  int64_t ptr=0;
  arb_set(biggest,entries[0].eps);
  for(uint64_t i=1;i<=next_entry;i++) // check the dummy one at the end
    {
      arb_sub(tmp,biggest,entries[i].eps,prec);
      if(arb_is_negative(tmp))
	{
	  arb_set(biggest,entries[i].eps);
	  ptr=i;
	}
      else
	if(!arb_is_positive(tmp)) // not enough precision to tell eps's apart
	  {
	    printf("Floating point comparison failed. Use more bits. Exiting.\n");
	    exit(0);
	  }
    }
  return ptr;
}

// add a prime p with log log_p
// lhs before was sigma(n)/n, change it to sigma(np)/np
// rhs before was log n, change it to log pn
void include_prime(entry_t *entry, uint64_t p, arb_t log_p, arb_t lhs, arb_t rhs, int64_t prec)
{
  static int init=false;
  static arb_t tmp1,tmp2;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
    }

  entry->p2=p;
  arb_add(rhs,rhs,log_p,prec); // log(n*p) = log n + log p
  arb_set_ui(tmp1,p+1);
  arb_div_ui(tmp2,tmp1,p,prec);
  arb_mul(lhs,lhs,tmp2,prec); // sigma(n*p)/(n*p) =sigma(n)/n * (1+p)/p
}

// increase the exponent of the single prime p1
// update lhs and rhs
void add_exp(entry_t *entry, arb_t lhs, arb_t rhs, int64_t prec)
{
  static int init=false;
  static arb_t tmp,tmp1,tmp2;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
    }
  entry->e++;
  arb_add(rhs,rhs,entry->log_p1,prec);
  arb_set_ui(tmp,entry->p1);
  eps_max(entry,prec);
  arb_pow_ui(tmp1,tmp,entry->e+1,prec); // p^(e+2)
  arb_sub_ui(tmp2,tmp1,1,prec); // p^(e+2)-1
  arb_sub_ui(tmp,tmp1,entry->p1,prec); // p^(e+2)-p
  arb_mul(tmp1,lhs,tmp2,prec);
  arb_div(lhs,tmp1,tmp,prec);
}

// replace sigma(p^e)/p^e with sigma(p^e+1)/p^e+1
void arb_sig_del(arb_t rhs, uint64_t p, uint64_t e, int64_t prec)
{
  static int init=false;
  static arb_t tmp3,tmp1,tmp2;
  if(!init)
    {
      init=true;
      arb_init(tmp3);
      arb_init(tmp1);
      arb_init(tmp2);
    }
  arb_set_ui(tmp1,p);
  arb_pow_ui(tmp2,tmp1,e+2,prec);
  arb_sub_ui(tmp1,tmp2,1,prec);
  arb_sub_ui(tmp3,tmp2,p,prec);
  arb_mul(tmp2,rhs,tmp1,prec);
  arb_div(rhs,tmp2,tmp3,prec);
}

int main(int argc, char **argv)
{
  // echo the command line
  printf("Command line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  // check command line parameters
  if(argc!=4)
    {
      printf("Usage:- %s <filename> <prec> <n its>\n",argv[0]);
      return 0;
    }
  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Failed to open file %s for input. Exiting.\n",argv[1]);
      return 0;
    }
  pari_init(500000,0);
  
  uint64_t prec=atol(argv[2]); // working prec for arb.
  uint64_t n_its=atol(argv[3]); // how many iterations

  // initialise database
  ca_t ca;
  ca.entries=(entry_t *)malloc(sizeof(entry_t)*MAX_ENTRIES);
  ca.next_entry=0;
  ca.n_entries=MAX_ENTRIES;
  for(uint64_t i=0;i<MAX_ENTRIES;i++)
    {
      arb_init(ca.entries[i].eps);  
      arb_init(ca.entries[i].log_p1);
    }
  uint64_t p1,p2,e;

  // setup primesieve
  primesieve_iterator it;
  primesieve_init(&it);
  
  arb_t pe,pe1,lhs,tmp1,tmp2,tmp3,rhs;
  arb_init(pe);arb_init(pe1);arb_init(lhs);arb_init(rhs);
  arb_set_ui(lhs,1); // will hold prod sigma(p^e)/p^e
  arb_set_ui(rhs,0); // a nop, but makes me happy. will hold log n
  arb_init(tmp1);arb_init(tmp2);arb_init(tmp3);
  uint64_t tp=primesieve_next_prime(&it); // 2

  // read in the entries for the database from the file
  while(fscanf(infile,"%lu %lu %lu\n",&p1,&p2,&e)==3)
    {
      if(ca.next_entry==ca.n_entries) // database full
	{
	  printf("Ran out of room. Exiting.\n");
	  exit(0);
	}
      ca.entries[ca.next_entry].p1=p1;
      arb_log_ui(ca.entries[ca.next_entry].log_p1,p1,prec);
      ca.entries[ca.next_entry].p2=p2;
      ca.entries[ca.next_entry].e=e;
      ca.next_entry++;
      
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
	  tp=primesieve_next_prime(&it);
	}
    }
  if(ca.next_entry==ca.n_entries) // datbase full
    {
      printf("Ran out of room. Exiting.\n");
      exit(0);
    }

  // add dummy entry at end, p 0 0, where p is next prime
  ca.entries[ca.next_entry].p1=tp;
  ca.entries[ca.next_entry].p2=0;
  ca.entries[ca.next_entry].e=0;
  arb_log_ui(ca.entries[ca.next_entry].log_p1,tp,prec);

  eps_ext(ca.entries+ca.next_entry,prec); // compute epsilon for dummy
  eps_inc(ca.entries+ca.next_entry-1,prec); // compute epsilon for last
  for(uint64_t i=0;i<ca.next_entry-1;i++) // compute other epsilons
    eps_max(ca.entries+i,prec);
  ca.next_entry++;

  uint64_t report_its=n_its/100; // report every 1%
  if(report_its==0) report_its=1;
  // entries now set up, rhs contains log n, lhs constains sigma(n)/n
  for(uint64_t n=0;n<n_its;n++)
    {
      if((n%report_its)==0)
	{printf("%lu ",n);check_Robin(lhs,rhs,prec,1);} // verbose check
      else
	check_Robin(lhs,rhs,prec,0); // silent check
      uint64_t ptr=find_biggest(ca.entries,ca.next_entry,prec); 

      if(ca.entries[ptr].p1==ca.entries[ptr].p2) // just increment exponent
	{
	  add_exp(ca.entries+ptr,lhs,rhs,prec);
	  continue;
	}
      if(ptr==ca.next_entry-1) // need a new prime, add it to the last proper entry
	{

	  include_prime(ca.entries+ptr-1,ca.entries[ptr].p1,ca.entries[ptr].log_p1,lhs,rhs,prec);
	  tp=primesieve_next_prime(&it);
	  ca.entries[ptr].p1=tp;
	  ca.entries[ptr].p2=0;
	  ca.entries[ptr].e=0;
	  arb_log_ui(ca.entries[ptr].log_p1,tp,prec);
	  eps_ext(ca.entries+ptr,prec);
	  continue;
	}
      // going to split a prime off the start and increase its exponent
      // is the previous entry of the right exponent?
      if(ca.entries[ptr-1].e==ca.entries[ptr].e+1)
	{
	  arb_add(rhs,rhs,ca.entries[ptr].log_p1,prec); // increase log(n)
	  arb_sig_del(lhs,ca.entries[ptr].p1,ca.entries[ptr].e,prec);
	  ca.entries[ptr-1].p2=ca.entries[ptr].p1;
	  // we now want the next prime after entries[ptr].p1
	  // but the primesieve iterator could be pointing anywhere
	  // so we use pari's nextprime instead
	  // strictly, this returns a pseudoprime but there are no
	  // composite pseudoprimes < 2^64

	  ca.entries[ptr].p1=unextprime(ca.entries[ptr].p1+1); // call pari to do nextprime
	  arb_log_ui(ca.entries[ptr].log_p1,ca.entries[ptr].p1,prec);
	  if(ca.entries[ptr].e==1)
	    eps_inc(ca.entries+ptr,prec);
	  else
	    eps_max(ca.entries+ptr,prec);
	}
      else
	{
	  // first check there is room
	  if(ca.next_entry==ca.n_entries) // full
	    {
	      printf("Ran out of room. Exiting.\n");
	      exit(0);
	    }
	  // move everything up one to accomodate it.
	  for(uint64_t ptr1=ca.next_entry-1;ptr1>ptr;ptr1--)
	    {
	      ca.entries[ptr1+1].p1=ca.entries[ptr1].p1;
	      ca.entries[ptr1+1].p2=ca.entries[ptr1].p2;
	      ca.entries[ptr1+1].e=ca.entries[ptr1].e;
	      arb_swap(ca.entries[ptr1+1].log_p1,ca.entries[ptr1].log_p1);
	      arb_swap(ca.entries[ptr1+1].eps,ca.entries[ptr1].eps);
	    }
	  // now remove the first prime from the original record
	  //char buff[1024];
	  //sprintf(buff,"nextprime(%lu)",ca.entries[ptr].p1+1);
	  ca.entries[ptr+1].p1=unextprime(ca.entries[ptr].p1+1);
	  ca.entries[ptr+1].p2=ca.entries[ptr].p2;
	  ca.entries[ptr+1].e=ca.entries[ptr].e;
	  arb_log_ui(ca.entries[ptr+1].log_p1,ca.entries[ptr+1].p1,prec);
	  if(ca.entries[ptr+1].e==1)
	    eps_inc(ca.entries+ptr+1,prec);
	  else
	    eps_max(ca.entries+ptr+1,prec);
	  ca.entries[ptr].p2=ca.entries[ptr].p1;
	  eps_max(ca.entries+ptr,prec);
	  add_exp(ca.entries+ptr,lhs,rhs,prec);
	  ca.next_entry++;
	}      
    }
  printf("n=");
  arb_exp(tmp1,rhs,prec);
  arb_printd(tmp1,20);
  printf("\n");
  for(uint64_t i=0;i<ca.next_entry-1;i++)
    printf("%lu %lu %lu\n",ca.entries[i].p1,ca.entries[i].p2,ca.entries[i].e);

  pari_close();
  
  return 0;
}

