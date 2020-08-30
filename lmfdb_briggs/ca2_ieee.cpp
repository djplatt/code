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

Uses int_double package

D.J.Platt 2018
*/
  
#include "stdio.h"
#include "stdlib.h"
#include "primesieve.h" // use Kim Walisch's primesieve
#undef ulong
#include "pari/pari.h"
#include "pari.c"
#undef ulong
#include "../includes/int_double14.0.h"


// data structure for n
// p1 p2 e represents prod p in[p1,p2] p^e
// log_p1 is log(p1) so we can compute log_p easily
// eps is the epsilon in Briggs
typedef struct
{
  uint64_t p1;
  uint64_t p2;
  uint64_t e;
  int_double log_p1;
  int_double eps;
} entry_t;

// print one
void print_entry(entry_t *ent)
{
  printf("%lu %lu %lu ",ent->p1,ent->p2,ent->e);
  print_int_double(ent->log_p1);
  printf(" ");
  print_int_double(ent->eps);
  printf("\n");
}

#define MAX_ENTRIES (128)

typedef struct
{
  uint64_t n_entries;
  uint64_t next_entry;
  entry_t *entries;
} ca_t;

#define STACK_LEN (128)
int_double log_stack[STACK_LEN];
uint64_t log_ptr=0;
int_double sig_stack[STACK_LEN];
uint64_t sig_ptr=0;

int_double add_item(int_double *stack, int_double item, uint64_t k, uint64_t *ptr)
{
  stack[ptr[0]++]=item;
  while(k&1)
    {
      stack[ptr[0]-2]+=stack[ptr[0]-1];
      ptr[0]--;
      k>>=1;
    }
  int_double res=stack[0];
  for(uint64_t p=1;p<ptr[0];p++)
    res+=stack[p];
  return res;
}

// compute a new epsilon for a new prime = log(1+1/p)/log(p)
void eps_ext(entry_t *ent) // new prime (exponent == 0)
{
  static int init=false;
  ent->eps=log1p(int_double(1)/ent->p1)/ent->log_p1;
}

void eps_inc(entry_t *ent) // exponent == 1
{
  uint64_t p=ent->p1;
  int_double tmp=p+1;
  tmp*=p;
  ent->eps=log1p(int_double(1)/tmp)/ent->log_p1;
}

void eps_max(entry_t *ent) // exponent > 1
{
  uint64_t p=ent->p1;
  int_double pe=exp((ent->e+2)*ent->log_p1); // p^{e+2}
  int_double pep=pe-p; // p^{e+2}-p
  ent->eps=log1p(int_double(p-1)/pep)/ent->log_p1;
  return;
}

// rhs contains log n
// lhs contains sigma(n)/n
int check_Robin(int_double lhs, int_double rhs, int verbosep)
{
  static bool init=false;
  static int_double exp_gam;
  if(!init)
    {
      init=true;
      exp_gam=exp(d_gamma);
    }
  int_double lln=log(rhs); // log log n
  int_double eglln=lln*exp_gam;
  int_double del=lhs-eglln;
  if(verbosep)
    {printf("log log n = ");print_int_double(lln);
      printf(" exp(gamma) log log n = ");print_int_double(eglln);
      printf("sigma(n)/n = ");print_int_double(lhs);
      printf(" lhs-rhs=");print_int_double(del);printf("\n");fflush(stdout);}
  if(del.right<=0.0)
    {
      printf("Robin check failed.\n");
      printf("log log n = ");print_int_double(lln);
      printf(" exp(gamma) log log n = ");print_int_double(eglln);
      printf("sigma(n)/n = ");print_int_double(lhs);
      printf(" lhs-rhs=");print_int_double(del);printf("\n");fflush(stdout);
      return false;
    }
  return true;
}

// find the largest epsilon in the database
int64_t find_biggest(entry_t *entries, uint64_t next_entry)
{
  int64_t ptr=0;
  int_double biggest=entries[0].eps;
  for(uint64_t i=1;i<=next_entry;i++) // check the dummy one at the end
    {
      int_double del=biggest-entries[i].eps;
      if(del.right>0.0)
	{
	  biggest=entries[i].eps;
	  ptr=i;
	}
      else
	if(del.left<=0.0)
	  {
	    printf("Floating point comparison failed at i=%lu. Use more bits. Exiting.\n",i);
	    for(uint64_t j=0;j<next_entry;j++)
	      print_entry(entries+j);
	    exit(0);
	  }
    }
  return ptr;
}

// add a brand new prime p with log log_p
// lhs before was sigma(n)/n, change it to sigma(np)/np
// rhs before was log n, change it to log pn
void include_prime(entry_t *entry, uint64_t p, int_double log_p, int_double &lhs, int_double &rhs, uint64_t k)
{
  entry->p2=p;
  rhs=add_item(log_stack,log_p,k,&log_ptr); // log(n*p) = log n + log p
  lhs=add_item(sig_stack,lhs/p,k,&sig_ptr); // sigma(np)/np = sigma(n)/n*(p+1)/p
}

// increase the exponent of the single prime p1
// update lhs and rhs
void add_exp(entry_t *entry, int_double &lhs, int_double &rhs, uint64_t k)
{
  entry->e++;
  rhs=add_item(log_stack,entry->log_p1,k,&log_ptr);
  eps_max(entry);
  int_double p=entry->p1;
  lhs=add_item(sig_stack,lhs*(p-1)/(exp((entry->e+1)*entry->log_p1)-p),k,&sig_ptr);
}

// replace sigma(p^e)/p^e with sigma(p^e+1)/p^e+1
void arb_sig_del(int_double &lhs, uint64_t p, uint64_t e, uint64_t k)
{
  int_double pp=p;
  lhs=add_item(sig_stack,lhs*(pp-1)/(pow(pp,e+2)-pp),k,&sig_ptr);
}

int main(int argc, char **argv)
{
  // set sse registers to round down
  // needed for int_double
  _fpu_rndd();
  // echo the command line
  printf("Command line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  // check command line parameters
  if(argc!=3)
    {
      printf("Usage:- %s <filename> <n its>\n",argv[0]);
      return 0;
    }
  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    {
      printf("Failed to open file %s for input. Exiting.\n",argv[1]);
      return 0;
    }
  uint64_t n_its=atol(argv[2]); // how many iterations


  // initialise database
  ca_t ca;
  ca.entries=(entry_t *)malloc(sizeof(entry_t)*MAX_ENTRIES);
  ca.next_entry=0;
  ca.n_entries=MAX_ENTRIES;

  uint64_t p1,p2,e;

  // setup primesieve
  primesieve_iterator it;
  primesieve_init(&it);

  int_double pe,pe1,lhs,tmp1,tmp2,tmp3,rhs;
  lhs=1;rhs=0;
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
      ca.entries[ca.next_entry].log_p1=log(p1);
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
	  pe=pow(int_double(tp),e);
	  pe1=pe*tp;
	  tmp1=pe1-pe;
	  tmp2=pe1-1;
	  tmp3=tmp2/tmp1;
	  lhs*=tmp3;
	  tmp2=tp;
	  tmp1=log(tmp2);
	  tmp2=tmp1*e;
	  rhs+=tmp2;
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
  ca.entries[ca.next_entry].log_p1=log(int_double(tp));
  eps_ext(ca.entries+ca.next_entry); // compute epsilon for dummy record
  eps_inc(ca.entries+ca.next_entry-1); // compute epsilon for last record
  for(uint64_t i=0;i<ca.next_entry-1;i++) // compute other epsilons
    eps_max(ca.entries+i);
  ca.next_entry++;

  for(uint64_t i=0;i<ca.next_entry;i++)
    print_entry(ca.entries+i);

  uint64_t report_its=n_its/1000; // report every 0.1%
  if(report_its==0) report_its=1;
  // entries now set up, rhs contains log n, lhs constains sigma(n)/n
  int rc;uint64_t ptr;
  log_stack[0]=rhs;
  log_ptr=1;
  sig_stack[0]=lhs;
  sig_ptr=1;
  for(uint64_t n=1;n<=n_its;n++)
    {
      /*
      for(uint64_t i=0;i<ca.next_entry;i++)
	print_entry(ca.entries+i);
      */
      if((n%report_its)==0)
	{printf("%lu ",n);rc=check_Robin(lhs,rhs,1);
	  for(uint64_t i=0;i<ca.next_entry-1;i++)
	    printf("%lu %lu %lu\n",ca.entries[i].p1,ca.entries[i].p2,ca.entries[i].e);
	} // verbose check
      else
	rc=check_Robin(lhs,rhs,0); // silent check
      if(!rc)
	{
	  ca.entries[ptr].e--;
	  break;
	}

      ptr=find_biggest(ca.entries,ca.next_entry); 

      if(ca.entries[ptr].p1==ca.entries[ptr].p2) // just increment exponent
	{
	  add_exp(ca.entries+ptr,lhs,rhs,n);
	  continue;
	}
      if(ptr==ca.next_entry-1) // need a new prime, add it to the last proper entry
	{

	  include_prime(ca.entries+ptr-1,ca.entries[ptr].p1,ca.entries[ptr].log_p1,lhs,rhs,n);
	  tp=primesieve_next_prime(&it);
	  ca.entries[ptr].p1=tp;
	  ca.entries[ptr].p2=0;
	  ca.entries[ptr].e=0;
	  ca.entries[ptr].log_p1=log(int_double(tp));
	  eps_ext(ca.entries+ptr);
	  continue;
	}
      // going to split a prime off the start and increase its exponent
      // is the previous entry of the right exponent?
      if(ca.entries[ptr-1].e==ca.entries[ptr].e+1)
	{
	  /*
	  print_int_double_str("Splitting.\nrhs was ",rhs);
	  print_entry(ca.entries+ptr-1);
	  print_entry(ca.entries+ptr);
	  print_int_double_str("log p1 = ",ca.entries[ptr].log_p1);
	  */
	  rhs=add_item(log_stack,ca.entries[ptr].log_p1,n,&log_ptr); // increase log(n)
	  //print_int_double_str("rhs now ",rhs);
	  arb_sig_del(lhs,ca.entries[ptr].p1,ca.entries[ptr].e,n);
	  ca.entries[ptr-1].p2=ca.entries[ptr].p1;
	  char buff[1024];
	  sprintf(buff,"nextprime(%lu)",ca.entries[ptr].p1+1);
	  ca.entries[ptr].p1=parilong(buff);
	  ca.entries[ptr].log_p1=log(int_double(ca.entries[ptr].p1));
	  if(ca.entries[ptr].e==1)
	    eps_inc(ca.entries+ptr);
	  else
	    eps_max(ca.entries+ptr);
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
	      ca.entries[ptr1+1].log_p1=ca.entries[ptr1].log_p1;
	      ca.entries[ptr1+1].eps=ca.entries[ptr1].eps;
	    }
	  // now remove the first prime from the original record
	  char buff[1024];
	  sprintf(buff,"nextprime(%lu)",ca.entries[ptr].p1+1);
	  ca.entries[ptr+1].p1=parilong(buff);
	  ca.entries[ptr+1].p2=ca.entries[ptr].p2;
	  ca.entries[ptr+1].e=ca.entries[ptr].e;
	  ca.entries[ptr+1].log_p1=log(int_double(ca.entries[ptr+1].p1));
	  if(ca.entries[ptr+1].e==1)
	    eps_inc(ca.entries+ptr+1);
	  else
	    eps_max(ca.entries+ptr+1);
	  ca.entries[ptr].p2=ca.entries[ptr].p1;
	  eps_max(ca.entries+ptr);
	  add_exp(ca.entries+ptr,lhs,rhs,n);
	  ca.next_entry++;
	}      
    }
  if(rc)
    {
      printf("Successful completion....\n");
      check_Robin(lhs,rhs,1);
      for(uint64_t i=0;i<ca.next_entry-1;i++)
	printf("%lu %lu %lu\n",ca.entries[i].p1,ca.entries[i].p2,ca.entries[i].e);
    }

  return 0;
}

