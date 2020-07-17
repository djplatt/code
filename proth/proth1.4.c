//***********************************************************
//
// proth1.4.c
//
//
// HA Helfgott & DJ Platt 2013
// 
// Search for a ladder of Proth primes
//

#include "inttypes.h"
#include "gmp.h"
#include "pari.h"

// This is how far Binary Goldbach has been checked
#define STEP_SIZE (4000000000000000000L)

uint64_t PRIME_LIMIT; // we sieve by primes less than this

#define SMALL_PRIME_LIMIT (10) // small primes used to test Proth numbers

uint64_t ismall_primes[SMALL_PRIME_LIMIT]={2,3,5,7,11,13,17,19,23,29};

mpz_t one,small_primes[SMALL_PRIME_LIMIT],pp,pe,pq;

uint64_t test_masks[64],clear_masks[64];

void clear_bit(uint64_t pos, uint64_t *bits)
{
  bits[pos>>6]&=clear_masks[pos&0x3F];
}

uint64_t test_bit(uint64_t pos, uint64_t *bits)
{
  return(bits[pos>>6]&test_masks[pos&0x3F]);
}

void init_proth()
{
  pari_init(1<<20,1<<20);
  mpz_init(one);
  mpz_set_ui(one,1);
  test_masks[0]=1;
  uint64_t i;
  for(i=1;i<64;i++)
    test_masks[i]=test_masks[i-1]<<1;
  for(i=0;i<64;i++)
    clear_masks[i]=~test_masks[i];
  for(i=0;i<SMALL_PRIME_LIMIT;i++)
    {
      mpz_init(small_primes[i]);
      mpz_set_ui(small_primes[i],ismall_primes[i]);
    }
  mpz_init(pp);mpz_init(pe);mpz_init(pq);
}


// find smallest k>=0 such that r1+k*r2=0 (mod p)
// Should use extended Euclid but its only called once
uint64_t solve_mod(uint64_t r1, uint64_t r2, uint64_t p)
{
  uint64_t k,r;
  for(k=0,r=r1;k<p;k++,r+=r2)
    if(r%p==0)
      return(k);

  fprintf(stderr,"solve_mod broke. Exiting.\n");
  exit(0);
}

// is h*2^n+1 a Proth prime
int proth_p (uint64_t h, uint64_t n, FILE *outfile)
{
  mpz_set_ui(pe,h);
  mpz_mul_2exp(pe,pe,n-1);
  mpz_add(pp,pe,pe);
  mpz_add_ui(pp,pp,1);
  int ptr;
  for(ptr=0;ptr<SMALL_PRIME_LIMIT;ptr++)
    if(mpz_jacobi(small_primes[ptr],pp)==-1)
      break;
  if(ptr==SMALL_PRIME_LIMIT) // non of the small primes worked
    return(0);
  mpz_powm(pq,small_primes[ptr],pe,pp);
  mpz_add_ui(pq,pq,1);
  int res=mpz_cmp(pq,pp);
  if(res==0) // its a prime so save it for double checking
    {
      fwrite(ismall_primes+ptr,sizeof(uint64_t),1,outfile);
      fwrite(&h,sizeof(uint64_t),1,outfile);
      return(1);
    }
  else
    return(0);
}

int main(int argc, char **argv)
{

  if(argc!=6)
    {
      printf("Usage:- %s <n> <h0> <h1> <sieve prime limit> < outfile>.\n",argv[0]);
      exit(0);
    }

  uint64_t zero=0;
  uint64_t n=atol(argv[1]);
  uint64_t h0=atol(argv[2]);
  uint64_t h1=atol(argv[3]);
  PRIME_LIMIT=atol(argv[4]);
  FILE *outfile=fopen(argv[5],"wb");
  if(!outfile)
    {
      fprintf(stderr,"Error opending %s for binary output. Exiting.\n",argv[5]);
      exit(0);
    }

  fwrite(&n,sizeof(uint64_t),1,outfile);
  fwrite(&h0,sizeof(uint64_t),1,outfile);
  fwrite(&h1,sizeof(uint64_t),1,outfile);

  mpz_t h0z,h1z;
  mpz_init(h0z);mpz_init(h1z);
  mpz_set_ui(h0z,h0);
  mpz_set_ui(h1z,h1);

  init_proth();

  mpz_t two_n;
  mpz_init(two_n);
  mpz_mul_2exp(two_n,one,n);
  if(mpz_cmp(h1z,two_n)>=0)
    {
      fprintf(stderr,"h1 exceeds 2^n. Exiting.\n");
      exit(0);
    }
 
  printf("n=%lu, 2^n=",n);mpz_out_str(NULL,10,two_n);printf("\n");

  printf("h0=");mpz_out_str(NULL,10,h0z);printf("\n");
  printf("h1=");mpz_out_str(NULL,10,h1z);printf("\n");

  mpz_t start,finish;
  mpz_init(start);mpz_init(finish);
  mpz_mul_2exp(start,h0z,n);
  mpz_add_ui(start,start,1);
  mpz_mul_2exp(finish,h1z,n);
  mpz_add_ui(finish,finish,1);

  printf("start=");mpz_out_str(NULL,10,start);printf("\n");
  printf("finish=");mpz_out_str(NULL,10,finish);printf("\n");


  uint64_t bit_array_len=((h1-h0)>>6)+1;

  uint64_t *bits=(uint64_t *)malloc(bit_array_len<<3);
  if(!bits)
    {
      fprintf(stderr,"Failed to allocate memory for bit array. Exiting.\n");
      exit(0);
    }

  uint64_t i;
  for(i=0;i<bit_array_len;i++)
    bits[i]=0xFFFFFFFFFFFFFFFFL;

  uint64_t pp;
  printf("enumerating primes to %ld\n",PRIME_LIMIT);
  uint8_t *primes=(uint8_t *) malloc(sizeof(uint8_t)*(PRIME_LIMIT+1));
  if(!primes)
    {
      fprintf(stderr,"Error allocating memory for primes. Exiting.\n");
      exit(0);
    }
  for(i=3;i<=PRIME_LIMIT;i+=2)
    primes[i]=1;
  for(i=4;i<=PRIME_LIMIT;i+=2)
    primes[i]=0;
  uint64_t ptr,ptr1;
  for(ptr=3;ptr<floor(sqrt(PRIME_LIMIT));ptr+=2)
    if(primes[ptr])
      for(ptr1=ptr*ptr;ptr1<=PRIME_LIMIT;ptr1+=ptr)
	primes[ptr1]=0;

  for(pp=3;pp<PRIME_LIMIT;)
    {
      // this could be done in integers
      uint64_t r1=mpz_fdiv_ui(start,pp);
      uint64_t r2=mpz_fdiv_ui(two_n,pp);
      uint64_t k=solve_mod(r1,r2,pp);
      uint64_t ptr=k;
      while(ptr<=(h1-h0))
	{
	  //printf("Clearing ptr %lu\n",ptr);
	  clear_bit(ptr,bits);
	  ptr+=pp;
	}
      pp+=2;
      while((pp<=PRIME_LIMIT)&&(!primes[pp]))
	pp+=2;
    }
  free(primes);

  printf("Searching for Proth primes...\n");
  ptr=0;
  uint64_t last_h=0,this_h;
  uint64_t first_h=STEP_SIZE>>(n+1);
  uint64_t del_h=STEP_SIZE>>n;
  GEN p;
  pari_sp sp;
  for(this_h=first_h;this_h<=h1-h0;)
    {
      while((!test_bit(this_h,bits))&&(this_h>last_h))
	this_h--;
      if(this_h==last_h)
	{
	  //printf("Failed to find proth prime after h=%lu.\n",last_h);
	  //printf("Using PARI library to find a general prime instead.\n");
	  last_h+=del_h;
	  sp=avma; // top of PARI stack
	  p=stoi(last_h+h0);
	  p=mpshift(p,n);
	  //p=addis(p,1);
	  p=precprime(p);
	  //output(p);
	  if(!isprime(p))
	    {
	      fprintf(stderr,"Pari produced a pseudoprime!.\n");
	      exit(0);
	    }
	  
	  else
	    {
	      pari_printf("%Pd Pari\n",p);
	    }
	  
	  //p=subis(p,1);
	  p=mpshift(p,-n);
	  last_h=itou(p);
	  fwrite(&zero,sizeof(uint64_t),1,outfile);
	  fwrite(&last_h,sizeof(uint64_t),1,outfile);
	  last_h-=h0;
	  avma=sp; // return Pari stack to virgin
	  this_h=last_h+del_h;
	  //printf("Resuming from h=%lu\n",this_h);
	  continue;
	}
      //printf("Proth testing for h=%lu\n",this_h+h0);
      if(!proth_p(this_h+h0,n,outfile))
	this_h--; // try the previous potential Proth prime
      else
	{
	  //printf("Proth prime found at %lu*2^%lu+1=",this_h+h0,n);
	  //printf("%lu Proth\n",this_h+h0);
	  //printf("\nChange in h=%lu\n",this_h-last_h);
	  last_h=this_h;
	  this_h=this_h+del_h;
	}
    }

  // the last rung has to be within STEP_SIZE/2 of the top
  if(last_h+del_h/2<h1-h0) // not close enough to top
    {
      sp=avma; // top of PARI stack
      p=stoi(h1);
      p=mpshift(p,n); // p=h1*2^n
      //p=addis(p,1);
      p=precprime(p); // p<h1*2^n
      //output(p);
      if(!isprime(p))
	{
	  fprintf(stderr,"Pari produced a pseudoprime!.\n");
	  exit(0);
	}
      
      else
	pari_printf("%Pd Pari\n",p);
      
      //p=subis(p,1);
      p=mpshift(p,-n); // pretend this prime was at a Proth number
      last_h=itou(p);
      fwrite(&zero,sizeof(uint64_t),1,outfile); // 0 means not Proth
      fwrite(&last_h,sizeof(uint64_t),1,outfile); // last_h*2^n+1 <= p 
      last_h-=h0;
      avma=sp; // return Pari stack to virgin
      if(last_h+del_h/2<h1-h0) // still not near enough to top
	fprintf(stderr,"Failed to find a prime near enough to top of range.\n");
    }
  fclose(outfile);
  return (0);
}
