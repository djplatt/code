// check_proth.c
//
#include "stdio.h"
#include "inttypes.h"
#include "stdlib.h"
#include "gmp.h"


#define STEP_SIZE (4000000000000000000LL)

mpz_t pp,pe,pq,a;
void init_proth()
{
  mpz_init(pp);mpz_init(pe);mpz_init(pq);mpz_init(a);
}

int proth_p (uint64_t h, uint64_t n, uint64_t aa)
{
  mpz_set_ui(pe,h);
  mpz_mul_2exp(pe,pe,n-1);
  mpz_add(pp,pe,pe);
  mpz_add_ui(pp,pp,1);
  mpz_set_ui(a,aa);
  mpz_powm(pq,a,pe,pp);
  mpz_add_ui(pq,pq,1);
  return(mpz_cmp(pq,pp)==0);
}

int main(int argc, char **argv)
{

  if(argc!=2)
    exit(0);

  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    {
      printf("Failed to open %s for binary input. Exiting.\n");
      exit(0);
    }

  init_proth();

  uint64_t n,h0,h1,a,h;
  fread(&n,sizeof(uint64_t),1,infile);
  fread(&h0,sizeof(uint64_t),1,infile);
  fread(&h1,sizeof(uint64_t),1,infile);

  fread(&a,sizeof(uint64_t),1,infile);
  fread(&h,sizeof(uint64_t),1,infile);

  uint64_t max_delta=STEP_SIZE>>n;
  uint64_t first_delta=max_delta>>1;
  uint64_t last_h=h;
  if(h>first_delta)
    {
      printf("First h chosen=%lu was too far from starting h=%lu,\n",h+h0,h0);
      exit(0);
    }
  if(a!=0) // not a PARI prime
    if(!proth_p(h+h0,n,a))
      {
	printf("Proth number at %lu*2^%lu+1 with a=%lu failed.\n",h+h0,n,a);
	exit(0);
      }

       while(fread(&a,sizeof(uint64_t),1,infile)==1)
	 {
	   if(fread(&h,sizeof(uint64_t),1,infile)!=1)
	     {
	       printf("Error reading infile.\n");
	       exit(0);
	     }
	   if(h-last_h>max_delta)
	     {
	       printf("Gap from h=%lu to h=%lu too large for n=%lu.\n",last_h+h0,h+h0,n);
	       exit(0);
	     }
	   last_h=h;
	   if(a!=0) // not a PARI prime
	     if(!proth_p(h+h0,n,a))
	       {
		 printf("Proth number at %lu*2^%lu+1 with a=%lu failed.\n",h+h0,n,a);
		 exit(0);
	       }
	 }

       if(h1-first_delta>last_h+h0)
	 {
	   printf("Final prime at or after h=%lu was too far from end at h=%lu for n=%lu.\n",last_h+h0,h1,n);
	   exit(0);
	 }
       fclose(infile);
       return(0);
}
