#include "inttypes.h"
#include "stdio.h"
#include "stdlib.h"

#define LOG_HASH_LEN (30)
#define HASH_LEN (1LL<<LOG_HASH_LEN)
// hash is least sig bits of zero
#define HASH_MASK (HASH_LEN-1)

typedef struct
{
  uint64_t *z;
  uint64_t *N;
  uint64_t *c;
  uint8_t *t;
} hash_t;

// if top bit set is nth, then L function is degree n.
// degree 1
#define ZETA (2)
#define DIRICHLET (3)


// returns top 63 bits of 104 bit zero
// top bit resrved for Zeta zero at 14.
uint64_t read_zero(FILE *infile)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
  uint64_t res;

  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  // we want all 8 bits of c, all 32 bits of b and 23 bits of a 
  res=((uint64_t) c << 55)|((uint64_t) b << 23)|(a >> 41);
  uint64_t rest=a&((1LL<<42)-1);
  //printf("%LX %LX %LX %LX\n",res,a,rest,1LL<<40);
  if(rest>=(uint64_t)(1LL<<40))
    return res+1;
  else
    return res;
}

inline uint64_t hash_me(uint64_t zero)
{
  return zero&HASH_MASK;
}

void print_zero(uint64_t z)
{
  long double zz=z;
  zz/=(1LL<<60);
  printf("%l20.18e",zz);
}

void print_type(uint8_t t)
{
  switch(t)
    {
    case ZETA: printf("Zeta");break;
    case DIRICHLET: printf("Dirichlet");break;
    default: printf("unknown L-function type. Exiting.\n");exit(0);
    }
}

void print_z(uint64_t z, uint64_t N, uint64_t c, uint8_t t)
{
  print_zero(z);
  printf(" N:%lu char:%lu ",N,c);
  print_type(t);
  printf("\n");
}

bool insert_zero(uint64_t z, uint64_t N, uint64_t c, uint8_t t, hash_t *hash)
{
  uint64_t h=hash_me(z);
  while(hash->z[h])
    {
      if(z==hash->z[h])
	{
	  printf("duplicate found.\nwas adding ");
	  print_z(z,N,c,t);
	  printf("table contains ");
	  print_z(hash->z[h],hash->N[h],hash->c[h],hash->t[h]);
	  return false;
	}
      h++;
      h&=HASH_MASK;
    }
  hash->z[h]=z;
  hash->N[h]=N;
  hash->c[h]=c;
  hash->t[h]=t;
  return true;
}


