//
// Version Original
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "../primesieve/src/soe/PrimeSieve.h"

#define START 0x400000000000L // 2^46
#define LEN 0x200000L

#define debug printf("Reached line number %d.\n",__LINE__)

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}

void print_usage()
{
  printf("usage:- sieve <sieve_num> <num_its> <cache in Kbytes> <outfile>\n");
  printf("sieve_num 0 starts at %lu\n",START);
  exit(0);

}

PrimeSieve ps;
long unsigned int prime_count;

inline void sieve_zero(long unsigned int p)
{
  prime_count++;
}

inline void sieve (FILE *outfile, unsigned long int start)
{
  prime_count=0;
  ps.generatePrimes(start+1,start+LEN,sieve_zero); // sieve the "small" primes  
  fprintf(outfile,"Primes in %lu to %lu %lu\n",start+1,start+LEN,prime_count);
}

int main(int argc, char **argv)
{
  printf("Command line:- ");
  for(int i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=5)
    {
      printf("argc was %ld\n",argc);
      print_usage();
    }

  ps.setSieveSize(atoi(argv[3]));
  unsigned long int start=START+atol(argv[1])*LEN;
  unsigned long int num_its=atol(argv[2]);
  FILE *outfile=fopen(argv[4],"w");
  for(unsigned long int i=0;i<num_its;i++,start+=LEN)
      sieve(outfile,start);      
  return(0);
}
