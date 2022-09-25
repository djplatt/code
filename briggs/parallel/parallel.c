//#include <pari/pari.h>
#include "primesieve.h"
#include "arb.h"
#include "inttypes.h"

int main(int argc, char** argv)
{
  printf("Command line:- ");
  int i;
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  
  if(argc!=4)
    {
      printf("Usage:- %s <prec> <start> <end>\n",argv[0]);
      return 0;
    }
  
  //pari_init(500000,0);
  uint64_t prec=atol(argv[1]);
  uint64_t start=atol(argv[2]);
  uint64_t end=atol(argv[3]);

  // setup primesieve
  primesieve_iterator it;
  primesieve_init(&it);
  primesieve_skipto(&it,start,end);

  uint64_t p0,p,p1;
  arb_t n,sig;
  arb_init(n);arb_init(sig);
  arb_set_ui(sig,1);arb_set_ui(n,1);
  p0=primesieve_next_prime(&it); // first prime used
  p=p0;
  while(p <= end)
    {
      p1=p; // last prime used
      arb_mul_ui(n,n,p,prec); // n <- prod p in [p0,p1]
      arb_mul_ui(sig,sig,p+1,prec); // sig <- sigma(n)
      p=primesieve_next_prime(&it);
    }
  
  printf("%lu %lu\n",p0,p1); // <first prime> <last prime>
  arb_dump_file(stdout,n);
  printf("\n");
  arb_dump_file(stdout,sig);
  printf("\n");
  //arb_printd(n,40);printf("\n");
  //arb_printd(sig,40);printf("\n");
  arb_clear(n);
  arb_clear(sig);
  //pari_close();
  return 0;
}
