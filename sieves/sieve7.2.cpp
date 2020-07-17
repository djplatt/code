// 7.0 modified to use primesieve in place of primegen
// 7.2 modified to use bit sets
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"
#include <bitset>
#include "../includes/pi_x.h"
#include "../primesieve/src/soe/PrimeSieve.h"

#define debug printf("Reached line number %d.\n",__LINE__)

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}

void print_usage()
{
  printf("usage:- sieve7.0 <sieve_num> <num_its> <outfile> <cache in Kbytes>\n");
  printf("sieve_num 0 starts at -num_sieves\n");
  exit(0);

}

std::bitset <TARGET_LEN/2> target;

// returns x%y
// x has at most 80 significant bits (1.2e24)
// y has at most 48 (2.8e14)
inline __uint64_t rem_80_48(__uint128_t x, __uint64_t y)
{
  
  __uint64_t w=x&0xFFFF,z=x>>16;
  return((((z%y)<<16)+w)%y);
  
  /*
  // inline assembler is actually slower!
  __uint64_t temp,*_x,x1,x2;
  _x=(__uint64_t *) &x;
  //printf("rem_80_40 called with x=");print_bigint(x);printf(" and y=%lu\n",y);
  x1=_x[0];
  x2=_x[1];
  //printf("x1=%lu x2=%lu y=%lu\n",x1,x2,y);
  __asm__("movq %2,%%rax\n\t"
	  "xor %%rdx,%%rdx\n\t"
	  "divq %3\n\t"
	  "movq %1,%%rax\n\t"
	  "divq %3\n\t"
	  "movq %%rdx,%0\n\t"
	  :"=r" (temp)
	  :"r" (x1), "r" (x2), "r" (y)
	  :"rax", "rdx");
  //printf("rem_80_40 returning %lu\n",temp);
  return(temp);
  */
}
/*
inline ptype get_index(ptype p, bigint start)
{
  ptype rem=rem_80_48(start,p);
  if(rem==0)
    return(1);
  rem=p-rem+1;
  if(rem&1)
    return(rem);
  else
    return(rem+p);
}
*/

inline ptype get_index(ptype p, bigint start)
{
  //printf("p=%lu start=");print_bigint(start);printf("\n");
  ptype rem=p-rem_80_48(start,p);
  if(rem==p)
    return(0);
  //print_bigint(start);printf(" mod %lu=%lu\n",p,rem);
  if(rem&1)
    return(rem+p);
  else
    return(rem);
  /*
  ptype rem_bar=(~rem)&1;
  //printf("rem_bar=%lu\n",rem_bar);
  int64_t rem64=-((int64_t) rem_bar);
  //printf("rem64=%ld\n",rem64);
  ptype del_p=p&rem64;
  return(rem+del_p); // if rem is even, add p to it
  */
}

bigint target_start;

//ptype index_sum=0;
inline void sieve_zero(ptype q) // q<=length/2 so there must be one to sieve
{                               // and there might be lots.

  //index_sum+=q;
  //return;
  //ptype two_q=q<<1;
  //printf("Clearing multiples of %lu\n",q);
  ptype ptr=get_index(q,target_start)>>1;
  //printf("ptr=%lu\n",ptr);
  while(ptr<TARGET_LEN/2)
    {
      target.reset(ptr);
      ptr+=q;
    }
  /*
  if(q==19)
    {
      for(int i=0;i<100;i++)
	if(target.test(i))
	  printf("1");
	else
	  printf("0");
      printf("\n");
      exit(0);
    }
  */
}

inline void sieve_one(ptype q) // q>length/2 so there can be at most one to sieve
{
  //index_sum+=q;
  //return;

  //printf("In sieve_one with q=%lu\n",q);
  ptype ptr=get_index(q,target_start);
  /*
  ptype p1=ptr>>(log_num_bits+1);
  __builtin_prefetch(&target[p1],1,0); // go fetch it as soon as poss
  ptype p2=mask[ptr&all_ones];
  if(ptr<=target_len)
    target[p1]&=p2;
  */
  if(ptr<TARGET_LEN)
    target.reset(ptr>>1);
}


PrimeSieve ps;
// primes to sieve are all odd
// assumes target length < 2 sqrt(target end)
inline void sieve1 ()
{
  printf("In sieve1 with start=");print_bigint(target_start);printf("\n");
  bigint target_end=target_start+TARGET_LEN-1;
  ptype max_p=sqrt(target_end);

  // we can't trust double prec square root at large height
  // so just use it as a guide
  printf("max_p was %ld\n",max_p);
  while((bigint) max_p*max_p<target_end)
    max_p++;
  while((bigint) max_p*max_p>target_end)
    max_p--;
  printf("max_p set to %ld\n",max_p);
  ps.generatePrimes(3,TARGET_LEN>>1,sieve_zero); // sieve the "small" primes  
  ps.generatePrimes((TARGET_LEN>>1)+1,max_p,sieve_one); // sieve the large primes
}

int main(int argc, char **argv)
{
  long int sieve_num,num_its;
  ptype i,j,it;
  bigint del2_sum,X0;
  ptype count,tl,tl2,segs_per_sieve=SEGS_PER_SIEVE,pix;
  __int64_t del,del2,del_sum;
  FILE *outfile;
  
  if(argc!=5)
    {
      printf("argc was %ld\n",argc);
      print_usage();
    }
  sieve_num=atol(argv[1]);
  
  if(sieve_num<0)
    {
      printf("sieve_num was %ld\n",sieve_num);
      print_usage();
    }
  
  num_its=atol(argv[2]);
  if(num_its<1)
    {
      printf("num_its was %ld\n",num_its);
      print_usage();
    }
  if(!(outfile=fopen(argv[3],"wb")))
    fatal_error("Failed to open outfile for binary write.\n");

  ps.setSieveSize(atoi(argv[4]));


  fwrite(&num_its,sizeof(long int),1,outfile);
  fwrite(&segs_per_sieve,sizeof(ptype),1,outfile);

  for(i=0,X0=1;i<LOG_10_X0;i++,X0*=10);

  target_start=X0-XI*SEGS_PER_SIEVE*NUM_SIEVES/2;
  target_start+=XI*SEGS_PER_SIEVE*sieve_num+1;
  for(it=sieve_num;it<num_its+sieve_num;it++,target_start+=TARGET_LEN)
    {
      printf("Running sieve starting at ");print_bigint(target_start);printf("\n");
      fwrite(&it,sizeof(ptype),1,outfile);
      fwrite(&target_start,sizeof(bigint),1,outfile);

      target.set();
      sieve1();

      //printf("index_sum=%lu\n",index_sum);
      //exit(0);

      for(tl2=XI,tl=XI<<1,i=1;i<=TARGET_LEN;tl2+=XI<<1,tl+=XI<<1)
	{
	  printf("from ");
	  print_bigint(target_start+i-1);
	  printf(" to ");
	  print_bigint(target_start+tl-1);
	  printf(" with t0=");
	  print_bigint(target_start+tl2-1);
	  printf("\n");
	  for(count=0,del_sum=0,del2_sum=0;i<=tl;i+=2)
	    if(target.test(i>>1))
	      {
		count++;
		del=i-tl2;
		del_sum+=del;
		del2=del*del;
		del2_sum+=del2;
		/*
		  printf("prime was ");print_bigint(target_start+i+1);
		  printf("\ndel was %ld\ndel^2 was %ld\n",del,del2);
		  printf("del_sum=%ld\ndel2_sum=",del_sum);
		  print_bigint(del2_sum);
		  printf("\n\n");
		  if(count==20)
		  return(0);
		*/
	      }

	  fwrite(&count,sizeof(ptype),1,outfile);
	  fwrite(&del_sum,sizeof(__int64_t),1,outfile);
	  fwrite(&del2_sum,sizeof(bigint),1,outfile);
	  printf("we have %lu primes\nsigma (p-t0)=%ld\nsigma (p-t0)^2=",count,del_sum);
	  print_bigint(del2_sum);
	  printf("\n");
	}
    }
  return(0);
}
