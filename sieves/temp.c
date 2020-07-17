#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"
#include "../includes/pi_x.h"
#include "../primegen/primegen.h"
#include "../primegen/primegen.c"
#include "../primegen/primegen_next.c"
#include "../primegen/primegen_init.c"


#define debug printf("Reached line number %d.\n",__LINE__)


#define num_bits (sizeof(ptype)<<3)
#define log_num_bits (6)
#define all_ones ((num_bits<<1)-1)

#define primep(p,v) (v[p>>(log_num_bits+1)]&mask1[p&all_ones])
//#define primep(p,v) (v[p/(num_bits*2)]&mask1[p&all_ones])
#define clear_prime(p,v) (v[p>>(log_num_bits+1)]&=mask[p&all_ones])
//#define clear_prime(p,v) (v[p/(num_bits*2)]&=mask[p&all_ones])

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}

void print_usage()
{
  printf("usage:- sieve1.9 <sieve_num> <num_its> <outfile>\n");
  printf("sieve_num 0 starts at -num_sieves\n");
  exit(0);

}

ptype mask1[num_bits<<1],mask[num_bits<<1];

inline void set_primes(ptype *vec, const ptype len)
{
  ptype i;
  //printf("in set_primes with len=%lu = %lu words.\n",len,len>>(log_num_bits+1));
  for(i=0;i<=(len>>(log_num_bits+1));i++)
    {
      vec[i]=0;
      vec[i]=~vec[i];
    }
}


// returns x%y
// x has at most 80 significant bits (1.2e24)
// y has at most 48 (2.8e14)
/*
inline __uint64_t rem_80_48(__uint128_t x, __uint64_t y)
{
  __uint64_t w=x&0xFFFF,z=x>>16;
  return((((z%y)<<16)+w)%y);
}
*/

inline __uint64_t rem_80_48(__uint128_t x, __uint64_t y)
{
  __uint64_t temp;
  /*
  __asm("movq %1,%%rax\n\t"
	"movq %2,%%rdx\n\t"
	"divq %3\n\t"
	"movq %%rax,%0\n\t"
	:"=m" (temp)
	:"m" (x), "m" (x), "m" (y)
	:"rax", "rdx");
  */
  __asm("movq $123,%%rax\n\t"
	"movq $1,%%rdx\n\t"
	"movq %1,%%rcx\n\t"
	"divq %%rcx\n\t"
	"movq %%rax,%0\n\t"
	:"=r" (temp)
	:"m" (y)
	:"rax","rdx","ecx");
  return(temp);
}

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



// primes to sieve are all odd
// assumes target length < 2 sqrt(target end)
inline void sieve1 (ptype *target, bigint target_start, const ptype target_len, primegen *pg)
{
  bigint target_end=target_start+target_len-1;
  ptype q,ptr,two_q,max_p=sqrt(target_end)/100.0;

  printf("max_p set to %ld\n",max_p);
  // now max_p is the largest integer s.t. max_p*max_p <= target end
  primegen_init(pg);
  primegen_next(pg); // 2
  for(q=primegen_next(pg),two_q=q<<1;two_q<=target_len;q=primegen_next(pg),two_q=q<<1)
    {
      //printf("prime=%lu\n",q);
#ifdef TEN_20
      if(q>max_p)
	return;
#endif
      ptr=get_index(q,target_start); // there must be at least one odd multiple in target
      for(;ptr<=target_len;ptr+=two_q)
	clear_prime(ptr,target);
    }
  for(;q<=max_p;q=primegen_next(pg))
    {
      ptr=get_index(q,target_start);
      if(ptr<=target_len) // might not even be an odd multiple in target
	clear_prime(ptr,target);
    }
}

int main(int argc, char **argv)
{
  ptype *target;
  long int sieve_num,num_its;
  ptype i,j,target_len=TARGET_LEN,it;
  bigint target_start,del2_sum,X0;
  ptype count,tl,tl2,segs_per_sieve=SEGS_PER_SIEVE,pix;
  __int64_t del,del2,del_sum;
  FILE *outfile;
  primegen pg[1];

  for(i=1,j=1;i<(num_bits<<1);i+=2,j<<=1)
    {
      mask1[i]=j;
      mask[i]=0;
      mask[i]=~mask[i];
      mask[i]^=j;
    }

  if(argc!=4)
    {
      printf("argc was %ld\n",argc);
      print_usage();
    }
  sieve_num=atoi(argv[1]);
  
  if(sieve_num<0)
    {
      printf("sieve_num was %ld\n",sieve_num);
      print_usage();
    }
  
  num_its=atoi(argv[2]);
  if(num_its<1)
    {
      printf("num_its was %ld\n",num_its);
      print_usage();
    }
  if(!(outfile=fopen(argv[3],"wb")))
    fatal_error("Failed to open outfile for binary write.\n");

  fwrite(&num_its,sizeof(long int),1,outfile);
  fwrite(&segs_per_sieve,sizeof(ptype),1,outfile);

  for(i=0,X0=1;i<LOG_10_X0;i++,X0*=10);

  print_bigint(X0);printf("\n");
  printf("%123=%lu\n",rem_80_48(X0,123));
  exit(0);

  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for target sieve.");

  target_start=X0-XI*SEGS_PER_SIEVE*NUM_SIEVES/2;
  target_start+=XI*SEGS_PER_SIEVE*sieve_num+1;
  for(it=sieve_num;it<num_its+sieve_num;it++,target_start+=target_len)
    {
      printf("Running sieve starting at ");print_bigint(target_start);printf("\n");
      fwrite(&it,sizeof(ptype),1,outfile);
      fwrite(&target_start,sizeof(bigint),1,outfile);

      set_primes(target,target_len);
      sieve1(target,target_start,target_len,pg);
      
      for(tl2=XI,tl=XI<<1,i=1;i<=target_len;tl2+=XI<<1,tl+=XI<<1)
	{
	  printf("from ");
	  print_bigint(target_start+i-1);
	  printf(" to ");
	  print_bigint(target_start+tl-1);
	  printf(" with t0=");
	  print_bigint(target_start+tl2-1);
	  printf("\n");
	  for(count=0,del_sum=0,del2_sum=0;i<=tl;i+=2)
	    if(primep(i,target))
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
