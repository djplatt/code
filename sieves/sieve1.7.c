//
// Uses precomputed prime database
//

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"
#include "../code/includes/pi_x.h"

//#define DEBUG

#define num_bits (sizeof(ptype)<<3)
#define log_num_bits (6)
#define all_ones ((num_bits<<1)-1)

#define primep(p,v) (v[p>>(log_num_bits+1)]&mask1[p&all_ones])
#define primep1(p,q,v) (v[q]&mask1[p&all_ones])
//#define primep(p,v) (v[p/(num_bits*2)]&mask1[p&all_ones])
#define clear_prime(p,v) (v[p>>(log_num_bits+1)]&=mask[p&all_ones])
#define clear_prime1(p,q,v) (v[q]&=mask[p&all_ones])
//#define clear_prime(p,v) (v[p/(num_bits*2)]&=mask[p&all_ones])

#define FNAME "/home/dave/work/data/primes.dat"
#define BUFF_SIZE (8192)
#define V_BIG_PRIME ((ptype) 0xFFFFFFFFFFFFFFFEL)
#define ALL_ONES ((ptype) 0xFFFFFFFFFFFFFFFFL)

unsigned char buff[BUFF_SIZE];
int buff_ptr=BUFF_SIZE; // forces a read straight away

// returns x%y
// x has at most 80 significant bits (1.2e24)
// y has at most 48 (2.8e14)
inline __uint64_t rem_80_48(__uint128_t x, __uint64_t y)
{
  __uint64_t w=x&0xFFFF,z=x>>16;
  return((((z%y)<<16)+w)%y);
}

unsigned char next_char(FILE *infile)
{
  int dummy;
  if(buff_ptr==BUFF_SIZE)
    {
      buff_ptr=0;
      dummy=fread(buff,sizeof(unsigned char),BUFF_SIZE,infile);
    }
  return(buff[buff_ptr++]);
}

ptype last_prime=3; // first gap in file is from 3 to 5

ptype next_prime(FILE *infile)
{
  ptype gap=next_char(infile);
  if(gap==0)
    {
      gap=next_char(infile);
      if(gap==0) // end of file
	return(V_BIG_PRIME); // this is even and too big
      gap=(gap<<8)+next_char(infile);
      last_prime+=gap<<1;
      return(last_prime);
    }
  last_prime+=gap<<1;
  return(last_prime);
}

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}

ptype mask1[num_bits<<1],mask[num_bits<<1];

inline void set_masks ()
{
  ptype i,j;
  for(i=1,j=1;i<(num_bits<<1);i+=2,j<<=1)
    {
      mask1[i]=j;
      mask[i]=ALL_ONES;
      mask[i]^=j;
    }
}

inline void set_primes(ptype *vec, const ptype len)
{
  ptype i;
  for(i=0;i<=(len>>(log_num_bits+1));i++)
    vec[i]=ALL_ONES;
}

int rem_count=0;

inline ptype get_index(ptype p, bigint start)
{
  ptype rem=rem_80_48(start,p);//start%p;
  if(rem==0)
    return(1);
  rem=p-rem+1;
  if(rem&1)
    return(rem);
  else
    return(rem+p);
}

// for each prime p, find the first odd multiple of p within the target
// then cross it out and every 2p thereafter
inline void sieve (ptype *target, bigint target_start, const ptype target_len, FILE *infile)
{
  ptype i,p,q,ptr,cache_ptr,next_ptr,next_cache_ptr;
  ptype max_p=sqrt(target_start+target_len-1);
  ptype lim1=target_len>>1;
  
  if(lim1>max_p)
    fatal_error("must have L/2<=sqrt(end of sieve).");
  
  //printf("maximum prime to sieve=%lu\n",max_p);
  // p in [3,L/2)
  // we can get np,(n+2)p,(n+4)p... in the target range
  // or (n+1)p,(n+3)p... (if n is even)
  // so we must have >=1 to cross out
  for(p=3;p<=lim1;p=next_prime(infile))
    {
      //printf("p=%lu\n",p);
      ptr=get_index(p,target_start);
      cache_ptr=ptr>>(log_num_bits+1);
      __builtin_prefetch(&target[cache_ptr],1,0);
      q=p<<1;
      
#ifdef DEBUG
      if(ptr>target_len)
	{
	  printf("Prime missed target array altogether in first loop.\np=");
	  print_bigint(p);
	  printf("\n");
	  continue;
	}
#endif
      for(;;)
	{
	  next_ptr=ptr+q;
	  next_cache_ptr=next_ptr>>(log_num_bits+1);
	  __builtin_prefetch(&target[next_cache_ptr],1,0);
	  if(ptr>target_len)
	    break;
	  clear_prime1(ptr,cache_ptr,target);
	  ptr=next_ptr+q;
	  cache_ptr=ptr>>(log_num_bits+1);
	  __builtin_prefetch(&target[cache_ptr],1,0);
	  if(next_ptr>target_len)
	    break;
	  clear_prime1(next_ptr,next_cache_ptr,target);
	}
    }
#ifdef DEBUG
  printf("Done small primes < %lu\n",lim1);
#endif
  // since p > len/2, we cannot get np and (n+2)p in one target
  // so we must have <=1 to cross out

  // no point being cache clever here
  // only a small percentage are "hits" anyway
   for(;p<=max_p;p=next_prime(infile)) // <= one prime in target
    {
      ptr=get_index(p,target_start);
      if(ptr<=target_len)
	clear_prime(ptr,target);
#ifdef DEBUG
      if(ptr+(p<<1)<=target_len)
	{
	  printf("missed one or more prime multiples in end loop.\np=");
	  print_bigint(p);
	  printf("\n");
	}
#endif
    }
}


int main()
{
  __int64_t del,del_sum,del2;
  ptype *target,tl,tl2;
  ptype i,j,count,target_len=TARGET_LEN;
  bigint target_start,del2_sum;
  FILE *infile;
  ptype num_its=1;

  if(!(infile=fopen(FNAME,"rb")))
    fatal_error("Failed to open infile for binary input.");

  set_masks();

  for(i=0,target_start=1;i<22;i++,target_start*=10);
  target_start+=1; // must be odd
  target_start-=target_len*num_its;

  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for target sieve.");

  for(j=0;j<num_its;j++,target_start+=target_len)
    {  
      set_primes(target,target_len);
      sieve(target,target_start,target_len,infile);

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
	  printf("we have %lu primes\nsigma (p-t0)=%ld\nsigma (p-t0)^2=",count,del_sum);
	  print_bigint(del2_sum);
	  printf("\n");
	}
      
      
      rewind(infile);      // back to start of prime gaps
      buff_ptr=BUFF_SIZE;  // force an immediate read
      last_prime=3;        // first gap is 3 to 5
    }
  
  fclose(infile);
  return(0);
}
