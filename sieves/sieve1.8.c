//
// Uses precomputed prime database
// Works mod 2,3,5

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"

//#define DEBUG
#define TIME

#ifdef TIME
#include "time.h"
#endif


#define debug printf("Reached line number %d.\n",__LINE__)

typedef __uint64_t ptype;
typedef __uint128_t bigint;

// can fit 240 flags into one 64 bit word

ptype mask1[240],mask[240];


inline ptype primep(ptype p, ptype *v) 
{
  ptype q,r; 
  q=p/240; r=p%240; // I hope these get done atomically!
  /*
  __asm("xorq %%rdx,%%rdx\n\t"
	"movq %2,%%rax\n\t"
	"divq %3\n\t"
	"movq %%rax,%0\n\t"
	"movq %%rdx,%1\n\t"
	: "=r" (q), "=r" (r)
	: "r" (p), "r" (t)
	: "rdx","rax");
  */
  return(v[q]&mask1[r]);
}

inline void clear_prime(ptype p, ptype *v) 
{
  ptype q,r; 
  q=p/240; r=p%240;
  /*
  __asm("xorq %%rdx,%%rdx\n\t"
	"movq %2,%%rax\n\t"
	"divq %3\n\t"
	"movq %%rax,%0\n\t"
	"movq %%rdx,%1\n\t"
	: "=r" (q), "=r" (r)
	: "r" (p), "r" (t)
	: "rdx","rax");
  */
  v[q]&=mask[r];
}

#define clear_prime1(q,r,v) (v[q]&=mask[r])



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


/*
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
*/

inline ptype gcd (ptype a, ptype b)
// Euclid algorithm gcd
// best if a<=b ?
{
  unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
}


inline void set_masks()
{
  ptype i,j;
  for(i=0,j=1;i<240;i++)
    {
      if(gcd(30,i)==1)
	{
	  mask1[i]=j;
	  mask[i]=ALL_ONES;
	  mask[i]^=j;
	  j<<=1;
	}
      else
	{
	  mask1[i]=0;
	  mask[i]=ALL_ONES;
	}
    }
}

inline void set_primes(ptype *vec, const ptype len)
{
  ptype i;
  for(i=0;i<=(len/240);i++)
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

time_t this_t,last_t;

// for each prime p, find the first odd multiple of p within the target
// then cross it out and every 2p thereafter
inline void sieve (ptype *target, bigint target_start, const ptype target_len, FILE *infile)
{
  ptype i,p,q,ptr,cache_ptr,next_ptr,next_cache_ptr;
  ptype max_p=sqrt(target_start+target_len-1);
  ptype lim1=target_len>>1,cache_rem,next_cache_rem;
  
  if(lim1>max_p)
    fatal_error("must have L/2<=sqrt(end of sieve).");
  
  //printf("maximum prime to sieve=%lu\n",max_p);
  // p in [3,L/2)
  // we can get np,(n+2)p,(n+4)p... in the target range
  // or (n+1)p,(n+3)p... (if n is even)
  // so we must have >=1 to cross out
  p=next_prime(infile);
  p=next_prime(infile);
  for(;p<=lim1;p=next_prime(infile))
    {
      //printf("p=%lu\n",p);
      ptr=get_index(p,target_start);
      cache_ptr=ptr/240;
      cache_rem=ptr%240;
      /*
  __asm("xorq %%rdx,%%rdx\n\t"
	"movq %2,%%rax\n\t"
	"divq %3\n\t"
	"movq %%rax,%0\n\t"
	"movq %%rdx,%1\n\t"
	: "=r" (cache_ptr), "=r" (cache_rem)
	: "r" (ptr), "r" (t)
	: "rdx","rax");
      */
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
	  next_cache_ptr=next_ptr/240;
	  next_cache_rem=next_ptr%240;
	  /*
  __asm("xorq %%rdx,%%rdx\n\t"
	"movq %2,%%rax\n\t"
	"divq %3\n\t"
	"movq %%rax,%0\n\t"
	"movq %%rdx,%1\n\t"
	: "=r" (next_cache_ptr), "=r" (next_cache_rem)
	: "r" (next_ptr), "r" (t)
	: "rdx","rax");
	  */
	  __builtin_prefetch(&target[next_cache_ptr],1,0);
	  if(ptr>target_len)
	    break;
	  clear_prime1(cache_ptr,cache_rem,target);
	  ptr=next_ptr+q;
	  cache_ptr=ptr/240;
	  cache_rem=ptr%240;
	  /*
  __asm("xorq %%rdx,%%rdx\n\t"
	"movq %2,%%rax\n\t"
	"divq %3\n\t"
	"movq %%rax,%0\n\t"
	"movq %%rdx,%1\n\t"
	: "=r" (cache_ptr), "=r" (cache_rem)
	: "r" (ptr), "r" (t)
	: "rdx","rax");
	  */

	  __builtin_prefetch(&target[cache_ptr],1,0);
	  if(next_ptr>target_len)
	    break;
	  clear_prime1(next_cache_ptr,next_cache_rem,target);
	}
    }
#ifdef TIME
  printf("Done small primes < %lu\n",lim1);
  this_t=time(NULL);
  printf("time elapsed=%f\n",difftime(this_t,last_t));
  last_t=this_t;
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
#ifdef TIME
  this_t=time(NULL);
  printf("finished sieving.\n");
  printf("time elapsed=%f\n",difftime(this_t,last_t));
  last_t=this_t;
#endif

}

#define TARGET_POW_2 (36) // will be 30*2^(TP-5) long

#define TARGET_LEN ((((long unsigned int) 1) << (TARGET_POW_2-5))*30)

// target will get split into 2^TS segments, each of length <2^32 bits
// if TP2<=33, TS=0
#if (TARGET_POW_2>33)
#define TARGET_SHIFT (TARGET_POW_2-1-32) 
#else
#define TARGET_SHIFT (0)
#endif

int main()
{
  long int del,del_sum,del2;
  ptype *target,tl,tl2,k;
  ptype i,j,count,target_len=TARGET_LEN;
  bigint target_start,del_sum2;
  FILE *infile;
  ptype num_its=1;
  if(!(infile=fopen(FNAME,"rb")))
    fatal_error("Failed to open infile for binary input.");

  set_masks();

  for(i=0,target_start=1;i<22;i++,target_start*=10);
  target_start-=target_len*num_its;

  while((target_start%30)!=1)
    target_start++;
  printf("target_len=%lu\n",target_len);
  printf("Going to allocate %lu bytes.\n",(1+(target_len/240)*8));
  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len/240)))))
    fatal_error("Failed to allocate memory for target sieve.");

  last_t=time(NULL);
  for(j=0;j<num_its;j++,target_start+=target_len)
    {  
      set_primes(target,target_len);
      sieve(target,target_start,target_len,infile);
#ifdef DEBUG
      for(i=1,count=0;count<20;i++)
	if(primep(i,target))
	  {
	    count++;
	    print_bigint(i+target_start-1);
	    printf(" is prime\n");
	  }
#endif
      
      for(k=0,tl=target_len>>TARGET_SHIFT,tl2=target_len>>(TARGET_SHIFT+2),i=1;
	  k<(1<<TARGET_SHIFT);
	  k++,tl+=target_len>>TARGET_SHIFT,tl2+=target_len>>(TARGET_SHIFT+1))
	{
	  printf("summing from %lu to %lu with centre %lu\n",i,tl,tl2<<1);
	  for(count=0,del_sum=0,del_sum2=0;i<tl;i+=2)
	    if(primep(i,target))
	      {
		count++;
		del=(i>>1)-tl2; // i not even, 1->0, 3->1 etc.
		del_sum+=del;
		del2=del*del; // this will be 31*31 signed -> 64 signed
		del_sum2+=del2; // 128 unsigned + 64 signed
	      }
	  printf("sigma 1 = %lu\n",count);
	  printf("sigma i-B/2 = %ld\n",del_sum);
	  printf("sigma (i-B/2)^2 =");
	  print_bigint(del_sum2);
	  printf("\n");
	}
      
      rewind(infile);      // back to start of prime gaps
      buff_ptr=BUFF_SIZE;  // force an immediate read
      last_prime=3;        // first gap is 3 to 5
    }
  
  fclose(infile);
  return(0);
}
