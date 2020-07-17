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
//#define primep(p,v) (v[p/(num_bits*2)]&mask1[p&all_ones])
#define clear_prime(p,v) (v[p>>(log_num_bits+1)]&=mask[p&all_ones])
//#define clear_prime(p,v) (v[p/(num_bits*2)]&=mask[p&all_ones])

#define FNAME "/home/dave/work/data/primes.dat"
#define BUFF_SIZE (8192)
#define V_BIG_PRIME ((ptype) 0xFFFFFFFFFFFFFFFEL)
#define ALL_ONES ((ptype) 0xFFFFFFFFFFFFFFFFL)

unsigned char buff[BUFF_SIZE];
int buff_ptr=BUFF_SIZE; // forces a read straight away

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

inline ptype get_index(ptype p, bigint start)
{
  ptype rem=start%p;
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
  ptype i,p,q,ptr;
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
      for(;ptr<=target_len;ptr+=q) // skip by 2p each time
	clear_prime(ptr,target);
    }
#ifdef DEBUG
  printf("Done small primes < %lu\n",lim1);
#endif
  // since p > len/2, we cannot get np and (n+2)p in one target
  // so we must have <=1 to cross out
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
  ptype *target,tl,tl2,k;
  ptype i,j,count,target_len=TARGET_LEN;
  bigint target_start,del_sum2;
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
      
      for(tl2=XI,tl=XI<<1,i=1;
	  i<=target_len;
	  tl2+=XI<<1,tl+=XI<<1)
	{
#ifdef DEBUG
	  printf("summing from %lu to %lu with centre %lu\n",i,tl,tl2);
#endif
	  for(count=0,del_sum=0,del_sum2=0;i<tl;i+=2)
	    if(primep(i,target))
	      {
		count++;
		del=i-tl2;
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
