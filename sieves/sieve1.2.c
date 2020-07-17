#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"

#define debug printf("Reached line number %d.\n",__LINE__)

typedef long unsigned int ptype;
typedef __uint128_t bigint;

#define num_bits (sizeof(ptype)<<3)
#define log_num_bits (6)
#define all_ones ((num_bits<<1)-1)

#define primep(p,v) v[p>>(log_num_bits+1)]&mask1[p&all_ones]
#define clear_prime(p,v) v[p>>(log_num_bits+1)]&=mask[p&all_ones]

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
};


inline void print_bigint(bigint i)
{
  if(i<10)
    printf("%1lu",(long unsigned int) i);
  else
    {
      print_bigint(i/10);
      printf("%lu",(long unsigned int) (i%10));
    }
}

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}

ptype mask1[num_bits<<1],mask[num_bits<<1];

inline void set_primes(ptype *vec, const ptype len)
{
  long unsigned int i;
  //printf("in set_primes with len=%lu = %lu words.\n",len,len>>(log_num_bits+1));
  for(i=0;i<=(len>>(log_num_bits+1));i++)
    {
      vec[i]=0;
      vec[i]=~vec[i];
    }
}

// create base primes 2,3,5,...
inline void erat_sieve (ptype *erat, const unsigned long int len)
{
  long unsigned int i,j,sqrt_len;
  for(i=1,j=1;i<(num_bits<<1);i+=2,j<<=1)
    {
      mask1[i]=j;
      mask[i]=0;
      mask[i]=~mask[i];
      mask[i]^=j;
    }
  
  //printf("mask1[%lu]=%lX\nmask[%lu]=%lX\n",9,mask1[9],9,mask[9]);
  set_primes(erat,len);
  //printf("mask1[%lu]=%lX\nmask[%lu]=%lX\n",9,mask1[9],9,mask[9]);
  sqrt_len=sqrt(len);
  for(i=3;i<=sqrt_len;i+=2)
    {
      //printf("target[0]=%lX\n",erat[0]);
      if(primep(i,erat))
	{
	  //printf("prime %lu\n",i);
	  for(j=i*i;j<=len;j+=i+i)
	    {
	      //printf("   clearing %lu\n",j);
	      clear_prime(j,erat);
	    }
	}
    }
}

inline ptype get_index(ptype p, bigint start)
  ptype rem=start%p;
  if(rem==0)
    return(1);
  rem=p-rem+1;
  if(rem&1)
    return(rem);
  else
    return(rem+p);
}
*/

// sieving using base primes 2,3,5... held in small_primes
inline void sieve (ptype *target, bigint target_start, const ptype target_len, const ptype *small_primes, const ptype pix)
{
  ptype i,j,p,two_p,ptr;
  ptype thirty_p,start,start_rem,two_p_rem;
  ptype max_p=sqrt(target_start+target_len-1);
  for(i=1;i<pix;i++) // start at 3
      {
	p=small_primes[i]; // 7,11,13,...
	//printf("p=%lu\n",p);
	if(p>max_p)
	  return;
	two_p=p<<1;
	start=get_index(p,target_start);	
	thirty_p=p*30;
	if(primep(start,target))
	  for(ptr=start;ptr<=target_len;ptr+=thirty_p)
	    clear_prime(ptr,target);
	for(j=1;i<15;j++)
	  {
	    start+=two_p; // we don't store evens
	    if(start>target_len)
	      break;
	    if(primep(start,target))
	      for(ptr=start;ptr<=target_len;ptr+=thirty_p)
		clear_prime(ptr,target);
	  }
      }
}

// primes to sieve are all odd
inline void sieve1 (ptype *target, bigint target_start, const ptype target_len, const ptype *source, const ptype source_start, const ptype source_len)
{
  ptype p,q,i,thirty_q,start,ptr,count=0;
  ptype source_offset=source_start-1,two_q;
  ptype max_p=sqrt(target_start+target_len-1);
  ptype start_rem,two_q_rem;
  for(p=1;p<=source_len;p+=2) // skip by 2 as all primes are odd
    if(primep(p,source))
      {
	q=p+source_offset; // q is the prime
	if(q>max_p)
	  return;
	start=get_index(q,target_start);
	// there is no guarentee that this ptr will fall within target
	if(start>target_len)
	  continue;
	thirty_q=q*30;
	if(primep(start,target))
	  for(ptr=start;ptr<=target_len;ptr+=thirty_q)
	    clear_prime(ptr,target);
	two_q=(q<<1);
	for(i=1;i<15;i++)
	  {
	    start+=two_q;
	    if(start>target_len)
	      break;
	    if(primep(start,target))
	      for(ptr=start;ptr<=target_len;ptr+=thirty_q)
		clear_prime(ptr,target);
	  }
      }
}

int main()
{
  ptype *erat,*source,*target,*small_primes;
  ptype i,j,erat_len,source_len=100000,target_len=400000;
  bigint target_start;
  ptype count,source_start,sqrt_end;


  for(i=0,target_start=1;i<12;i++,target_start*=10);
  target_start+=1;
  sqrt_end=sqrt(target_start+target_len-1);
  erat_len=sqrt(sqrt_end);
  if(erat_len&1)
    erat_len++;
  if(source_len&1)
    fatal_error("source length must be even.\n");
  if(target_len<=erat_len)
    fatal_error("Target length must be large enough to ensure it contains at least one small prime multiple.\n");

  if(!(erat=(ptype *) malloc(sizeof(ptype)*(1+(erat_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for erat sieve.");
  if(!(source=(ptype *) malloc(sizeof(ptype)*(1+(source_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for source sieve.");
  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for target sieve.");
  erat_sieve(erat,erat_len);
  for(i=3,count=1;i<=erat_len;i+=2)
      if(primep(i,erat))
	{
	  //printf("%lu is prime\n",i);
	  count++;
	}
  printf("Pi(%lu)=%lu\n",erat_len,count);
  if(!(small_primes=(ptype *) malloc(sizeof(ptype)*count)))
    fatal_error("Failed to alloacte memory for small_primes.");
  small_primes[0]=2;
  for(i=3,j=1;j<count;i+=2)
    if(primep(i,erat))
      small_primes[j++]=i;

  // remove the small primes from the target
  set_primes(target,target_len);
  debug;  
  sieve(target,target_start,target_len,small_primes,count);
  debug;
  for(source_start=erat_len+1;source_start<sqrt_end;source_start+=source_len)
    {
      
      // remove small primes from source
      set_primes(source,source_len);
      debug;
      sieve(source,source_start,source_len,small_primes,count);
      // remove source primes from target
      printf("sieving %lu to %lu\n",(ptype) source_start,(ptype) source_start+source_len-1);
      sieve1(target,target_start,target_len,source,source_start,source_len);
      debug;
    }
  for(i=1,count=0;i<=target_len;i+=2)
    if(primep(i,target))
      {
	count++;
	print_bigint(target_start+i-1);
	printf(" is prime\n");
	if(count==40)
	  return(0);
      }
  printf("Pi(");print_bigint(target_start+target_len-1);printf(")-Pi(");
  print_bigint(target_start);printf(")=%lu\n",count);
  
  return(0);
}
