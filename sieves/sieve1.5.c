#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"

#define debug printf("Reached line number %d.\n",__LINE__)

typedef long int ptype;
typedef __int128_t bigint;

#define num_bits (sizeof(ptype)<<3)
#define log_num_bits (6)
#define all_ones ((num_bits<<1)-1)

#define primep(p,v) (v[p>>(log_num_bits+1)]&mask1[p&all_ones])
//#define primep(p,v) (v[p/(num_bits*2)]&mask1[p&all_ones])
#define clear_prime(p,v) (v[p>>(log_num_bits+1)]&=mask[p&all_ones])
//#define clear_prime(p,v) (v[p/(num_bits*2)]&=mask[p&all_ones])

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
  ptype i;
  //printf("in set_primes with len=%lu = %lu words.\n",len,len>>(log_num_bits+1));
  for(i=0;i<=(len>>(log_num_bits+1));i++)
    {
      vec[i]=0;
      vec[i]=~vec[i];
    }
}

// create base primes 2,3,5,...
inline void erat_sieve (ptype *erat, const ptype len)
{
  ptype i,j,sqrt_len;
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


// sieving using base primes 2,3,5... held in small_primes
inline void sieve (ptype *target, bigint target_start, const ptype target_len, const ptype *small_primes, const ptype pix)
{
  ptype i,p,q,ptr;
  ptype max_p=sqrt(target_start+target_len-1);
  //printf("in sieve with max_p=%lu\n",max_p);
  for(i=1;i<pix;i++) // start at 3
      {
	p=small_primes[i]; // 3,5,7...
	
	if(p>max_p)
	  return;
	
	ptr=get_index(p,target_start);
	q=p<<1;
	for(;ptr<=target_len;ptr+=q) // skip by 2p each time
	  
	    clear_prime(ptr,target);
	  
      }
  //num_prims(target,target_start,target_len);
}

// primes to sieve are all odd
inline void sieve1 (ptype *target, bigint target_start, const ptype target_len, const ptype *source, const ptype source_start, const ptype source_len)
{
  ptype p,q,ptr;
  ptype source_offset=source_start-1,rem,two_q;
  ptype max_p=sqrt(target_start+target_len-1);
  for(p=1;p<=source_len;p+=2) // skip by 2 as all primes are odd
    if(primep(p,source))
      {
	q=p+source_offset; // q is the prime
	
	if(q>max_p)
	  return;
	
	ptr=get_index(q,target_start);
	// there is no guarentee that this ptr will fall within target
	if(ptr>target_len)
	  continue;
	two_q=(q<<1);
	for(;ptr<=target_len;ptr+=two_q) // skip by 2q
	  clear_prime(ptr,target);
      }
}
// primes to sieve are all odd
// used when there can be at most one value in target per prime
// i.e. when source_start>target_len
inline void fast_sieve (ptype *target, bigint target_start, const ptype target_len, const ptype *source, const ptype source_start, const ptype source_len)
{
  ptype p,q,ptr;
  ptype source_offset=source_start-1;
  ptype max_p=sqrt(target_start+target_len-1);
  for(p=1;p<=source_len;p+=2) // skip by 2 as all primes are odd
    if(primep(p,source))
      {
	q=p+source_offset; // q is the prime
	
	if(q>max_p)
	  return;
	
	ptr=get_index(q,target_start);
	// there is no guarentee that this ptr will fall within target
	if(ptr<=target_len)
	  clear_prime(ptr,target);
      }
}

int main()
{
  ptype *erat,*source,*target,*small_primes;
  ptype i,j,erat_len,source_len=4444440,target_len=16000000000;
  bigint target_start;
  ptype count,source_start,sqrt_end;


  for(i=0,target_start=1;i<22;i++,target_start*=10);
  target_start+=1;
  sqrt_end=sqrt(target_start+target_len-1);
  erat_len=sqrt(sqrt_end);
  erat_len=8000000;
  if(erat_len&1)
    erat_len++; // make sure it's even
  if(source_len&1)
    fatal_error("source length must be even.\n");

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
  sieve(target,target_start,target_len,small_primes,count);

  for(i=0,source_start=erat_len+1;source_start<sqrt_end;i++,source_start+=source_len)
    {
      if((i%100)==0)
	printf("sieving from %lu to %lu\n",source_start,source_start+source_len-1);
      // remove small primes from source
      set_primes(source,source_len);
      sieve(source,source_start,source_len,small_primes,count);
      // remove source primes from target
      if(source_start>target_len)
	fast_sieve(target,target_start,target_len,source,source_start,source_len);
      else
	sieve1(target,target_start,target_len,source,source_start,source_len);
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
