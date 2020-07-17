#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"

#define debug printf("Reached line number %d.\n",__LINE__)

typedef long unsigned int ptype;
typedef __uint128_t bigint;

#define num_bits (sizeof(ptype)<<3)
#define log_num_bits (6)
#define all_ones (num_bits-1)

#define primep(p,v) v[p>>log_num_bits]&mask1[p&all_ones]
#define clear_prime(p,v) v[p>>log_num_bits]&=mask[p&all_ones]

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

ptype mask1[num_bits],mask[num_bits];

void num_prims(const ptype *source, const long unsigned int source_start, const long unsigned int source_len)
{
  long unsigned int i,count;
  for(i=1,count=0;i<=source_len;i++)
    if(primep(i,source))
      count++;
  printf("Pi(%lu)-Pi(%lu)=%lu\n",source_start+source_len-1,source_start,count);
}


inline void set_primes(ptype *vec, const unsigned long int len)
{
  long unsigned int i;
  for(i=0;i<=(len>>log_num_bits);i++)
    {
      vec[i]=0;
      vec[i]=~vec[i];
    }
}

// create base primes 2,3,5,...
inline void erat_sieve (ptype *erat, const unsigned long int len)
{
  long unsigned int i,j,sqrt_len;
  for(i=0,j=1;i<num_bits;i++,j<<=1)
    {
      mask1[i]=j;
      mask[i]=0;
      mask[i]=~mask[i];
      mask[i]^=j;
      //printf("mask1[%lu]=%lu\nmask[%lu]=%lu\n",i,mask1[i],i,mask[i]);
    }
  set_primes(erat,len);
  clear_prime(1,erat);
  sqrt_len=sqrt(len);
  for(i=4;i<=len;i+=2)
    clear_prime(i,erat);
  for(i=3;i<=sqrt_len;i+=2)
    if(primep(i,erat))
      for(j=i*i;j<=len;j+=i+i)
	clear_prime(j,erat);
}

// sieving using base primes 2,3,5... held in small_primes
inline void sieve (ptype *target, bigint target_start, const ptype target_len, const ptype *small_primes, const ptype pix)
{
  ptype i,p,q,ptr;
  ptype rem;
  ptype max_p=sqrt(target_start+target_len-1);
  
  // cross out evens
  for(ptr=1+(target_start&1);ptr<=target_len;ptr+=2)
    clear_prime(ptr,target);
  for(i=1;i<pix;i++) // start at 3
      {
	p=small_primes[i]; // 3,5,7...
	if(p>max_p)
	  break;
	q=p<<1;
	rem=target_start%p;
	if(rem==0) rem=p;
	ptr=p-rem+1;
	if(!((ptr+target_start-1)&1)) // first n*prime is even
	  {
	    clear_prime(ptr,target); // clear it
	    ptr+=p;                  // point to (n+1)*prime
	  }
	for(;ptr<=target_len;ptr+=q) // skip by 2p each time
	  clear_prime(ptr,target);
      }
  //num_prims(target,target_start,target_len);
}

// primes to sieve are all odd
inline void sieve1 (ptype *target, bigint target_start, const long unsigned int target_len, const ptype *source, const long unsigned int source_start, const unsigned long int source_len)
{
  long unsigned int p,q,ptr,source_end=source_start+source_len-1;
  long unsigned int source_offset=source_start-1,rem,two_q;
  long unsigned int max_p=sqrt(source_end)-source_offset;
  if(source_len<max_p)
    max_p=source_len;
  //printf("in sieve with start=%lu len=%lu\n",source_start,source_len);
  // don't need to go beyond p>sqrt(source_end)-source_offset
  //find first prime
  for(p=1;!primep(p,source);p++);
  for(;p<=max_p;p+=2) // skip by 2 as all primes are odd ?
    if(primep(p,source))
      {
	q=p+source_offset; // q is the prime
	two_q=(q<<1);
	rem=target_start%q;
	if(rem==0) rem=q;
	ptr=q-rem+1;
	// there is no guarentee that this ptr will fall within target
	if(ptr>target_len)
	  continue;
	if(!((ptr+target_start-1)&1)) // its 2m*p
	  {
	    clear_prime(ptr,target);
	    ptr+=q;                   // 
	    // ditto
	    if(ptr>target_len)
	      continue;
	  }
	for(;ptr<=target_len;ptr+=two_q) // skip by 2q
	  clear_prime(ptr,target);
      }
}

int main()
{
  ptype *erat,*source,*target,*small_primes;
  ptype i,j,erat_len,source_len=1000000,target_len=4000000000;
  bigint target_start=1;
  ptype count=0,source_start,sqrt_end;


  for(i=0,target_start=1;i<20;i++,target_start*=10);

  sqrt_end=sqrt(target_start+target_len-1);
  erat_len=sqrt(sqrt_end);

  if(!(erat=(ptype *) malloc(sizeof(ptype)*(1+(erat_len>>log_num_bits)))))
    fatal_error("Failed to allocate memory for erat sieve.");
  if(!(source=(ptype *) malloc(sizeof(ptype)*(1+(source_len>>log_num_bits)))))
    fatal_error("Failed to allocate memory for source sieve.");
  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len>>log_num_bits)))))
    fatal_error("Failed to allocate memory for target sieve.");
  erat_sieve(erat,erat_len);
  set_primes(target,target_len);
  for(i=2;i<=erat_len;i++)
      if(primep(i,erat))
	count++;
  printf("Pi(%lu)=%lu\n",erat_len,count);
  if(!(small_primes=(ptype *) malloc(sizeof(ptype)*count)))
    fatal_error("Failed to alloacte memory for small_primes.");
  small_primes[0]=2;
  for(i=3,j=1;j<count;i++)
    if(primep(i,erat))
      small_primes[j++]=i;
  
  sieve(target,target_start,target_len,small_primes,count);

  for(source_start=erat_len+1;source_start<sqrt_end;source_start+=source_len)
    {
      set_primes(source,source_len);
      //printf("sieving source %lu\n",source_start);
      sieve(source,source_start,source_len,small_primes,count);
      //printf("sieving target\n");
      sieve1(target,target_start,target_len,source,source_start,source_len);
    }
  
  for(i=1,count=0;i<=target_len;i++)
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
