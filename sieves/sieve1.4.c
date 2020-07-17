#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"

// as sieve1.1 but working mod 30
// so get flags for 240 ints into one 64 bit ptype

#define debug printf("Reached line number %d.\n",__LINE__)

typedef long int ptype;
typedef __int128_t bigint;

// hardwired for mod 2*3*5=30
#define primep(p,v) (v[p/240]&mask1[p%240])
//#define primep(p,v) (v[p/(num_bits*2)]&mask1[p&all_ones])
#define clear_prime(p,v) (v[p/240]&=mask[p%240])
//#define clear_prime(p,v) (v[p/(num_bits*2)]&=mask[p&all_ones])

inline ptype gcd (ptype a, ptype b)
// Euclid algorithm gcd
// best if a<=b ?
{
  ptype c;
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

ptype mask1[240],mask[240];

inline void set_primes(ptype *vec, const ptype len)
{
  ptype i;
  //printf("in set_primes with len=%lu = %lu words.\n",len,len>>(log_num_bits+1));
  for(i=0;i<=((len-1)/240);i++)
    {
      vec[i]=0;
      vec[i]=~vec[i];
    }
}

// create base primes 2,3,5,...
inline void erat_sieve (ptype *erat, const ptype len)
{
  ptype i,j,sqrt_len;
  for(i=0,j=1;i<240;i++)
    {
      if(gcd(i,30)==1)
	{
	  mask1[i]=j;
	  mask[i]=0;
	  mask[i]=~mask[i];
	  mask[i]^=j;
	  j<<=1;
	}
      else
	{
	  mask1[i]=0;
	  mask[i]=0;
	  mask[i]=~mask[i];
	}
    }
  
  set_primes(erat,len);
  sqrt_len=sqrt(len);
  for(i=3;i<=sqrt_len;i+=2)
    if(primep(i,erat))
      for(j=i*i;j<=len;j+=i+i)
	clear_prime(j,erat);
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

int main()
{
  ptype *erat,*source,*target,*small_primes;
  ptype i,j,erat_len,source_len=4444440,target_len=30000000000;
  bigint target_start;
  ptype count,source_start,sqrt_end;

  while((source_len%30)!=0)
    source_len++; // must be multiple of 30
  for(i=0,target_start=1;i<22;i++,target_start*=10);
  while((target_start%30)!=1)
    target_start++; // must start at 1 mod 30
  sqrt_end=sqrt(target_start+target_len-1);
  erat_len=sqrt(sqrt_end);
  while((erat_len%30)!=0)
    erat_len++; // must be a multiple of 30

  if(!(erat=(ptype *) malloc(sizeof(ptype)*(1+(erat_len/240)))))
    fatal_error("Failed to allocate memory for erat sieve.");
  if(!(source=(ptype *) malloc(sizeof(ptype)*(1+(source_len/240)))))
    fatal_error("Failed to allocate memory for source sieve.");
  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len/240)))))
    fatal_error("Failed to allocate memory for target sieve.");
  erat_sieve(erat,erat_len);
  for(i=7,count=3;i<=erat_len;i+=2)
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
