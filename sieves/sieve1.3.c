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

#define primep(p,v) (v[p>>(log_num_bits+1)]&mask1[p&all_ones])
#define clear_prime(p,v) (v[p>>(log_num_bits+1)]&=mask[p&all_ones])
#define many_ones ((ptype) ~0)

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

ptype mask1[num_bits<<1],mask[num_bits<<1];

inline void set_primes(ptype *vec, const ptype len)
{
  long unsigned int i;
  //printf("in set_primes with len=%lu = %lu words.\n",len,len>>(log_num_bits+1));
  for(i=0;i<=(len>>(log_num_bits+1));i++)
    vec[i]=many_ones;
}

// create base primes 2,3,5,...
inline void erat_sieve (ptype *erat, const ptype len)
{
  long unsigned int i,j,sqrt_len;
  for(i=1,j=1;i<(num_bits<<1);i+=2,j<<=1)
    {
      mask1[i]=j;
      mask[i]=0;
      mask[i]=~mask[i];
      mask[i]^=j;
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
inline void sieve (ptype *target, bigint target_start, const ptype target_len, const ptype *erat, const ptype erat_len)
{
  ptype p,two_p,ptr;
  ptype thirty_p,start,start_rem,two_p_rem;
  for(p=3;p<=erat_len;p+=2) // start at 3
    if(primep(p,erat))
      {
	two_p=p<<1;
	for(ptr=get_index(p,target_start);ptr<=target_len;ptr+=two_p)
	  clear_prime(ptr,target);
      }
}


ptype fill_Ws(ptype *erat) // assumes there are enough primes in erat
{
  ptype p,W;
  printf("in fill_Ws\n");
  p=3;
  W=3;
  for(p+=2;!primep(p,erat);p+=2);
  while(many_ones/W>p+1)
    {
      //printf("%lu/%lu=%lu vs p=%lu\n",many_ones,W[w],many_ones/W[w],p);
      //printf("muliplying by %lu\n",p);
      W*=p;
      for(p+=2;!primep(p,erat);p+=2);
    }
  return(W);
  //for(w=0;w<NO_WS;w++) printf("W[%lu] set to %lu\n",w,W[w]);
}

int main()
{
  ptype *erat,*target;
  ptype i,j,erat_len,target_len=4000000000;
  bigint target_start,p;
  ptype count,sqrt_end,A=2,max_p;
  ptype lastgcd,q,w,two_q,ptr,W;  


  for(i=0,target_start=1;i<18;i++,target_start*=10);
  target_start+=1;
  max_p=sqrt(target_start+target_len);
  sqrt_end=sqrt(target_start+target_len-1);
  erat_len=sqrt_end/A;
  if(!(erat_len&1))
    erat_len++;
  if(!(erat=(ptype *) malloc(sizeof(ptype)*(1+(erat_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for erat sieve.");

  if(!(target=(ptype *) malloc(sizeof(ptype)*(1+(target_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for target sieve.");

  printf("Doing sieve of Eratostenes up to %lu\n",erat_len);
  erat_sieve(erat,erat_len);


  for(i=3,count=1;i<=erat_len;i+=2)
      if(primep(i,erat))
	{
	  //printf("%lu is prime\n",i);
	  count++;
	}
  printf("Pi(%lu)=%lu\n",erat_len,count);

  // remove the small primes from the target
  set_primes(target,target_len);
  printf("Removing Eratostenes primes from target.\n");  
  sieve(target,target_start,target_len,erat,erat_len);

  /*
  for(ptr=1,count=0;ptr<=target_len;ptr++)
    if(primep(ptr,target))
      {
	count++;
	print_bigint(target_start+ptr-1);printf(" might be prime.\n");
	if(count==40)
	  return(0);
      }
  */


  W=fill_Ws(erat);

  printf("sieving target (from %lu to %lu) by wheel.\n",erat_len,max_p);
  for(q=erat_len;q<=max_p;q+=2)
    {
      if(gcd(q,W)!=1)
	continue;
      two_q=q<<1;
      ptr=get_index(q,target_start);

      if(ptr>target_len)
	continue;
      for(;ptr<=target_len;ptr+=two_q)
	clear_prime(ptr,target);
    }

  for(ptr=1,count=0;ptr<=target_len;ptr++)
    if(primep(ptr,target))
      {
	count++;
	if(count<=40)
	  {print_bigint(target_start+ptr-1);printf(" is prime.\n");}
      }
  printf("Pi(");print_bigint(target_start+target_len-1);printf(")-Pi(");
  print_bigint(target_start);printf(")=%lu\n",count);
  
  return(0);
}
