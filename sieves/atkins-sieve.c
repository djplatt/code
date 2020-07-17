// as yet this is bugged due to conversion
// of __int128_t to doubles which is a
// rounded conversion e.g.
// x=n;
// x=sqrt(x);
// may yield x^2 in [n-1,n+1] (I think)
// need to fix atkins_sieve(...) which
// assumes x^2 in [n-1,n]
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

#define primep(p,v) v[p>>(log_num_bits+1)]&mask1[p&all_ones]
#define clear_prime(p,v) v[p>>(log_num_bits+1)]&=mask[p&all_ones]
#define toggle_prime(p,v) v[p>>(log_num_bits+1)]^=mask1[p&all_ones]

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
  long int i;
  //printf("in set_primes with len=%lu = %lu words.\n",len,len>>(log_num_bits+1));
  for(i=0;i<=(len>>(log_num_bits+1));i++)
    {
      vec[i]=0;
      vec[i]=~vec[i];
    }
}
inline void unset_primes(ptype *vec, const ptype len)
{
  long int i;
  //printf("in set_primes with len=%lu = %lu words.\n",len,len>>(log_num_bits+1));
  for(i=0;i<=(len>>(log_num_bits+1));i++)
      vec[i]=0;
}

// note we don't trust the first entry of target (unless it is np^2)
// because we always miss A=x^2+y^2 etc. due to rounding of 128 bits into
// double prior to square root.
void atkins_sieve (ptype *target, const bigint target_start, const ptype target_len)
{
  ptype step_count=0;
  bigint n,x2,twox2,y2,target_offset=target_start-1,
    target_end=target_offset+target_len;
  ptype xlim,xmin,x,y,ptr,ymax,inc;
  double xd;

  //printf("target_len=%lu\n",target_len);
  //printf("target_end=");print_bigint(target_end);printf("\n");

  printf("In atkins sieve for ");
  print_bigint(target_start+2);
  printf(" to ");
  print_bigint(target_end);
  printf("\n");

  /*  
  printf("doing x^2+y^2=1 mod 4, x>y>0, => x!=y mod 2\n");
  x=sqrt(target_start>>1);
  x2=x*x;
  if(x2+1<target_start)
    {
      x++;
      x2+=(x<<1)-1;
    }

  xlim=sqrt(target_end-1);
  for(;x<=xlim;x++,x2+=(x<<1)-1)
    {
      if(x2>=target_start)
	y=1;
      else
	y=sqrt(target_start-x2)+1;
      ymax=sqrt(target_end-x2);
      
      if(((x^y)&1)==0)
	y++;
      
      if(ymax>=x)
	ymax=x-1;
      y2=y;
      y2*=y;
      n=x2+y2;
      for(;y<=ymax;)
	{
	  if(n>target_end)
	    break;
	  if((n&3)==1)
	    {
	      ptr=n-target_offset;
	      toggle_prime(ptr,target);
	    }
	  y+=2; // maintain relative parities
	  inc=(y<<2)-4;
	  y2+=inc;
	  n+=inc;
	}
    }
  
  //printf("Stepped on %ld times.\n",step_count);
  printf("doing 2x^2+y^2=3 mod 8, x>0, y>0 => x=y=1 mod 2\n");
  xlim=sqrt((target_end>>1)-1);
  for(x=1,twox2=2;x<=xlim;)
    {
      //printf("x=%ld\n",x);
      y=sqrt(target_start-twox2)+1;
      if(!(y&1))
	y++;
      y2=y;
      y2*=y;
      n=twox2+y2;
      for(;n<=target_end;)
	{
	  //printf("y=%ld\n",y);
	  //printf("n=");print_bigint(n);printf("\n");
	  if((n&7)==3)
	    {
	      ptr=n-target_offset;
	      toggle_prime(ptr,target);
	    }
	  y+=2;
	  inc=(y<<2)-4;
	  y2+=inc;
	  n+=inc;
	}
      x+=2;
      twox2+=(x<<3)-8;
    }
  */

  printf("doing 2x^2-y^2=7 mod 8, x>y>0 => x=0, y=1 mod 2\n");
  x=sqrt((target_start+1)>>1)+1;
  if((x&1)==1)
    x++;
  xlim=sqrt(target_end+2)-1;
  twox2=x;
  twox2*=x;
  twox2<<=1;
  for(;x<xlim;)
    {
      xd=twox2-target_start;
      xd=sqrt(xd);
      y=xd;
      if(y>=x)
	y=x-1;
      else if ((y&1)==0)
	y--;
      y2=y*y;
      n=twox2-y2;
      if((x>=2505394894))
	{
	  printf("x=%ld y=%ld xd=%20.18e\n",x,y,xd);
	  printf("twox2=");
	  print_bigint(twox2);
	  printf("\n");
	}
      for(;y>0;)
	{
	  if((x==2628169654))
	    printf("x=%ld y=%ld\n",x,y);

	  if(n>target_end)
	    break;
	  if((n&7)==7)
	    {
	      ptr=n-target_offset;
	      toggle_prime(ptr,target);
	    }
	  y-=2;
	  inc=(y<<2)+4;
	  n+=inc;
	}
      x+=2;
      twox2+=(x<<3)-8;
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

// sieving using base primes 2,3,5... held in small_primes
inline void sieve (ptype *target, bigint target_start, const ptype target_len, const ptype *small_primes, const ptype pix)
{
  ptype i,p,q,ptr;
  ptype rem;
  ptype max_p=sqrt(target_start+target_len-1);
  //printf("in sieve with max_p=%lu\n",max_p);
  for(i=1;i<pix;i++) // start at 3
      {
	p=small_primes[i]; // 3,5,7...
	if(p>max_p)
	  return;
	q=p<<1;
	rem=target_start%p;
	if(rem==0) rem=p;
	ptr=p-rem+1;
	if((ptr&1)==0)
	  ptr+=p;
	//printf("ptr=%lu\n",ptr);

	for(;ptr<=target_len;ptr+=q) // skip by 2p each time
	  {
	    //printf("ptr=%lu\n",ptr);
	    //printf("clearing %lu at %lu from target.\n",ptr+target_start-1,ptr);
	    clear_prime(ptr,target);
	  }
      }
  //num_prims(target,target_start,target_len);
}
// sieving using base primes 2,3,5... held in small_primes
inline void sieve_small_squares (ptype *target, bigint target_start, const ptype target_len, const ptype *small_primes, const ptype pix)
{
  ptype i,ptr;
  ptype twop2,p2; // p^2 will not exceed 64 bits
  ptype rem;
  for(i=3;i<pix;i++) // start at 7
      {
	p2=small_primes[i]*small_primes[i]; // 3,5,7...
	twop2=p2+p2;
	rem=target_start%p2;
	if(rem==0) rem=p2;
	ptr=p2-rem+1;
	if((ptr&1)==0)
	  ptr+=p2;
	if(ptr>target_len)
	  continue;

	for(;ptr<=target_len;ptr+=twop2)
	  clear_prime(ptr,target);
      }
}

// primes to sieve are all odd
inline void sieve_squares (ptype *target, bigint target_start, const ptype target_len, const ptype *source, const ptype source_start, const ptype source_len)
{
  ptype p;
  bigint rem,q,ptr,twop2;
  ptype source_offset=source_start-1;
  for(p=1;p<=source_len;p+=2) // skip by 2 as all primes are odd
    if(primep(p,source))
      {
	q=p+source_offset; // q is the prime
	q*=q; // might be > 64 bits
	rem=target_start%q; // ditto
	if(rem==0) rem=q;
	ptr=q-rem+1;
	if(!(ptr&1))
	  ptr+=q;
	//printf("ptr=%lu\n",ptr);
	// there is no guarentee that this ptr will fall within target
	if(ptr>target_len)
	  continue;
	twop2=q+q;
	for(;ptr<=target_len;ptr+=twop2) // skip by 2q^2
	  clear_prime(ptr,target);
      }
}

int main()
{
  ptype *erat,*source,*target,*small_primes;
  ptype i,j,erat_len,source_len=4444440,target_len=4000000000;
  bigint target_start;
  ptype count=0,source_start,sqrt_end;


  for(i=0,target_start=1;i<19;i++,target_start*=10);
  target_start+=1;
  sqrt_end=sqrt(target_start+target_len-1);
  erat_len=sqrt(sqrt_end);
  if(erat_len&1)
    erat_len++;
  printf("erat length = %lu\n",erat_len);
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

  unset_primes(target,target_len);

  atkins_sieve(target,target_start,target_len); // this leaves 3's,5's and mp^2
  
  sieve_small_squares(target,target_start,target_len,small_primes,count);
  

  for(source_start=erat_len+1;source_start<sqrt_end;source_start+=source_len)
    {
      // remove small primes from source
      set_primes(source,source_len);
      
      sieve(source,source_start,source_len,small_primes,count);
      
      // remove source prime squares from target
      sieve_squares(target,target_start,target_len,source,source_start,source_len);
    }

  // remove 3's
  i=target_start%3;
  if(i==0) i=3;
  i=4-i;
  if(!(i&1))
    i+=3;

  for(;i<=target_len;i+=6)
    clear_prime(i,target);

  // remove 5's
  i=target_start%5;
  if(i==0) i=5;
  i=6-i;
  if(!(i&1))
    i+=5;

  for(;i<=target_len;i+=10)
    clear_prime(i,target);

  // don't trust the first number due to sqrt rounding
  for(i=3,count=0;i<=target_len;i+=2)
    if(primep(i,target))
      {
	count++;
	print_bigint(target_start+i-1);printf(" is prime.\n");
	if(count==40)
	  return(0);
      }
  printf("Pi(");print_bigint(target_start+target_len-1);printf(")-Pi(");
  print_bigint(target_start+2);printf(")=%lu\n",count);
  
  return(0);
}
