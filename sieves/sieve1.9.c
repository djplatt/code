#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"
#include "../includes/pi_x.h"

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

// returns x%y
// x has at most 80 significant bits (1.2e24)
// y has at most 48 (2.8e14)
inline __uint64_t rem_80_48(__uint128_t x, __uint64_t y)
{
  __uint64_t w=x&0xFFFF,z=x>>16;
  return((((z%y)<<16)+w)%y);
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


// sieving using base primes 2,3,5... held in small_primes
// all primes here are "small" <=x^(0.25)
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
  printf("in sieve1 with source_start=%lu\n",source_start);
  for(p=1,q=p+source_offset;p<=source_len;p+=2,q+=2) // skip by 2 as all primes are odd
    if(primep(p,source))
      {
	//printf("prime = %lu\n",q);
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

int main(int argc, char **argv)
{
  ptype *erat,*source,*target,*small_primes;
  long int sieve_num,num_its;
  ptype i,j,erat_len,source_len=SOURCE_LEN,target_len=TARGET_LEN,it;
  bigint target_start,del2_sum,X0;
  ptype count,source_start,sqrt_end,tl,tl2,segs_per_sieve=SEGS_PER_SIEVE,pix;
  __int64_t del,del2,del_sum;
  FILE *outfile;

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

  erat_len=sqrt(sqrt(X0));
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
  for(i=3,pix=1;i<=erat_len;i+=2)
      if(primep(i,erat))
	{
	  //printf("%lu is prime\n",i);
	  pix++;
	}
  printf("Pi(%lu)=%lu\n",erat_len,pix);
  if(!(small_primes=(ptype *) malloc(sizeof(ptype)*pix)))
    fatal_error("Failed to alloacte memory for small_primes.");
  small_primes[0]=2;
  for(i=3,j=1;j<pix;i+=2)
    if(primep(i,erat))
      small_primes[j++]=i;



  target_start=X0-target_len*(NUM_SIEVES/2-sieve_num)+1;
  for(it=sieve_num;it<num_its+sieve_num;it++,target_start+=target_len)
    {
      printf("Running sieve starting at ");print_bigint(target_start);printf("\n");
      fwrite(&it,sizeof(ptype),1,outfile);
      fwrite(&target_start,sizeof(bigint),1,outfile);
      // remove the small primes from the target      
      set_primes(target,target_len);
      sieve(target,target_start,target_len,small_primes,pix);

      sqrt_end=sqrt(target_start);

      for(i=0,source_start=erat_len+1;source_start<sqrt_end;i++,source_start+=source_len)
	{
	  // remove small primes from source
	  set_primes(source,source_len);
	  sieve(source,source_start,source_len,small_primes,pix);
	  // remove source primes from target
	  sieve1(target,target_start,target_len,source,source_start,source_len);
	}
      
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
