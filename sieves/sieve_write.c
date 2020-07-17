#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "inttypes.h"

#define debug printf("Reached line number %d.\n",__LINE__)

typedef long unsigned int ptype;

#define num_bits (sizeof(ptype)<<3)
#define log_num_bits (6)
#define all_ones ((num_bits<<1)-1)

#define primep(p,v) (v[p>>(log_num_bits+1)]&mask1[p&all_ones])
//#define primep(p,v) (v[p/(num_bits*2)]&mask1[p&all_ones])
#define clear_prime(p,v) (v[p>>(log_num_bits+1)]&=mask[p&all_ones])
//#define clear_prime(p,v) (v[p/(num_bits*2)]&=mask[p&all_ones])

#define BUFF_SIZE (8192)
#define FNAME "primes.dat"

inline void fatal_error(const char *str)
{
  fputs(str,stderr);
  fputs(" Exiting.\n",stderr);
  abort();
}


unsigned char buff[BUFF_SIZE];
unsigned long int buff_ptr=0;
ptype last_prime=3;

// write two zeros and force write
void flush_buff(FILE *outfile)
{
  buff[buff_ptr++]=0;
  if(buff_ptr==BUFF_SIZE)
    {
      buff_ptr=0;
      fwrite(buff,sizeof(unsigned char),BUFF_SIZE,outfile);
    }
  buff[buff_ptr]=0;
  fwrite(buff,sizeof(unsigned char),BUFF_SIZE,outfile);
  fclose(outfile);
}

void write_prime(ptype p, FILE *outfile)
{
  ptype pdiff=(p-last_prime)>>1,pdiff2;
  if(pdiff>255) // write zero then two bytes for gap
    {
      printf("writing prime %lu.\n",p);
      printf("previous prime %lu\n",last_prime);
      printf("Difference between primes %lu\n",pdiff<<1);
      buff[buff_ptr++]=0;
      if(buff_ptr==BUFF_SIZE)
	{
	  buff_ptr=0;
	  fwrite(buff,sizeof(unsigned char),BUFF_SIZE,outfile);
	}
      pdiff2=pdiff>>8; // get msb
      if(pdiff2>255)
	fatal_error("Prime gap exceeded two bytes!.");
      buff[buff_ptr++]=(unsigned char) pdiff2; // msb
      if(buff_ptr==BUFF_SIZE)
	{
	  buff_ptr=0;
	  fwrite(buff,sizeof(unsigned char),BUFF_SIZE,outfile);
	}
      buff[buff_ptr++]=(unsigned char) pdiff; // lsb
      if(buff_ptr==BUFF_SIZE)
	{
	  buff_ptr=0;
	  fwrite(buff,sizeof(unsigned char),BUFF_SIZE,outfile);
	}
    }
  else
    {
      buff[buff_ptr++]=(unsigned char) pdiff;
      if(buff_ptr==BUFF_SIZE)
	{
	  buff_ptr=0;
	  fwrite(buff,sizeof(unsigned char),BUFF_SIZE,outfile);
	}
    }
  last_prime=p;
}

ptype mask1[num_bits<<1],mask[num_bits<<1];

inline void set_primes(ptype *vec, const ptype len)
{
  ptype i;
  //printf("in set_primes with len=%lu = %lu words.\n",len,len>>(log_num_bits+1));
  for(i=0;i<=(len>>(log_num_bits+1));i++)
    vec[i]=0XFFFFFFFFFFFFFFFFL;
}

// create base primes 2,3,5,...
inline void erat_sieve (ptype *erat, const ptype len)
{
  ptype i,j,sqrt_len;
  for(i=1,j=1;i<(num_bits<<1);i+=2,j<<=1)
    {
      mask1[i]=j;
      mask[i]=0XFFFFFFFFFFFFFFFFL;
      mask[i]^=j;
    }
  
  //printf("mask1[%lu]=%lX\nmask[%lu]=%lX\n",9LU,mask1[9],9LU,mask[9]);
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

inline ptype get_index(ptype p, ptype start)
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
inline void sieve (ptype *target, ptype target_start, const ptype target_len, const ptype *small_primes, const ptype pix)
{
  ptype i,p,q,ptr;
  ptype max_p=floor(sqrt(target_start+target_len-1));
  //printf("in sieve with max_p=%lu\n",max_p);
  for(i=1;i<pix;i++) // start at 3
      {
	p=small_primes[i]; // 3,5,7...
	
	if(p>max_p)
	  {
	    //printf("Returning at p=%ld\n",p);
	    return;
	  }
	
	ptr=get_index(p,target_start); // there must be some in here
	q=p<<1;
	for(;ptr<=target_len;ptr+=q) // skip by 2p each time
	    clear_prime(ptr,target);
      }
}


// hard wired for x=10^24, 100,000 iterations
int main()
{
  ptype *erat,*source,*target,*small_primes;
  ptype i,j,erat_len=1000000,source_len=9999990;
  ptype count,source_start,sqrt_end,ptr,p;
  FILE *outfile;

  if(!(outfile=fopen(FNAME,"wb")))
    fatal_error("Failed to open outfile.");

  for(sqrt_end=1,i=0;i<12;i++,sqrt_end*=10);

  if(erat_len&1)
    erat_len++; // make sure it's even
  if(source_len&1)
    fatal_error("source length must be even.\n");

  if(!(erat=(ptype *) malloc(sizeof(ptype)*(1+(erat_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for erat sieve.");
  if(!(source=(ptype *) malloc(sizeof(ptype)*(1+(source_len>>(log_num_bits+1))))))
    fatal_error("Failed to allocate memory for source sieve.");

  erat_sieve(erat,erat_len);
  for(i=5,count=2;i<=erat_len;i+=2)
      if(primep(i,erat))
	{
	  write_prime(i,outfile);
	  count++;
	}
  printf("Pi(%lu)=%lu\n",erat_len,count);

  if(!(small_primes=(ptype *) malloc(sizeof(ptype)*count)))
    fatal_error("Failed to alloacte memory for small_primes.");

  small_primes[0]=2;
  for(i=3,j=1;j<count;i+=2)
    if(primep(i,erat))
      small_primes[j++]=i;

  // will sieve to 10^12+10^6 using primes to 10^6
  // primes to 10^6 are sufficient
  for(i=0,source_start=erat_len+1;i<=100000;i++,source_start+=source_len)
    {
      if((i%1000)==0)
	printf("sieving from %lu to %lu\n",source_start,source_start+source_len*1000-1);
      // remove small primes from source
      set_primes(source,source_len);
      sieve(source,source_start,source_len,small_primes,count);
      for(ptr=1,p=source_start;ptr<source_len;ptr++,p++)
	  if(primep(ptr,source))
	    write_prime(p,outfile);
    }
  flush_buff(outfile);  
  return(0);
}
