#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#define two_64 ((double) 18446744073709551616.0) // 2^64
#define two_32 ((double) 4294967296.0) // 2^32
#define two_101 ((double) 2535301200456458802993406410752.0) // 2^101

// we store the difference between ordinates of zeros
// each difference is encoded as unsigned 64 bit (a) 32 bit (b) and 8 bit (c)
// a difference is then (a+b<<64+c<<96)>>101

class zero_t{
public:
  uint64_t a; // note this may overflow when adding
  uint64_t b; // b and c will not unless we add unfeasibly many of them
  uint64_t c;

  zero_t ()
  {
  }

  zero_t(const uint64_t a1, const uint32_t b1, const uint64_t c1)
  {
    a=a1;b=b1;c=c1;
  }

  // convert to a double.
  operator double() const
  {
    double res=this->c*two_32;
    res+=this->b;
    res*=two_64;
    res+=this->a;
    res/=two_101;
    return res;
  }

    friend zero_t operator + (const zero_t &lhs, const zero_t &rhs);
};

// a is the only field that can overflow as we have less than 2^32
// entries per block. Do adc anyway.
  inline zero_t operator + (const zero_t &lhs, const zero_t &rhs)
{
  zero_t res;
__asm__("movq %3,%0\n\t"
        "addq %4,%0\n\t"
	"movq %5,%1\n\t"
        "adcq %6,%1\n\t"
        "movq %7,%2\n\t"
	"adcq %8,%2\n\t"
	:"=r" (res.a), "=r" (res.b), "=r" (res.c)
	: "m" (lhs.a), "m" (rhs.a), "m" (lhs.b), "m" (rhs.b), "m" (lhs.c), "m" (rhs.c)
	  :);      
 return res;
}
/*
zero_t operator + (zero_t &lhs, const zero_t &rhs)
{
  zero_t res;

  if(__builtin_uaddl_overflow(lhs.a,rhs.a,&res.a))
    {
      printf("overflow happened.\n");
      res.b=lhs.b+rhs.b+1;
    }
  else
    res.b=lhs.b+rhs.b;
  res.c=lhs.c+rhs.c;
  return res;
}
*/
// read a zero from zfile
zero_t get_zero(FILE *zfile)
{
  zero_t res;
  uint32_t b;
  uint8_t c;
  
  if(fread(&res.a,sizeof(uint64_t),1,zfile)!=1)
    {
      printf("Error reading 64 bit part of zero from file. Exiting.\n");
      exit(0);
    }

  if(fread(&b,sizeof(uint32_t),1,zfile)!=1)
    {
      printf("Error reading 32 bit part of zero from file. Exiting.\n");
      exit(0);
    }
  res.b=b; // convert to 64 bit

  if(fread(&c,sizeof(uint8_t),1,zfile)!=1)
    {
      printf("Error reading 8 bit part of zero from file. Exiting.\n");
      exit(0);
    }
  res.c=c; // convert to 64 bit

  return res;
}


int main(int argc, char **argv)
{

  if(argc!=2)
    {
      printf("Usage:- %s <zeros file>.\n",argv[0]);
      exit(0);
    }


  FILE* zfile=fopen(argv[1],"rb");
  if(!zfile)
    {
      printf("Error opening file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  uint64_t num_blocks;
  if(fread(&num_blocks,sizeof(uint64_t),1,zfile)!=1)
    {
      printf("Error reading number of blocks from zeros file. Exiting.\n");
      exit(0);
    }
  //printf("Going to do %ld blocks.\n",num_blocks);

  double st[2];
  uint64_t zs[2],z,it;
  zero_t sum_zeros;
  double t;
  uint64_t z_rec[12];
  for(it=0;it<num_blocks;it++)
    {
      sum_zeros=zero_t(0,0,0);
      fread(st,sizeof(double),2,zfile);
      fread(zs,sizeof(uint64_t),1,zfile);
      if(st[0]==0.0)
	{
	  printf("st[0] was 0.0. Exiting.\n");
	  exit(0);
	}
      //printf("Block starts at %e and ends at %e.\n",st[0],st[1]);
      fread(zs+1,sizeof(uint64_t),1,zfile);
      for(z=zs[0]+1;z<=zs[1];z++)
	{
	  zero_t del_t=get_zero(zfile); // exact
	  if((del_t.a==0)&&(del_t.b==0)&&(del_t.c==0))
	    {
	      printf("get_zero returned 0. Exiting.\n");
	      exit(0);
	    }

	  sum_zeros=sum_zeros+del_t;
	  t=st[0]+(double) sum_zeros;
	  //
	  // t is now (approx) the ordinate of the z'th zero
	  // do something with it.
	  //
	  printf("Zero number %ld is at %f + %20.18e.\n",z,st[0],(double) sum_zeros);
	  //
	}
    }
  
  fclose(zfile);
  return 0;
}
