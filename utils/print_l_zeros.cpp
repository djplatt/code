//
// Read an LMFDB file containing zeros of Dirichlet L-functions
//
// DJ Platt Oct/2026
//
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
/*
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
*/
zero_t operator + (zero_t &lhs, const zero_t &rhs)
{
  zero_t res;

  if(__builtin_uaddl_overflow(lhs.a,rhs.a,&res.a))
    res.b=lhs.b+rhs.b+1;
  else
    res.b=lhs.b+rhs.b;
  res.c=lhs.c+rhs.c;
  return res;
}

// extended GCD
void XGCD(int64_t *d, 
	  int64_t *s, 
	  int64_t *t, 
	  int64_t a, 
	  int64_t b)
{
   int64_t  u, v, u0, v0, u1, v1, u2, v2, q, r;

   int64_t aneg = 0, bneg = 0;

   if (a < 0) {
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d[0] = u;
   s[0] = u1;
   t[0] = v1;
}
   
// taken from NTL:
// compute a^-1 mod n
int64_t InvMod(int64_t a, int64_t n)
{
   int64_t d, s, t;

   XGCD(&d, &s, &t, a, n);
   if (d != 1) return -1;
   if (s < 0)
      return s + n;
   else
      return s;
}

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

  uint64_t q;
  if(fread(&q,sizeof(uint64_t),1,zfile)!=1)
    {
      printf("Error reading conductor. Exiting.\n");
      exit(0);
    }
  printf("Doing modulus %lu\n",q);
  zero_t gamma_0=zero_t(0,0,0);
  if(q==3) // these are offset by 8 because first zero is so high
    gamma_0.c=8*32;
    
  uint64_t ind;
  while(fread(&ind,sizeof(uint64_t),1,zfile)==1)
    {
      printf("Doing index %lu\n",ind);
      uint64_t nzeros;
      if(fread(&nzeros,sizeof(uint64_t),1,zfile)!=1)
	{
	  printf("Error reading zeros count. Exiting.\n");
	  exit(0);
	}
      printf("Going to read %lu zeros.\n",nzeros);
      zero_t gap,gamma;
      gamma=gamma_0;
      for(uint64_t z=0;z<nzeros;z++)
	{
	  gap=get_zero(zfile);
	  gamma=gamma+gap;
	  //
	  // at this point, gamma is the mid point of an interval
	  // of width 2^(-101) containing the ordinate of a zero
	  // For rigorous work, convert it into an ARB or similar
	  // and add an error of +/- 2^(-102)
	  // We will just approximate it with a double
	  //
	  printf("There is a zero near %10.8e.\n",(double) gamma);
	}
      uint64_t exp_ind=InvMod(ind,q);
      if(ind==exp_ind) // a real character
	continue;
      if(fread(&ind,sizeof(uint64_t),1,zfile)!=1)
	{
	  printf("Error reading compositive index. Exiting.\n");
	  exit(0);
	}
      if(ind!=exp_ind)
	{
	  printf("Index found (%lu) not as expected (%lu). Exiting.\n",ind,exp_ind);
	  exit(0);
	}
      printf("Doing composite index %lu\n",ind);
      if(fread(&nzeros,sizeof(uint64_t),1,zfile)!=1)
	{
	  printf("Error reading zeros count. Exiting.\n");
	  exit(0);
	}
      printf("Going to read %lu zeros.\n",nzeros);
      gamma=gamma_0;
      for(uint64_t z=0;z<nzeros;z++)
	{
	  gap=get_zero(zfile);
	  gamma=gamma+gap;
	  printf("There is a zero near %10.8e.\n",(double) gamma);
	}
    }
  fclose(zfile);
  return 0;
}
