#include "stdio.h"
#include "stdlib.h"
#include "inttypes.h"
#include "math.h"

/* Integer square root by Halleck's method, with Legalize's speedup */
inline uint64_t isqrt (uint64_t x){
  uint64_t   squaredbit, remainder, root;

   if (!x) return 0;
  
   /* Load the binary constant 01 00 00 ... 00, where the number
    * of zero bits to the right of the single one bit
    * is even, and the one bit is as far left as is consistant
    * with that condition.)
    */
   squaredbit  = (uint64_t) 0x4000000000000000;
   //((((uint64_t) ~0L) >> 1) & 
   //                   ~(((uint64_t) ~0L) >> 2));
   /* This portable load replaces the loop that used to be 
    * here, and was donated by  legalize@xmission.com 
    */

   /* Form bits of the answer. */
   remainder = x;  root = 0;
   while (squaredbit > 0) {
     if (remainder >= (squaredbit | root)) {
         remainder -= (squaredbit | root);
         root >>= 1; root |= squaredbit;
     } else {
         root >>= 1;
     }
     squaredbit >>= 2; 
   }

   return root;
}

// OK so long as x<=53 bits
uint64_t isqrt2(uint64_t x)
{
  double y=x;
  y=floor(sqrt(y));
  return(y);
}

inline uint64_t issquare1(uint64_t x)
{
  uint64_t is=isqrt(x);
  if(is*is==x)
    return is;
  else
    return 0;
}

uint64_t issquare2(uint64_t n)
{
    int h = n & 0xF; // last hexidecimal "digit"
    if (h > 9)
        return 0; // return immediately in 6 cases out of 16.

    // Take advantage of Boolean short-circuit evaluation
    if ( h != 2 && h != 3 && h != 5 && h != 6 && h != 7 && h != 8 )
    {
      // take square root if you must
      uint64_t t = (uint64_t) floor( sqrt((double) n));
      if(t*t == n)
	return t;
      else
	return 0;
    }
    return 0;
}
/*
inline uint64_t issquare2(uint64_t n)
{
  if((n&3)>1)
    return 0;
  if(__builtin_popcountl(n%7)!=1)
    return 0;
  uint64_t t = (uint64_t) floor( sqrt((double) n));
  return t*t == n;
} 
*/
#define LIM1 10000LL

int main(int argc, char**argv)
{
  uint64_t a,b,c,r;

  for(a=1;a<=LIM1;a++)
    for(b=a+2;b<=LIM1;b++)
      {
	//a=atol(argv[1]);b=atol(argv[2]);
	r=issquare1(a*b+1);
	if(r)
	  {
	    c=a+b-(r<<1);
	    if((!issquare1(a*c+1))||(!issquare1(b*c+1)))
	      {
		printf("Failed with %lu %lu %lu\n",a,b,c);
		exit(0);
	      }
	    printf("%lu %lu %lu\n",a,b,c);
	    fflush(stdout);
	  }
      }
  return 0;
}
