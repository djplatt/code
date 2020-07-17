#include "stdio.h"
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

inline int issquare1(uint64_t x)
{
  uint64_t is=isqrt2(x);
  return(is*is==x);
}

uint64_t issquare2(uint64_t n)
{
    int h = n & 0xF; // last hexidecimal "digit"
    if (h > 9)
        return 0; // return immediately in 6 cases out of 16.

    // Take advantage of Boolean short-circuit evaluation
    if ( h != 2 && h != 3 && h != 5 && h != 6 && h != 7 && h != 8 )
    {
      if(__builtin_popcountl(n%7l)!=1)
	return 0;
      // take square root if you must
      uint64_t t = (uint64_t) floor( sqrt((double) n));
      return t*t == n;
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
#define LIM1 64LL
#define LIM2 ((uint64_t) LIM1*LIM1*LIM1*LIM1*LIM1)

int main()
{
  uint64_t a,b,c;

  //for(a=1;a<=LIM1;a++)
  //for(b=a+2;b<=LIM1;b++)
  a=1;b=3;
      if(issquare2(a*b+1))
	{
	  c=b*b;
	  c*=c;c*=b; //b^5
	  for(;c<=LIM2;c++)
	    if(issquare2(b*c+1))
	      if(issquare2(a*c+1))
		printf("%lu %lu %lu\n",a,b,c);
	}
  return 0;
}
