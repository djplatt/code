#include <math.h>

long double
nextafterl (long double x, long double y)
{
  union {
      long double ld;
      struct {
        unsigned long long mantissa;
        unsigned short expn;
        unsigned short pad;
      } __attribute__ ((packed)) parts; 
  } u;

  /* The normal bit is explicit for long doubles, unlike
     float and double.  */
  static const unsigned long long normal_bit = 0x8000000000000000ull;

  if (isnan (y) || isnan (x))
    return x + y;

  if (x == y )
     /* nextafter (0.0, -O.0) should return -0.0.  */
     return y;

  u.ld = x;
  if (x == 0.0L)
    {
      u.parts.mantissa = 1ull;
      return y > 0.0L ? u.ld : -u.ld;
    }

  if (((x > 0.0L) ^ (y > x)) == 0)
    {
      u.parts.mantissa++;
      if ((u.parts.mantissa & ~normal_bit) == 0ull)
	u.parts.expn++;
    }
  else
    {
      if ((u.parts.mantissa & ~normal_bit) == 0ull)
	u.parts.expn--;
      u.parts.mantissa--;
    }

  /* If we have updated the expn of a normal number,
     or moved from denormal to normal, [re]set the normal bit.  */ 
     
  if (u.parts.expn & 0x7fff)
    u.parts.mantissa |=  normal_bit;

  return u.ld;
}
