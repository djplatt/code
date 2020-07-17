#ifndef PI_X
#define PI_X

#define PREC (200)
//#define LOG_10_X0 (22) // end of sieve is 10^log_10_x0
#define LOG_10_X0 (13)

//#define LAMBDA ((double) 7699934548453755.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/512.0) // for x = 10^22 - a bit
//#define LAMBDA ((double) 1.0/1024.0) // for x = 10^6
//#define LAMBDA ((double) 45497558225465.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0) // for x = 10^11
#define LAMBDA ((double)   4646592975839.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0) // for x = 10^13

#define XI ((long unsigned int) 1 << 31) // so that (1-XI)^2 < 64 bit
//#define NUM_SIEVES ((long unsigned int) 273747) // this many sieves either side for 10^22
#define NUM_SIEVES (1)
#define SEGS_PER_SIEVE ((long unsigned int) 8) // must be a power of two
#define TARGET_LEN ((long unsigned int) XI*SEGS_PER_SIEVE) // 16G = 1Gbyte
#define LOG_2_SOURCE_LEN (22)
#define SOURCE_LEN ((long unsigned int) 1<<LOG_2_SOURCE_LEN)

typedef __uint64_t ptype;
typedef __uint128_t bigint;

inline void print_bigint(bigint i)
{
  if(i<10)
    printf("%1lu",(long unsigned int) i);
  else
    {
      print_bigint(i/10);
      printf("%1lu",(long unsigned int) (i%10));
    }
}

bigint calc_x()
{
  bigint x;
  ptype i;

  for(i=0,x=1;i<LOG_10_X0;i++,x*=10);
  return(x);
  return(x-NUM_SIEVES*XI*2+1);
}

#endif
