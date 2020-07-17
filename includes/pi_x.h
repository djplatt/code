#ifndef PI_X
#define PI_X

#define PREC (200)

//#define LAMBDA ((double) 1.0/1024.0) // for x = 10^6
//#define LAMBDA ((double) 45497558225465.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0) // for x = 10^11

#define TEN_24
//#define TEN_23
//#define TEN_22
//#define TEN_20
//#define TEN_20A
//#define TEN_20B
//#define TEN_20C
//#define TEN_20D
//#define TEN_20E
//#define TEN_15
//#define TEN_13

// lambda for T1=20,950,046,000
#ifdef TEN_24
#define LOG_10_X0 (24) // end of sieve is 10^log_10_x0
#define LAMBDA ((double) 6273445730170391.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/16.0) // *2^(-84)

// change 2 back to 10 or something like it
#define XI ((long unsigned int) 1 << (31-10)) // so that (1-XI)^2 < 64 bit
#define NUM_SIEVES ((long unsigned int) 352800) // number of 2^34 width sieves to do, 1/2 each side of X0
#define SEGS_PER_SIEVE ((long unsigned int) 8<<10) // each big sieve is XI*SEGS wide
#endif

#ifdef TEN_23
#define LOG_10_X0 (23)
#define LAMBDA ((double) 6224003264759175.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/8.0) // 2^-83
#define XI ((long unsigned int) 1 << (31-10)) // so that (1-XI)^2 < 64 bit
// each job does 240*7 2^34 sieves so this is the next even multiple up
#define NUM_SIEVES ((long unsigned int) 70560) // number of XI*SEGS=2^34 width sieves to do, 1/2 each side of X0
#define SEGS_PER_SIEVE ((long unsigned int) 8<<10) // each big sieve is XI*SEGS wide
#endif

#ifdef TEN_22
#define LOG_10_X0 (22) // end of sieve is 10^log_10_x0
#define LAMBDA ((double) 7879885163686215.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0) // 7699934548453755.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/512.0) // for x = 10^22 - a bit
#define XI ((long unsigned int) 1 << 31) // so that (1-XI)^2 < 64 bit
#define NUM_SIEVES ((double) 1000)//68436.75)
#define SEGS_PER_SIEVE ((long unsigned int) 8)
#endif

// used to run a short sieve, many zeros test
// need to reduce XI drastically to make Taylor error in phi OK
#ifdef TEN_20
#define LOG_10_X0 (20) // end of sieve is 10^log_10_x0
#define LAMBDA ((double) 5417744701091588.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/8.0) // 2^-83
#define XI ((long unsigned int) 1 << (31))//-10)) // so that (1-XI)^2 < 64 bit
#define NUM_SIEVES ((double) 144)
#define SEGS_PER_SIEVE ((long unsigned int) (8))//<<10))
#endif
// used to run a short sieve, many zeros test
// need to reduce XI drastically to make Taylor error in phi OK
#ifdef TEN_20A
#define LOG_10_X0 (20) // end of sieve is 10^log_10_x0
#define LAMBDA ((double) 4628846093743418.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/64.0) // 2^-76
#define XI ((long unsigned int) 1 << (31-10)) // so that (1-XI)^2 < 64 bit
#define NUM_SIEVES ((double) 7056)
#define SEGS_PER_SIEVE ((long unsigned int) 8<<10)
#endif

#ifdef TEN_20B
#define LOG_10_X0 (20) // end of sieve is 10^log_10_x0
#define LAMBDA ((double) 7523497452684932.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0) // 2^-80
#define XI ((long unsigned int) 1 << (31-10)) // so that (1-XI)^2 < 64 bit
#define NUM_SIEVES ((double) 7056)
#define SEGS_PER_SIEVE ((long unsigned int) 8<<10)
#endif

#ifdef TEN_20C
#define LOG_10_X0 (20) // end of sieve is 10^log_10_x0
#define LAMBDA ((double) 6032309500618333.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/4.0) // 2^-82
#define XI ((long unsigned int) 1 << (31-10)) // so that (1-XI)^2 < 64 bit
#define NUM_SIEVES ((double) 7056)
#define SEGS_PER_SIEVE ((long unsigned int) 8<<10)
#endif

#ifdef TEN_20D
#define LOG_10_X0 (20) // end of sieve is 10^log_10_x0
#define LAMBDA ((double) 5021307034833028.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/2.0) // 2^-81
#define XI ((long unsigned int) 1 << (31-10)) // so that (1-XI)^2 < 64 bit
#define NUM_SIEVES ((double) 7056)
#define SEGS_PER_SIEVE ((long unsigned int) 8<<10)
#endif

#ifdef TEN_20E
#define LOG_10_X0 (20) // end of sieve is 10^log_10_x0
#define LAMBDA ((double) 7524486719454186.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/2.0) // 2^-81
#define XI ((long unsigned int) 1 << (31-10)) // so that (1-XI)^2 < 64 bit
#define NUM_SIEVES ((double) 7056)
#define SEGS_PER_SIEVE ((long unsigned int) 8<<10)
#endif

#ifdef TEN_15
#define LOG_10_X0 (15)
#define LAMBDA ((double) 5069816557000022/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0/2.0) // for x=10^15
#define XI ((long unsigned int) 1<<20) // for 10^15
#define SEGS_PER_SIEVE ((long unsigned int) 1<<14)
#define NUM_SIEVES (2) // total sieve width=NUM_SIEVES*XI*SEGS, split either side of X
#endif

#ifdef TEN_13
#define LOG_10_X0 (13)
#define LAMBDA ((double)   4646592975839.0/1024.0/1024.0/1024.0/1024.0/1024.0/1024.0) // for x = 10^13
#define XI ((long unsigned int) 1<<16) // for 10^13
#define SEGS_PER_SIEVE ((long unsigned int) 1<<18)
#define NUM_SIEVES (2) // total sieve width=NUM_SIEVES*XI*SEGS, split either side of X
#endif


#define TARGET_LEN ((long unsigned int) XI*SEGS_PER_SIEVE) // 2^34=16G => 1Gbyte memory mod 2
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
  return(x); // no offset
  //return(x-NUM_SIEVES*XI*SEGS_PER_SIEVE/2+1);
}

#endif
