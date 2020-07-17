#ifndef UPSAMDEFS
#define UPSAMDEFS

#ifndef MAX_Q
#define MAX_Q (300000)
#endif

#define POS 0
#define NEG 1
#define UNK 2
#define UP 0
#define DOWN 1
#define UP_SAMPLE_IN_WIDTH (512)
#define UP_SAMPLE_OUT_WIDTH (UP_SAMPLE_IN_WIDTH*UP_SAMPLE_RATE)
#define UP_SAMPLE_H (90.0) // controls decay of Gaussian, small h, fast decay
#define UP_SAMPLE_USABLE_WIDTH (UP_SAMPLE_IN_WIDTH/2)
#define UP_SAMPLE_SACRIFICE ((UP_SAMPLE_IN_WIDTH-UP_SAMPLE_USABLE_WIDTH)/2)
#define FFTW_PLAN_MODE FFTW_MEASURE // we only plan once, so get the best plan poss

#define INTER_H (7.0/32.0)
#define TWO_H_SQR (2.0*INTER_H*INTER_H)
#define twoB_num (64)
#define twoB_den (5)
#define one_over_two_B ((double) ((double) twoB_den/(double) twoB_num))

//#define one_over_two_B (5.0/64.0)
#define INTER_N (20)
//
// use this for q in[80001,100000]
// we want t0=1250 so Turing zone =[1250,1260]
// plus 30 in case
// we need to shift TZ
// for q<=100,000 t<=1280
//#define d_inter_err ((double) 6.53e-8) // new bound 6.48
//
// use this for q in [50001,80000]
// for q<=80,000 t <=2030
//#define d_inter_err ((double) 7.81e-8) // new bound 7.78
//
// use this for q in [40001,50000]
// for q<=50000 t<=2530
//#define d_inter_err ((double) 5.65e-8) // new bound 5.63
//
// use this for q in [30001,40000]
// for q<=40000 t<=3370
//#define d_inter_err ((double) 5.74e-8) // best bound
//
// use this for q in [25001,30000]
// t<=4030
//#define d_inter_err ((double) 4.95e-8) // best bound
//
// use this for q in [20001,25000]
// t<=5030
//#define d_inter_err ((double) 4.99e-8) // best bound
//
// use this for q in [16001,20000]
// t<=6280
//#define d_inter_err ((double) 4.83e-8) // best bound
//
// use this for q in [13101,16000]
// t<=7680
//#define d_inter_err ((double) 4.59e-8) // best_bound
//
// q in [10001,13100]
// t<=10040
// #define d_inter_err ((double) 4.75e-8) // best bound
// q in [8001,10000]
// t<=12560
// #define d_inter_err ((double) 4.42e-8) // best bound
//
// q in [6541,8000]
// t<15360
// #define d_inter_err ((double) 4.24e-8) // best bound
//
// q in [5001,6540]
// t<20040
//#define d_inter_err ((double) 4.40e-8) // best bound
// use this q in[3,100000], qt<=162400000
//#define d_inter_err ((double) 1.14e-7)
//
// q in [4001,5000]
// t<=25030
// #define d_inter_err ((double) 4.14e-8)
// q in [3001,4000]
// t<=33360
//#define d_inter_err ((double) 4.30e-8) // best bound
//
// q in [2001,3000]
// t<=50030
//#define d_inter_err ((double) 4.66e-8) // best bound
//
// q in [1001,2000]
// t<=100030
//#define d_inter_err ((double) 5.79e-8) // best bound
//
// q in [401,1000]
// t<= 249410
//#define d_inter_err ((double) 6.79e-8) // best bound
//
// q in [201,400]
// t<= 497550
//#define d_inter_err ((double) 5.84e-8) // best bound
//
// q in [101,200]
// t<=990130
//#define d_inter_err ((double) 6.13e-8) // best bound
//
// q in [69,100]
// t<=1449310
//#define d_inter_err ((double) 5.40e-8) // best bound
//
// q in [45,68]
//t<=2222260
//#define d_inter_err ((double) 5.87e-8) // best bound
// q in [24,47]
// t<=4166700
//#define d_inter_err ((double) 7.36e-8) // best bound
#ifdef __cplusplus
int_double d_inter_err;//=int_double(-1e10,1e10);
#endif

#define BAD_DIFF (INT_MAX)
#define FALSE_MAX (0)
#define FALSE_MIN (1)
#define CHECK_PROB (2)
#define RE_Z_0_PROB (3)
#define sign_t char
#define dir_t char


typedef struct
{
  double gap;
  double omega[4];
  //double im_s;

  unsigned int q;
  int index;
  unsigned int rate;
  unsigned int n0;
  unsigned int n02;
  unsigned int offset1;
  unsigned int offset2;
	
  char neg_one;
  char type;
  sign_t exp_sign;

} q_state;

void print_qs(q_state qs)
{
  printf("gap= %20.18e\nomega[0]=%20.18e\nomega[1]=%20.18e\nomega[2]=%20.18e\nomega[3]=%20.18e\n",qs.gap,qs.omega[0],qs.omega[1],qs.omega[2],qs.omega[3]);
  printf("im_s=%30.28e\nq=%d\nindex=%d\n",0.0/*qs.im_s*/,qs.q,qs.index);
  printf("n0=%d..%d\noffset=%d..%d\n",qs.n0,qs.n02,qs.offset1,qs.offset2);
}

inline void save_it (q_state qs, FILE *outfile)
{
//	printf("type=%d q=%d index=%d neg_one=%d n0=%d offset1=%d offset2=%d expected sign=%d\n",
//		qs.type,qs.q,qs.index,qs.neg_one,qs.n0,qs.offset1,qs.offset2,qs.exp_sign);
//	print_int_complex_str("omega=",qs.omega);
	fwrite(&qs,sizeof(q_state),1,outfile);
}
	
#endif
