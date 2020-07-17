#define one_over_A ((double) 5.0/64.0)
#define LN_SMALL ((double) -230.0) // exp(LN_SMALL) is tiny but a normaised float 
#define TAYLOR_TERMS (30)  // how many terms to use in Taylor Series approx
#define MIN_M_FOR_FFT (50) // if no. f_hat terms required less than this, use simple sum
#define TARGET_ACC (100) // aiming for F_hat_err < exp(-TARGET_ACC) used 40 for main GRH run
//#define QT ((double) 1.0e7) // was 1.0e8 for the main event
double QT;
#define H ((double) 100.0)
