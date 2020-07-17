#include "inttypes.h"
#include "acb.h"

// structure describing L function family
// see p 387 col 2 of Booker 2006
// and conditions for Lemma 5.2, 5,5
typedef struct{
  //double dc;
  uint64_t r; // number of Gamma_R 0 means error
  double *mus; // the r mu_j's
  acb_t mu; // -1/2+1/r(1+sum mus[j])
  //int64_t m; // number of poles, all on Re(s)=1 -1 means error
  //arb_t *lambdas; // the imaginary parts of the poles
  //acb_t *residues; // residues at those poles
  arb_t *nus; // nu[j]=(Re mu[j]-1)/2
  double C; 
  double alpha; // |a_n| <= Cn^alpha
  arb_t c; // Re mu +1/2+alpha
  arb_t c_dash; // max(cr/2-1,0)
  double one_over_B;
  arb_t two_pi_by_B; // spacing of the Gs in the file
  int64_t low_i; // lowest value of i for which G^(k)(2 pi i/max_B) is present
                 // = round(log(1/sqrt(N))*max_B/2/pi)
  int64_t hi_i; // last i present
  uint64_t max_K; // max G^(k) provided
  arb_t **Gs; // the values of G^(k)(i/B)/k! read to be used
  uint64_t G_len; // number of G values to use =N/2-offset+1
} L_family_t;

// structure with data for this particular L function
typedef struct{
  uint64_t N; // conductor 0 means error
  double dc;
  bool self_dual_p; // L(1/2-it)=L(1/2+it)
  uint64_t M;
  acb_t *ans; // the coefficients a_n/sqrt(n)
  arb_t sum_ans; // sum over |an/sqrt(n)|
  acb_t epsilon; // the square root of the root number
  // we'll compute epsilons from the output.
} L_func_t;

// structure with parameters for the computation
typedef struct{
  uint64_t N; // length of FFT we are going for a positive power of 2 please
  int64_t offset; // lowest value of i for which G^(k)(i/B) is needed
  acb_t *w; // e(i/N) for i=0..N/2 the DFT
  double A; // N/B
  arb_t arb_A;
  arb_t one_over_A;
  uint64_t K; // truncation bound for Taylor approx 5-9 and 5-10
  double eta; // in (0,1) ARB does (-1,1)
  arb_t delta; // pi/2*(1-eta) (|eta| in ARB but see above)
  int64_t *ms; // ms[m]=round(log(m/sqrt(N))*B)
  arb_t *ums; // ms[m]/B
  arb_t *sks; // (log(m/sqrt(N))-um) m=1..M
  acb_t *G,*kres,*res,*skm;
} L_comp_t;

typedef struct{
  arb_t lem54;
  arb_t lem56;
  arb_t lem57;
  arb_t eq59;
} L_error_t;


void print_family(L_family_t *L)
{
  printf("# Gamma factors=%lu\n",L->r);
  for(uint64_t i=0;i<L->r;i++)
    printf("   mu_%lu=%10.8e\n",i,L->mus[i]);
  for(uint64_t i=0;i<L->r;i++)
    {
      printf("   nu_%lu=",i);
      arb_printd(L->nus[i],10);
      printf("\n");
    }
  printf("mu=");acb_printd(L->mu,10);
  /*
  printf("\n#poles=%lu\n",L->m);
  for(uint64_t i=0;i<L->m;i++)
    {
      printf("   lambda_%lu=",i);
      arb_printd(L->lambdas[i],10);
      printf(" with residue ");
      acb_printd(L->residues[i],10);
      printf("\n");
    }
  */
  printf("C=%10.8e\nalpha=%10.8e\n",L->C,L->alpha);
  printf("c=");arb_printd(L->c,10);
  printf("\nc'=");arb_printd(L->c_dash,10);printf("\n");
}

typedef struct{
  arb_t H;
  arb_t inv_2H2; // -1.0/2H^2
  uint64_t N; // number of samples to use each side of t
  arb_t *exps; // exp(-(n/A)^2/(2H^2)) n= 0 .. N-1

} L_upsample_t;
