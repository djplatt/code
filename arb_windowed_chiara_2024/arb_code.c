#include <inttypes.h>
#include <flint/fmpz.h>
#include <flint/acb_dirichlet.h>

#define N ((slong) 65536)
#define A ((slong) 16)
#define B (N/A)

int main()
{

  fmpz_t T;
  fmpz_init(T);
  fmpz_set_ui(T,100000000000ll);
  arb_t res[N];
  for(uint64_t i=0;i<N;i++)
    arb_init(res[i]);
  
  acb_dirichlet_platt_scaled_lambda_vec(res[0],T,A,B,50);

  return 0;
}
