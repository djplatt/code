#ifndef WIN_ZETA
#define WIN_ZETA

#include "inttypes.h"
#define OP_ACC (101) // Output is correct to +/- 2^(-OP_ACC-1)

// defined so that POS&NEG == 0 but no other combination
typedef int sign_t;
#define POS (1)
#define NEG (2)
#define UNK (3)

typedef int dir_t;
#define UP (1)
#define DOWN (2)

inline sign_t mpfr_sign(mpfr_ptr x)
{
  int s=mpfr_sgn(x);
  if(s<0)
    return(NEG);
  if(s>0)
    return(POS);
  return(UNK);
}

inline sign_t sign(mpfi_ptr x)
{
  if(mpfi_is_neg(x))
    return(NEG);
  if(mpfi_is_pos(x))
    return(POS);
  return(UNK);
}

inline dir_t mpfr_dir(mpfr_ptr left, mpfr_ptr right)
{
  int cp=mpfr_cmp(left,right);
  if(cp<0)
    return(UP);
  if(cp>0)
    return(DOWN);
  return(UNK);
}

inline dir_t dir(mpfi_ptr left, mpfi_ptr right)
{
  int cp=mpfi_cmp(left,right);
  if(cp<0)
    return(UP);
  if(cp>0)
    return(DOWN);
  return(UNK);
}

mpz_t out_z;
void init_out_bytes()
{
  mpz_init(out_z);
}

//long int ob_count=0;
void out_bytes(mpfr_ptr t, FILE *outfile)
{
  uint64_t a,i;
  uint32_t b;
  uint8_t c;

  mpfr_get_z(out_z,t,GMP_RNDN);
  a=mpz_get_ui(out_z); // just grabs least sig. long uns int (8)
  fwrite(&a,sizeof(uint64_t),1,outfile);
  mpz_fdiv_q_2exp(out_z,out_z,64);
  i=mpz_get_ui(out_z);
  b=i&0xFFFFFFFF;
  fwrite(&b,sizeof(uint32_t),1,outfile);
  i>>=32;
  if(i>255)
    {
      printf("Argument to out_bytes exceeds 13 bytes. Exiting.\n");
      exit(1);
    }
  c=i;
  fwrite(&c,sizeof(uint8_t),1,outfile);
}

mpfi_t pm1;
void init_in_bytes()
{
  mpfi_t tmp;
  mpfi_init(tmp);
  mpfi_init(pm1);
  mpfi_set_ui(pm1,1);
  mpfi_neg(tmp,pm1);
  mpfi_put(pm1,tmp);
  mpfi_div_2ui(pm1,pm1,OP_ACC+1);
  mpfi_clear(tmp);
}

void in_bytes(mpfi_ptr t, FILE *infile)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
  int res;

  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  mpfi_set_ui(t,c);
  mpfi_mul_2ui(t,t,32);
  mpfi_add_ui(t,t,b);
  mpfi_mul_2ui(t,t,64);
  mpfi_add_ui(t,t,a);
  mpfi_div_2ui(t,t,OP_ACC);
}
#endif 
