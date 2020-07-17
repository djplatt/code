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

inline sign_t sign(arb_t x)
{
  if(arb_is_negative(x))
    return(NEG);
  if(arb_is_positive(x))
    return(POS);
  return(UNK);
}

inline dir_t dir(arb_ptr left, arb_ptr right)
{
  if(arb_gt(left,right))
    return(DOWN);
  if(arb_lt(left,right))
    return(UP);
  return(UNK);
}

inline sign_t sign(acb_t x)
{
  if(arb_is_negative(acb_realref(x)))
    return(NEG);
  if(arb_is_positive(acb_realref(x)))
    return(POS);
  return(UNK);
}

inline dir_t dir(acb_ptr left, acb_ptr right)
{
  if(arb_gt(acb_realref(left),acb_realref(right)))
    return(DOWN);
  if(arb_lt(acb_realref(left),acb_realref(right)))
    return(UP);
  return(UNK);
}
/*
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
*/
void in_bytes(arb_t t, FILE *infile, int64_t prec)
{
  static bool init=false;
  static arb_t err;
  if(!init)
    {
      init=true;
      arb_init(err);
      arb_set_ui(err,1);
      arb_mul_2exp_si(err,err,-(int64_t)(OP_ACC+1));
    }
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
  arb_set_ui(t,c);
  arb_mul_2exp_si(t,t,32);
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64);
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-(int64_t)OP_ACC);
  arb_add_error(t,err);
}

double arb_get_d(arb_t x)
{
  static bool init=false;
  static arb_t tmp1;
  static mag_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      mag_init(tmp);
    }
  arb_get_mid_arb(tmp1,x);
  arb_get_mag(tmp,tmp1);
  if(arb_is_negative(tmp1))
    return(-mag_get_d(tmp));
  else
    return(mag_get_d(tmp));
}

inline void arb_mul_d(arb_t x, arb_t y, double d, int64_t prec)
{
  static bool init=false;
  static arb_t dd;
  if(!init)
    {
      arb_init(dd);
      init=true;
    }
  arb_set_d(dd,d);
  arb_mul(x,y,dd,prec);
}


inline void arb_div_d(arb_t x, arb_t y, double d, int64_t prec)
{
  static bool init=false;
  static arb_t dd;
  if(!init)
    {
      arb_init(dd);
      init=true;
    }
  arb_set_d(dd,d);
  arb_div(x,y,dd,prec);
}

inline void arb_add_d(arb_t x, arb_t y, double d, int64_t prec)
{
  static bool init=false;
  static arb_t dd;
  if(!init)
    {
      arb_init(dd);
      init=true;
    }
  arb_set_d(dd,d);
  arb_add(x,y,dd,prec);
}

inline void arb_sub_d(arb_t x, arb_t y, double d, int64_t prec)
{
  static bool init=false;
  static arb_t dd;
  if(!init)
    {
      arb_init(dd);
      init=true;
    }
  arb_set_d(dd,d);
  arb_sub(x,y,dd,prec);
}

int arb_cmp(arb_t x, arb_t y)
{
  if(arb_gt(x,y))
    return(1);
  if(arb_lt(x,y))
    return(-1);
  if(arb_eq(x,y))
    return(0);
  printf("fatal error in arb_cmp. Arguments overlap. Exiting.\n");
  exit(0);
}
#endif 
