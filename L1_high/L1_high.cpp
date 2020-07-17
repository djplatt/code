//
// l-func-lmfdb.cpp
//
// hi precision, rigorous computation of Dirichlet L-functions
//
// compute root number L(1) and L(1/2) for all primitive characters mod q
// output double precision and +/- 2^{-102} estimates
//
// Uses Bober's characters.h
//
#include <iostream>
#include <cstdlib>
#include "../characters-master/characters.h"
#include "../includes/mpfi_arb_hurwitz.h"
#include "../includes/mpfi_arb_psi.h"


using namespace std;

#define PREC (300) // mpfi working precision

mpfi_c_t ctemp,ctemp1;
mpfi_t itemp,itemp2,mpfi_2_pi;

mpfi_t o_x;
mpfr_t o_l;
mpz_t z_l,z_r;
bool out_initialised=false;

#define OP_ACC (101)

// write out_z rounded to double
// and to 2^{-OP_ACC}
// using 14 bytes
void out_bytes(mpz_t out_z, FILE *outfile)
{
  bool negative_p=false;
  if(mpz_sgn(out_z)<0)
    {
      negative_p=true;
      mpz_neg(out_z,out_z);
    }
  uint64_t a,i;
  uint32_t b;
  int16_t c;
  a=mpz_get_ui(out_z); // just grabs least sig. long uns int (8)
  fwrite(&a,sizeof(uint64_t),1,outfile);
  mpz_fdiv_q_2exp(out_z,out_z,64);
  i=mpz_get_ui(out_z);
  b=i&0xFFFFFFFF;
  fwrite(&b,sizeof(uint32_t),1,outfile);
  i>>=32;
  if(i>65535)
    {
      printf("Argument to out_bytes exceeds 14 bytes. Exiting.\n");
      exit(1);
    }
  if(negative_p)
    c=-i;
  else
    c=i;
  fwrite(&c,sizeof(int16_t),1,outfile);
  //printf("wrote %lu %u %d\n",a,b,c);
}


void out_real(mpfi_t x, FILE* outfile)
{
  //mpfi_print_str("Writing ",x);
  if(!out_initialised)
    {
      mpfi_init(o_x);
      mpfr_init(o_l);
      mpz_init(z_l);
      mpz_init(z_r);
      out_initialised=true;
    }
  mpfi_diam_abs(o_l,x);
  if(mpfr_cmp_ui_2exp(o_l,1,-OP_ACC)>0)
    {
      printf("End points more than 2^{-%lu} apart. Exiting.\n",OP_ACC);
      exit(0);
    }
  mpfi_get_fr(o_l,x);
  double d=mpfi_get_d(x);
  if(mpfr_mul_2exp(o_l,o_l,OP_ACC,GMP_RNDN)!=0)
    {
      printf("mpfr_mul_2exp in out_real was inexact. Exiting.\n");
      exit(0);
    }
  mpfr_get_z(z_l,o_l,GMP_RNDN); // this is the OP_ACC result
  // let's double check it is OK
  mpz_add(z_r,z_l,z_l);
  mpz_sub_ui(z_r,z_r,1);
  mpfi_set_z(o_x,z_r);
  mpz_add_ui(z_r,z_r,2);
  mpfi_put_z(o_x,z_r);
  mpfi_mul_2si(o_x,o_x,-OP_ACC-1);
  mpfi_sub(o_x,x,o_x);
  if(!mpfi_contains_zero(o_x))
    {
      printf("Strange error in out_real. Exiting.\n");
      exit(0);
    }
  fwrite(&d,sizeof(double),1,outfile);
  //printf("outputing ");mpz_out_str(stdout,10,z_l);printf("\n");
  out_bytes(z_l,outfile);
}

void out_complex(mpfi_c_t z, FILE *outfile)
{
  out_real(z->re,outfile);
  out_real(z->im,outfile);
}

// setup to compute root numbers
void do_epsilons(mpfi_c_t *omegas, mpfi_c_t *in_vec, uint64_t q, DirichletGroup *G)
{
  mpfi_div_ui(itemp,mpfi_2_pi,q);
  for(uint64_t i=1;i<q;i++)
    if(G->is_coprime_to_q(i))
      {
	mpfi_mul_ui(itemp2,itemp,i);
	mpfi_sin(in_vec[i]->im,itemp2);
	mpfi_cos(in_vec[i]->re,itemp2);
      }
  G->DFTsum(omegas,in_vec); // now contains Gauss sums
  mpfi_set_ui(itemp,1);
  mpfi_div_ui(itemp2,itemp,q);
  mpfi_sqrt(itemp,itemp2); // q^{-1/2}
  for(uint64_t i=2;i<q;i++)
    if(G->character(i).is_primitive())
      {
	mpfi_c_mul_i(omegas[i],omegas[i],itemp); // div by sqrt(q)
	if(!G->character(i).is_even())
	  {
	    mpfi_swap(omegas[i]->re,omegas[i]->im); // 
	    mpfi_neg(omegas[i]->im,omegas[i]->im); // conj()
	  }
      }
}

int main(int argc, char** argv)
{
  if(argc!=3)
    {
      printf("Usage:- %s <q> <outfile>. Exiting.\n",argv[0]);
      exit(0);
    }
 
  uint64_t q=atol(argv[1]);
  if((q&3)==2)
    {
      printf("moduli congruent to 2 mod 4 have no primitive characters. Exiting.\n");
      exit(0);
    }

  FILE *outfile=fopen(argv[2],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[3]);
      exit(0);
    }

  if(fwrite(&q,sizeof(uint64_t),1,outfile)!=1)
    {
      printf("Error writing q. Exiting.\n");
      exit(0);
    }

  // initialise
  dft_init(PREC);
  mpfi_c_init(ctemp);
  mpfi_c_init(ctemp1);
  mpfi_init(itemp);
  mpfi_init(itemp2);
  mpfi_init(mpfi_2_pi);
  mpfi_const_pi(mpfi_2_pi);
  mpfi_mul_ui(mpfi_2_pi,mpfi_2_pi,2);

  mpfi_c_t *in_vec=(mpfi_c_t *)malloc(q*sizeof(mpfi_c_t));
  mpfi_c_t *out_vec=(mpfi_c_t *)malloc(q*sizeof(mpfi_c_t));
  mpfi_c_t *omegas=(mpfi_c_t *)malloc(q*sizeof(mpfi_c_t));
  mpfi_c_t *L1s=(mpfi_c_t *)malloc(q*sizeof(mpfi_c_t));

  if(!in_vec||!out_vec||!omegas||!L1s)
    {
      printf("Error allocating memory for in_vec/out_vec/omegas/L1s. Exiting.\n");
      exit(0);
    }

  for(uint64_t i=0;i<q;i++)
    {
      mpfi_c_init(in_vec[i]);
      mpfi_c_init(out_vec[i]);
      mpfi_c_init(omegas[i]);
      mpfi_c_init(L1s[i]);
    }


  DirichletGroup G(q);

  bool *is_prim,*is_ev;
  is_prim=(bool *)malloc(sizeof(bool)*q);
  is_ev=(bool *)malloc(sizeof(bool)*q);
  if(!is_prim||!is_ev)
    {
      printf("Error allocating memory for is_prim/is_ev. Exiting.\n");
      exit(0);
    }
  for(uint64_t i=1;i<q;i++)
    {
      is_prim[i]=G.character(i).is_primitive();
      is_ev[i]=G.character(i).is_even();
    }

  // set up the root numbers, used to compute Lambda from L
  do_epsilons(omegas,in_vec,q,&G);

  mpfi_t alpha;
  mpfi_init(alpha);
  
  for(uint64_t i=1;i<q;i++)
    if(G.is_coprime_to_q(i))
      {
	mpfi_set_ui(itemp,i);
	mpfi_div_ui(alpha,itemp,q);
	mpfi_hurwitz(in_vec[i]->re,0.5,alpha); // zeta(0.5,i/q)
	mpfi_set_ui(in_vec[i]->im,0);
      }

  G.DFTsum(out_vec,in_vec);

  // Now do L(1)
  // L(1)=q^{-1} sum -chi(a)psi(a/q)
  //
  for(uint64_t i=1;i<q;i++)
    if(G.is_coprime_to_q(i))
      {
	mpfi_set_ui(itemp,i);
	mpfi_div_ui(alpha,itemp,q);
	mpfi_psi(in_vec[i]->re,alpha);
	mpfi_neg(in_vec[i]->re,in_vec[i]->re); // -psi(i/q)
	mpfi_set_ui(in_vec[i]->im,0);
      }

  G.DFTsum(L1s,in_vec);

  mpfi_t q_minus_half;
  mpfi_init(q_minus_half);
  mpfi_set_ui(q_minus_half,q);
  mpfi_sqrt(alpha,q_minus_half);
  mpfi_inv(q_minus_half,alpha);

  for(uint64_t i=1;i<q;i++)
    if(is_prim[i])
      {
	fwrite(&i,sizeof(uint64_t),1,outfile);
	// multiply sum chi(a)zeta(1/2,a/q) by q^{-1/2}
	mpfi_mul(out_vec[i]->re,out_vec[i]->re,q_minus_half);
	mpfi_mul(out_vec[i]->im,out_vec[i]->im,q_minus_half);
	// multiply sum -chi(a)psi(a/q) by q^{-1}
	mpfi_div_ui(L1s[i]->re,L1s[i]->re,q);
	mpfi_div_ui(L1s[i]->im,L1s[i]->im,q);

	out_complex(omegas[i],outfile);
	out_complex(out_vec[i],outfile);
	out_complex(L1s[i],outfile);
      }
  printf("Successful completion on modulus %lu\n",q);
  return(0);
}

