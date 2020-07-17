//
// check_zeros.cpp
//
// computes \Lambda_\chi(1/2_it)=\omega (q/\pi)^{it}\Gamma((1/2+it+a_\chi)/2)\exp(Pi t/4)L_\chi(1/2+it)
// where a_chi=(1-\chi(-1))/2
// omega=
//
#include <iostream>
#include <cstdlib>
#include "characters.h"
#include "mpfi_hurwitz.h"

#define OP_ACC (101)

using namespace std;

uint64_t PREC,num_checked=0;

mpfr_t a[4],b[4],c[4];

void mpfi_c_hurwitz(mpfi_c_ptr res, mpfi_c_ptr z, mpfi_c_ptr alpha)
{
  mpfi_get_left(b[0],z->re);
  mpfi_get_right(b[1],z->re);
  mpfi_get_left(b[2],z->im);
  mpfi_get_right(b[3],z->im);
  mpfi_get_left(c[0],alpha->re);
  mpfi_get_right(c[1],alpha->re);
  mpfi_get_left(c[2],alpha->im);
  mpfi_get_right(c[3],alpha->im);
  mpfi_c_hurwitz1(a,b,c,PREC);
  mpfi_set_fr(res->re,a[0]);
  mpfi_put_fr(res->re,a[1]);
  mpfi_set_fr(res->im,a[2]);
  mpfi_put_fr(res->im,a[3]);
  //mpfi_c_print_str("z=",z);
  //mpfi_c_print_str("alpha=",alpha);
  //mpfi_c_print_str("zeta=",res);
}
void mpfi_c_lngamma(mpfi_c_ptr res, mpfi_c_ptr z)
{
  mpfi_get_left(b[0],z->re);
  mpfi_get_right(b[1],z->re);
  mpfi_get_left(b[2],z->im);
  mpfi_get_right(b[3],z->im);
  mpfi_c_lngamma1(a,b,PREC);
  mpfi_set_fr(res->re,a[0]);
  mpfi_put_fr(res->re,a[1]);
  mpfi_set_fr(res->im,a[2]);
  mpfi_put_fr(res->im,a[3]);
  //mpfi_c_print_str("z=",z);
  //mpfi_c_print_str("lng=",res);
}


void read_null(FILE* infile)
{
  uint8_t buff[13];
  if(fread(buff,sizeof(uint8_t),13,infile)!=13)
    {
      printf("Fatal error reading null zero. Exiting.\n");
      exit(0);
    }
  uint64_t i;
  for(i=0;i<13;i++)
    if(buff[i]!=0)
      {
	printf("Zero was not null. Exiting.\n");
	exit(0);
      }
}

mpfi_t pm1;
bool in_bytes_initialised=false;

#define mpfi_is_exact(fn) (MPFI_FLAGS_BOTH_ENDPOINTS_EXACT==(fn))

void in_bytes(mpfi_ptr t, FILE *infile)
{
  if(!in_bytes_initialised)
    {
      mpfi_init(pm1);
      mpfi_set_ui(pm1,1);
      mpfi_div_2ui(pm1,pm1,OP_ACC+1);
      //mpfi_print_str("Acc=+/-",pm1);
      in_bytes_initialised=true;
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
  mpfi_set_ui(t,c);
  mpfi_mul_2ui(t,t,32);
  mpfi_add_ui(t,t,b);
  mpfi_mul_2ui(t,t,64);
  if(!mpfi_is_exact(mpfi_add_ui(t,t,a)))
    {
      printf("Error reading zeros. Endpoints not exact. Increase precison. Exiting.\n");
      exit(0);
    }
  mpfi_div_2ui(t,t,OP_ACC);
}


// read the next imaginary part of rho into del_t
inline void next_rho(mpfi_ptr del_t, FILE *infile)
{
  in_bytes(del_t,infile);
}

mpfi_c_t halpha,hhur,hchi;
bool sum_hur_initialised=false;

void sum_hur(mpfi_c_ptr res, mpfi_c_ptr s, uint64_t q, uint64_t index, DirichletGroup &G)
{
  if(!sum_hur_initialised)
    {
      mpfi_c_init(halpha);
      mpfi_set_ui(halpha->im,0);
      mpfi_c_init(hhur);
      mpfi_c_init(hchi);
      sum_hur_initialised=true;
    }
  mpfi_set_ui(res->re,0);
  mpfi_set_ui(res->im,0);
  for(uint64_t a=1;a<q;a++)
    {
      if(GCD(a,q)>1)
	continue;
      mpfi_set_ui(halpha->re,a);
      mpfi_div_ui(halpha->re,halpha->re,q);
      mpfi_c_hurwitz(hhur,s,halpha);
      G.chi(hchi,index,a);
      mpfi_c_mul(hhur,hhur,hchi);
      mpfi_c_add(res,res,hhur);
    }
}
#define POS (0)
#define NEG (1)
#define sign_t uint8_t

mpfi_c_t Lres,Lomega,lnt,lng,Lqs,Ls,Lg;
mpfi_t Lq,Larg,mpfi_pi;
bool L_init=false;
void L(mpfi_ptr t, uint64_t q, uint64_t index, DirichletGroup &G, bool real_p)
{
  if(!L_init)
    {
      mpfi_c_init(Lres);
      mpfi_c_init(Lomega);
      mpfi_c_init(lnt);
      mpfi_c_init(lng);
      mpfi_c_init(Lqs);
      mpfi_c_init(Ls);
      mpfi_c_init(Lg);
      mpfi_init(Lq);
      mpfi_init(Larg);
      mpfi_init(mpfi_pi);
      mpfi_const_pi(mpfi_pi);
      L_init=true;
    }
  bool even_p=G.character(index).is_even();

  mpfi_set_d(Ls->re,0.5);
  if(!mpfi_is_exact(mpfi_sub(Ls->im,t,pm1)))
    {
      printf("Endpoints not exact. Increase precision. Exiting.\n");
      exit(0);
    } // Ls=s0=1/2+i(t-delta)
  //printf("s=");mpfi_c_printn(Ls,150);
  mpfi_c_set(lnt,Ls);
  if(!even_p)
    mpfi_add_ui(lnt->re,lnt->re,1);
  mpfi_c_div_ui(lnt,lnt,2); // lnt = (s0+a)/2
  mpfi_c_lngamma(lng,lnt); // log gamma (s0+a)/2
  mpfi_c_exp(Lg,lng); // Lg=gamma((s0+a)/2)
  mpfi_c_div_d(lnt,Ls,-2.0); // lnt=-s0/2
  mpfi_set_ui(Lq,q);
  mpfi_c_pow_i_c(Lqs,Lq,Ls); // Lqs=q^s0
  //mpfi_c_print_str("q^(s)=",Lqs);
  sum_hur(Lres,Ls,q,index,G); // Lres=sum chi(a) zeta(s0,a/q)
  //printf("sum chi(a)*zeta(s,a/q)=");
  //mpfi_c_printn(Lres,150);
  mpfi_c_div(Lres,Lres,Lqs); // Lres=L(s0) 
  //printf("L_%lu_%lu(",q,index);
  //mpfi_c_print(Ls);
  //mpfi_c_print_str(")=",Lres);
  mpfi_c_mul(Lres,Lres,Lg); // Lres=Gam((s0+a)/2)*L(s0) 
  mpfi_div_ui(Lq,mpfi_pi,q);
  mpfi_c_pow_i_c(Lqs,Lq,lnt); // Lqs=(Pi/q)^(-s0/2)
  mpfi_c_mul(Lres,Lres,Lqs); // 
  if(real_p)
    {
      mpfi_c_set_ui(Lomega,1,0);
      mpfi_c_set(lnt,Lres);
    }
  else
    {
      mpfi_c_arg(Larg,Lres);
      //mpfi_print_str("arg1=",Larg);
      mpfi_sin(Lomega->re,Larg);
      mpfi_neg(Lomega->im,Lomega->re);
      mpfi_cos(Lomega->re,Larg);
      mpfi_c_mul(lnt,Lres,Lomega);
    }
  //printf("Lam_%lu_%lu(",q,index);
  //mpfi_print(Ls->im);
  //mpfi_c_print_str(")=",lnt);
  
  if(!mpfi_contains_zero(lnt->im))
    {
      printf("Error, Lambda does not contain zero. Exiting.\n");
      exit(0);
    }

  sign_t sign;
  if(mpfi_is_pos(lnt->re))
    sign=POS;
  else
    {
      if(mpfi_is_neg(lnt->re))
	sign=NEG;
      else
	{
	  printf("Error, Lambad is of indeterminate sign. Exiting.\n");
	  exit(0);
	}
    }

  mpfi_add(Ls->im,t,pm1);
  mpfi_c_set(lnt,Ls);
  if(!even_p)
    mpfi_add_ui(lnt->re,lnt->re,1);
  mpfi_c_div_ui(lnt,lnt,2);
  mpfi_c_lngamma(lng,lnt); // log gamma (s+a)/2
  mpfi_c_exp(Lg,lng);
  mpfi_c_div_d(lnt,Ls,-2.0); // -s/2
  mpfi_set_ui(Lq,q);
  mpfi_c_pow_i_c(Lqs,Lq,Ls);
  sum_hur(Lres,Ls,q,index,G); // sum chi(a) zeta(s,a/q)
  mpfi_c_div(Lres,Lres,Lqs);
  mpfi_c_mul(Lres,Lres,Lg); // Gam((s+a)/2) sum chi(a) zeta(s,a/q)
  mpfi_div_ui(Lq,mpfi_pi,q);
  mpfi_c_pow_i_c(Lqs,Lq,lnt);
  mpfi_c_mul(Lres,Lres,Lqs);
  mpfi_c_mul(lnt,Lres,Lomega);
  //printf("Lam_%lu_%lu(",q,index);
  //mpfi_print(Ls->im);
  //mpfi_c_print_str(")=",lnt);
 
  if(!mpfi_contains_zero(lnt->im))
    {
      printf("Error, Lambda does not contain zero. Exiting.\n");
      exit(0);
    }

  if(mpfi_is_pos(lnt->re))
    {
      if(sign==POS)
	{
	  printf("Error, both lambdas were positive. Exiting.\n");
	  exit(0);
	}
      else
	return;
    }
  if(mpfi_is_neg(lnt->re))
    {
      if(sign==NEG)
	{
	  printf("Error, both lambdas were negative. Exiting.\n");
	  exit(0);
	}
      else
	return;
    }

}


mpfi_t ddel_t,dt;
bool dz_init=false;
void do_zeros(uint64_t num_zeros, uint64_t q, FILE *infile, double prob, uint64_t index, DirichletGroup &G)
{
  if(!dz_init)
    {
      mpfi_init(ddel_t);
      mpfi_init(dt);
      dz_init=true;
    }
  printf("In do_zeros with index %lu",index);
  uint64_t z;
  uint8_t sign1,sign2;
  if(q==3)
    mpfi_set_ui(dt,8);
  else
    mpfi_set_ui(dt,0);
  for(z=0;z<num_zeros;z++)
    {
      next_rho(ddel_t,infile);
      if(!mpfi_is_exact(mpfi_add(dt,dt,ddel_t)))
	{
	  printf("Endpoints not exact. Increase precision. Exiting.\n");
	  exit(0);
	}
      if(((double) rand()/ (double) RAND_MAX)<=prob)
	{
	  num_checked++;
	  //printf("Checking zero at ");
	  //mpfi_printn(dt,120);
	  printf(".");
	  fflush(stdout);
	  L(dt,q,index,G,InvMod(index,q)==index);
	}
    }
  printf("\n");
}

void test_chis(uint64_t q, uint64_t index, DirichletGroup &G)
{
  mpfi_c_t chi,chi_sum;
  mpfi_c_init(chi);
  mpfi_c_init(chi_sum);
  mpfi_c_set_ui(chi_sum,0,0);
  for(uint64_t a=1;a<q;a++)
    if(GCD(a,q)==1)
      {
	G.chi(chi,index,a);
	mpfi_c_add(chi_sum,chi_sum,chi);
      }
  mpfi_c_print_str("Sum=",chi_sum);
  mpfi_c_clear(chi);
  mpfi_c_clear(chi_sum);
}

int main(int argc, char **argv)
{
  int i;
  printf("Command line:- ");
  for(i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=4)
    {
      printf("Fatal error in main: usage %s <prec> <zeros file> <prob>. Exiting.\n",argv[0]);
      exit(0);
    }
  PREC=atol(argv[1]);
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }

  double prob=atof(argv[3]);

  dft_init(PREC);
  for(uint64_t i=0;i<4;i++)
    {
      mpfr_init(a[i]);
      mpfr_init(b[i]);
      mpfr_init(c[i]);
    }

  uint64_t q,index,num_zeros,index1,num_zeros1,num_prims=0;
  if(fread(&q,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading q from %s. Exiting.\n",argv[2]);
      exit(0);
    }
  //printf("starting G setup\n");

  DirichletGroup G(q);


  //printf("G setup\n");

  while(fread(&index,sizeof(uint64_t),1,infile)==1)
    {
      //test_chis(q,index,G);
      if(fread(&num_zeros,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading num_zeros from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      do_zeros(num_zeros,q,infile,prob,index,G);
      uint64_t conj_index=InvMod(index,q);
      if(conj_index==index) // a real character
	continue;
      if(fread(&index1,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading index1 from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      if(index1!=conj_index)
	{
	  printf("Conj index read (%lu) does not match expected (%lu). Exiting.\n",index1,conj_index);
	  exit(0);
	}
      if(fread(&num_zeros1,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading num_zeros from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      do_zeros(num_zeros1,q,infile,prob,index1,G);
    }
  printf("We checked %lu zeros.\n",num_checked);
  return(0);
}
