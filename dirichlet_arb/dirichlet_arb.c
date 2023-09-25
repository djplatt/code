#include "acb_dirichlet.h"
#include "inttypes.h"
#include "stdbool.h"

#define Q ((uint64_t) 6323)
#define H1 ((double) 5.0)
#define STEP ((double) 0.0625)
#define POINTS ((uint64_t) H1/STEP+1)

uint64_t gcd (uint64_t a, uint64_t b)
/* Euclid algorithm gcd */
{
  uint64_t c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

int co_prime(uint64_t a, uint64_t b)
{
  return(gcd(a,b)==1);
};

// exp(2 pi i x)
void arb_e(acb_t res, arb_t x, uint64_t prec)
{
  static arb_t two_pi;
  static acb_t temp;
  static bool init=false;

  if(!init)
    {
      init=true;
      arb_init(two_pi);
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,1);
      acb_init(temp);
    }

  arb_mul(acb_imagref(temp),two_pi,x,prec);
  acb_exp(res,temp,prec);
}
  
void root_numbers(acb_t *res, dirichlet_group_t G, uint64_t q, uint64_t prec)
{
  acb_t *v;
  v=(acb_t *) malloc(sizeof(acb_t)*q);
  for(uint64_t i=0;i<q;i++)
    acb_init(v[i]);
  
  arb_t a;
  arb_init(a);

  arb_t sqrt_q;
  arb_init(sqrt_q);
  arb_sqrt_ui(sqrt_q,q,prec);
  
  acb_t temp;
  acb_init(temp);
  
  for(uint64_t i=1;i<q;i++)
    if(co_prime(i,q))
      {
	arb_set_ui(a,i);
	arb_div_ui(a,a,q,prec);
	arb_e(v[i],a,prec);
      }

  acb_dirichlet_dft((acb_ptr) res,(acb_srcptr) v,G,prec);
  // res now contains sum chi(a) e(a/q) = tau(chi)

  dirichlet_char_t ch;
  dirichlet_char_init(ch,G);
  dirichlet_char_first_primitive(ch,G);
  do
    {
      uint64_t i=dirichlet_char_exp(G,ch);
      acb_div_arb(temp,res[i],sqrt_q,prec);
      acb_inv(res[i],temp,prec);
      if(dirichlet_parity_char(G,ch)==1)
	acb_mul_onei(res[i],res[i]);
      // this is now the root number
      acb_sqrt(res[i],res[i],prec);
      if(!arb_is_positive(acb_realref(res[i])))
	acb_neg(res[i],res[i]);
      // take the square root with +ve real part.
    } while (dirichlet_char_next_primitive(ch,G)>=0);

  dirichlet_char_clear(ch);
  for(uint64_t i=0;i<q;i++)
    acb_clear(v[i]);
  arb_clear(a);
  arb_clear(sqrt_q);
  acb_clear(temp);
}

// compute L(s) for primitive chi into res
// will be called many times, hence use of static
// done via q^{-s} sum chi(a) zeta(s,a/q)
void arb_L(acb_t* res, acb_t s, dirichlet_group_t G, uint64_t q, uint64_t prec)
{
  static bool init=false;
  static arb_t log_q;
  static acb_t a,q_minus_s;
  static acb_t* temp;
  static dirichlet_char_t ch;
  if(!init)
    {
      init=true;
      temp=(acb_t *)malloc(sizeof(acb_t)*q);
      for(uint64_t i=0;i<q;i++)
	acb_init(temp[i]);
      acb_init(a);
      arb_init(log_q);
      acb_init(q_minus_s);
      arb_log_ui(log_q,q,prec);
      dirichlet_char_init(ch,G);
    }
  
  acb_mul_arb(q_minus_s,s,log_q,prec);
  acb_neg(q_minus_s,q_minus_s);
  acb_exp(q_minus_s,q_minus_s,prec);
  

  for(uint64_t i=1;i<q;i++)
    if(co_prime(i,q))
      {
	arb_set_ui(acb_realref(a),i);
	arb_div_ui(acb_realref(a),acb_realref(a),q,prec);
	acb_hurwitz_zeta(temp[i],s,a,prec);
      }

  acb_dirichlet_dft((acb_ptr) res,(acb_srcptr) temp,G,prec); // contains sum chi(a),zeta(s,a/q)
  
  dirichlet_char_first_primitive(ch,G);
  do
    {
      uint64_t i=dirichlet_char_exp(G,ch);
      acb_mul(res[i],res[i],q_minus_s,prec);
    } while (dirichlet_char_next_primitive(ch,G)>=0);
  
}

// compute Lambda(t)=epsilon (q/pi)^{it} * gamma(1/4+parity/2+it/2) * exp(pi t/4) * L(1/2+it)
void convert_L(acb_t *L, acb_t *roots, acb_t s, dirichlet_group_t G,
	       uint64_t q, uint64_t prec)
{
  static bool init=false;
  static arb_t log_q_pi;
  static acb_t s_by_two,s_one_by_two,gam_even,gam_odd,q_pi_it;
  static dirichlet_char_t ch;
  static arb_t pi_t_4,pi;
  if(!init)
    {
      init=true;
      arb_init(log_q_pi);
      acb_init(s_by_two);
      acb_init(s_one_by_two);
  
      acb_init(gam_even);
      acb_init(gam_odd);
      arb_init(log_q_pi);
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(pi_t_4);
      arb_div_ui(log_q_pi,pi,q,prec);
      arb_clear(pi);
      arb_log(log_q_pi,log_q_pi,prec);
      arb_neg(log_q_pi,log_q_pi); // log q/pi
      acb_init(q_pi_it);
      dirichlet_char_init(ch,G);
    }

  arb_mul(pi_t_4,pi,acb_imagref(s),prec);
  arb_mul_2exp_si(pi_t_4,pi_t_4,-2);
  arb_exp(pi_t_4,pi_t_4,prec);
  
  acb_add_ui(s_one_by_two,s,1,prec);
  acb_mul_2exp_si(s_one_by_two,s_one_by_two,-1); // (s+1)/2
  acb_mul_2exp_si(s_by_two,s,-1); // s/2
  acb_gamma(gam_even,s_by_two,prec); // gamma(s/2)
  acb_gamma(gam_odd,s_one_by_two,prec); // gamma ((s+1)/2)
  arb_mul(acb_imagref(q_pi_it),acb_imagref(s),log_q_pi,prec);
  arb_mul_2exp_si(acb_imagref(q_pi_it),acb_imagref(q_pi_it),-1);
  arb_zero(acb_realref(q_pi_it));
  acb_exp(q_pi_it,q_pi_it,prec); // q/pi^it/2
  
  acb_mul(gam_even,gam_even,q_pi_it,prec);
  acb_mul_arb(gam_even,gam_even,pi_t_4,prec);
  acb_mul(gam_odd,gam_odd,q_pi_it,prec);
  acb_mul_arb(gam_odd,gam_odd,pi_t_4,prec);
  dirichlet_char_first_primitive(ch,G);
  do
    {
      uint64_t i=dirichlet_char_exp(G,ch);
      acb_mul(L[i],L[i],roots[i],prec);
      if(dirichlet_parity_char(G,ch)==0)
	acb_mul(L[i],L[i],gam_even,prec);
      else
	acb_mul(L[i],L[i],gam_odd,prec);
      if(!arb_contains_zero(acb_imagref(L[i])))
	{
	  printf("Fatal error with q = %lu char = %lu. Non real Lambda.\n",
		 q,i);
	  printf("s=");acb_printd(s,10);
	  printf("\ngam_even=");acb_printd(gam_even,20);
	  printf("\ngam_odd=");acb_printd(gam_odd,20);
	  exit(0);
	}
    } while (dirichlet_char_next_primitive(ch,G)>=0);
}

void copy_L(arb_t **lambdas,acb_t *L,dirichlet_group_t G,uint64_t ptr)
{
  static bool init=false;
  static dirichlet_char_t ch;
  if(!init)
    {
      init=true;
      dirichlet_char_init(ch,G);
    }
  dirichlet_char_first_primitive(ch,G);
  do
    {
      uint64_t i=dirichlet_char_exp(G,ch);
      arb_set(lambdas[i][ptr],acb_realref(L[i]));
    } while (dirichlet_char_next_primitive(ch,G)>=0);
}


int main()
{
  uint64_t q=Q;
  uint64_t prec=200;
  acb_t* roots;
  roots=(acb_t*) malloc(q*sizeof(acb_t));
  acb_t* L;
  L=(acb_t*) malloc(q*sizeof(acb_t));
  
  for(uint64_t i=0;i<q;i++)
    {
      acb_init(roots[i]);
      acb_init(L[i]);
    }
  
  dirichlet_group_t G;

  dirichlet_group_init(G,q);
  uint64_t num_prims=dirichlet_group_num_primitive(G);

  arb_t **lambdas;
  lambdas=(arb_t **)malloc(sizeof(arb_t *)*q);
  for(uint64_t i=0;i<q;i++)
    {
      lambdas[i]=(arb_t *)malloc(sizeof(arb_t)*POINTS);
      for(uint64_t j=0;j<POINTS;j++)
	arb_init(lambdas[i][j]);
    }

  
  root_numbers(roots,G,q,prec);

  acb_t s;
  acb_init(s);
  arb_set_d(acb_realref(s),0.5);
  
  uint64_t ptr=0;
  for(double im_s=0.0;im_s<=H1;im_s+=STEP,ptr++)
    {    
      arb_set_d(acb_imagref(s),im_s);
      arb_L(L,s,G,q,prec); // L contains L(chi,s)
      convert_L(L,roots,s,G,q,prec); // L contains Lambda(chi,Im(s))
      copy_L(lambdas,L,G,ptr); // lambda[i][j] contains Lambda(chi_i,j*STEP)
      // check copy_L is working
    }

  return 0;
}
