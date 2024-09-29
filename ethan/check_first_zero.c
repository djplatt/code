#include "stdlib.h"
#include "flint/acb_dirichlet.h"
#include "inttypes.h"
#include "stdbool.h"
#define OP_ACC ((int64_t) 101)
#define ZERO_LEN ((uint64_t) 13)


void XGCD(int64_t *d, 
	  int64_t *s, 
	  int64_t *t, 
	  int64_t a, 
	  int64_t b)
{
   int64_t  u, v, u0, v0, u1, v1, u2, v2, q, r;

   int64_t aneg = 0, bneg = 0;

   if (a < 0) {
      a = -a;
      aneg = 1;
   }

   if (b < 0) {
      b = -b;
      bneg = 1;
   }

   u1=1; v1=0;
   u2=0; v2=1;
   u = a; v = b;

   while (v != 0) {
      q = u / v;
      r = u % v;
      u = v;
      v = r;
      u0 = u2;
      v0 = v2;
      u2 =  u1 - q*u2;
      v2 = v1- q*v2;
      u1 = u0;
      v1 = v0;
   }

   if (aneg)
      u1 = -u1;

   if (bneg)
      v1 = -v1;

   d[0] = u;
   s[0] = u1;
   t[0] = v1;
}
   
// taken from NTL:
int64_t InvMod(int64_t a, int64_t n)
{
   int64_t d, s, t;

   XGCD(&d, &s, &t, a, n);
   if (d != 1) return -1;
   if (s < 0)
      return s + n;
   else
      return s;
}


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
  
void in_bytes(arb_ptr t, FILE *infile, uint64_t prec)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;

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
  //printf("c=%u\n",c);
  arb_set_ui(t,c);
  arb_mul_2exp_si(t,t,32);
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64);
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-OP_ACC);
}

void do_rho(arb_t t, dirichlet_group_t G, uint64_t index, uint64_t q, uint64_t prec)
{
  static bool init=false;
  static dirichlet_char_t ch;
  static arb_t pm1,log_q_pi;
  static acb_t s0,s1,res0,res1,temp0,temp00,temp000,temp1,temp11,temp111,root;
  if(!init)
    {
      init=true;
      dirichlet_char_init(ch,G);
      arb_init(pm1);
      arb_set_ui(pm1,1);
      arb_mul_2exp_si(pm1,pm1,-OP_ACC-1);
      acb_init(s0);
      acb_init(s1);
      acb_init(res0);
      acb_init(res1);
      acb_init(temp0);
      acb_init(temp00);
      acb_init(temp000);
      acb_init(temp1);
      acb_init(temp11);
      acb_init(temp111);
      acb_init(root);
      arb_set_d(acb_realref(s0),0.5);
      arb_set_d(acb_realref(s1),0.5);
      arb_init(log_q_pi);
      arb_const_pi(log_q_pi,prec);
      arb_div_ui(log_q_pi,log_q_pi,q,prec);
      arb_log(log_q_pi,log_q_pi,prec);
      arb_neg(log_q_pi,log_q_pi);
      arb_mul_2exp_si(log_q_pi,log_q_pi,-1);
    }
  //printf("first zero for index %lu is at ",index);arb_printd(t,40);printf("\n");
  dirichlet_char_log(ch,G,index);
  //acb_set(root,roots[index]);
  acb_dirichlet_root_number_theta(root,G,ch,prec);
  arb_neg(acb_imagref(root),acb_imagref(root)); // arbitary, but works
  acb_sqrt(root,root,prec);
  //printf("sqrt(root) = ");acb_printd(root,20);printf("\n");
  arb_set(acb_imagref(s0),t);
  arb_set(acb_imagref(s1),t);
  arb_sub(acb_imagref(s0),acb_imagref(s0),pm1,prec);
  arb_add(acb_imagref(s1),acb_imagref(s1),pm1,prec);
  acb_dirichlet_l(res0,s0,G,ch,prec);
  acb_dirichlet_l(res1,s1,G,ch,prec);
  arb_mul(acb_imagref(temp0),acb_imagref(s0),log_q_pi,prec);
  arb_mul(acb_imagref(temp1),acb_imagref(s1),log_q_pi,prec);
  acb_exp(temp00,temp0,prec); // q/pi^it/2
  acb_exp(temp11,temp1,prec);
  //printf("q/pi^(it/2) = ");acb_printd(temp00,20);printf("\n");
  if(dirichlet_parity_char(G,ch)==1)
    {
      //printf("character %lu is odd.\n",index);
      acb_add_ui(temp000,s0,1,prec);
      acb_add_ui(temp111,s1,1,prec);
    }
  else
    {
      acb_set(temp000,s0);
      acb_set(temp111,s1);
    }
  //printf("s+a = ");acb_printd(temp000,20);printf("\n");
  acb_mul_2exp_si(temp000,temp000,-1);
  acb_mul_2exp_si(temp111,temp111,-1);
  acb_gamma(temp000,temp000,prec);
  acb_gamma(temp111,temp111,prec);
  //printf("gamma((s+a)/2) = ");acb_printd(temp000,20);printf("\n");
  acb_mul(res0,res0,root,prec); // * sqrt(omega)
  acb_mul(res0,res0,temp00,prec); // * (q/Pi)^(it/2)
  acb_mul(res0,res0,temp000,prec); // * gamma((s+a)/2)
  acb_mul(res1,res1,root,prec);
  acb_mul(res1,res1,temp11,prec);
  acb_mul(res1,res1,temp111,prec);
  if(!arb_contains_zero(acb_imagref(res0)))
    {
      printf("Error Lambda not real index %lu\n",index);
      acb_printd(res1,20);printf("\n");
      exit(0);
    }
  if(!arb_contains_zero(acb_imagref(res1)))
    {
      printf("Error Lambda not real index %lu\n",index);
      acb_printd(res1,20);printf("\n");
      exit(0);
    }
  arb_mul(acb_realref(res0),acb_realref(res0),acb_realref(res1),prec);
  if(!arb_is_negative(acb_realref(res0)))
    printf("Failed to bracket zero at q = %lu index = %lu\n",q,index);
}

int main(int argc, char** argv)
{
  printf("Command line:- ");
  for(uint64_t i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Fatal error in main: usage %s <prec> <zeros file>. Exiting.\n",argv[0]);
      exit(0);
    }
  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Fatal error in main: failed to open zeros file for binary input. Exiting.\n");
      exit(0);
    }
  uint64_t prec=atol(argv[1]);
  
  uint64_t q,index,num_zeros,index1,num_zeros1;
  if(fread(&q,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading q from %s. Exiting.\n",argv[2]);
      exit(0);
    }

  
  dirichlet_group_t G;

  dirichlet_group_init(G,q);

  arb_t t,del_t;
  arb_init(t);
  arb_init(del_t);

    while(fread(&index,sizeof(uint64_t),1,infile)==1)
    {
      if(fread(&num_zeros,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading num_zeros from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      uint64_t z;
      if(q==3)
	arb_set_ui(t,8);
      else
	arb_set_ui(t,0);
      in_bytes(del_t,infile,prec);
      arb_add(t,t,del_t,prec); // exact
      do_rho(t,G,index,q,prec); // check
      if(fseek(infile,(num_zeros-1)*ZERO_LEN,SEEK_CUR)!=0)
	{
	  printf("Fatal error skipping zeros %lu %lu. Exiting.\n",q,index);
	  exit(0);
	}
      uint64_t conj_index=InvMod(index,q);
      if(conj_index==index) // a real character
	  continue;
      if(fread(&index1,sizeof(uint64_t),1,infile)!=1)
	{
	  printf("Error reading index1 from %s. Exiting.\n",argv[2]);
	  exit(0);
	}
      //printf("Conjugate index read was %lu\n",index1);
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
      //printf("Doing %lu conjugate zeros.\n",num_zeros1);
      arb_set_ui(t,0);
      in_bytes(del_t,infile,prec);
      arb_add(t,t,del_t,prec); // exact
      do_rho(t,G,index1,q,prec); // check
      if(fseek(infile,(num_zeros1-1)*ZERO_LEN,SEEK_CUR)!=0)
	{
	  printf("Fatal error skipping zeros %lu %lu. Exiting.\n",q,index1);
	  exit(0);
	}
      
    }

  return 0;
}
