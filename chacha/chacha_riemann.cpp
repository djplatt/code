#include "inttypes.h"
#include "dirichlet.h"
#include "acb_dirichlet.h"
#include "quad.h"

#define GOOD_GO_LIMIT (128)

dirichlet_group_t dg;
dirichlet_char_t chi;
arb_t six,half;

// horizontal contour
void do_int1(arb_t res, acb_t z, arb_t dt, arb_t dt1,int64_t prec)
{
  static bool init=false;
  static acb_ptr tmp;
  static acb_t z1;
  if(!init)
    {
      init=true;
      tmp=(acb_ptr)malloc(sizeof(acb_t)*3);
      for(uint64_t n=0;n<3;n++)
	acb_init(tmp+n);
      acb_init(z1);
    }
  acb_dirichlet_l_jet(tmp,z,dg,chi,0,3,prec);
  acb_mul_2exp_si(tmp+2,tmp+2,1);
  acb_div(tmp,tmp+2,tmp+1,prec);
  arb_mul(res,acb_imagref(tmp),dt1,prec);
}

// vertical contour
void do_int2(arb_t res, acb_t z, arb_t dt1,int64_t prec)
{
  static bool init=false;
  static acb_ptr tmp;
  static acb_t z1;
  if(!init)
    {
      init=true;
      tmp=(acb_ptr)malloc(sizeof(acb_t)*3);
      for(uint64_t n=0;n<3;n++)
	acb_init(tmp+n);
      acb_init(z1);
    }
  arb_get_mid_arb(acb_realref(z1),acb_realref(z));
  arb_set(acb_imagref(z1),acb_imagref(z1));

  acb_dirichlet_l_jet(tmp,z,dg,chi,0,3,prec);
  acb_mul_2exp_si(tmp+2,tmp+2,1);
  acb_div(tmp,tmp+2,tmp+1,prec);
  arb_mul(res,acb_realref(tmp),dt1,prec);
}


int main(int argc, char **argv)
{

  if(argc!=4)
    {
      printf("Usage:- %s <q> <prec> <-log2 step>\n",argv[0]);
      return 0;
    }

  uint64_t q=atol(argv[1]);
  if((q%4)==2)
    {
      printf("There are no primitive characters for q = 2 mod 4.\n");
      return 0;
    }
  int64_t prec=atol(argv[2]);
  int64_t l2step=-atol(argv[3]);

  fmpz_t nz;
  fmpz_init(nz);

  arb_t two_pi;
  arb_init(two_pi);
  arb_const_pi(two_pi,prec);
  arb_mul_2exp_si(two_pi,two_pi,1);

  dirichlet_group_init(dg,q);

  dirichlet_char_init(chi,dg);
  dirichlet_char_first_primitive(chi,dg);

  acb_t z; // piece of the contour
  arb_t dt,dt1; // delta t for bottom and right
  arb_init(dt);arb_init(dt1);
  acb_init(z);

  arb_t tmp;
  arb_init(half);arb_init(six);arb_init(tmp);
  arb_set_d(half,0.5);

  arb_t rest,resb,resl,resr;
  arb_init(resb);
  arb_init(rest);
  arb_init(resl);
  arb_init(resr);

  switch(q)
    {
    case 3: arb_set_d(six,6.0);break;
    case 4: arb_set_d(six,5.0);break;
    case 5: 
    case 7: arb_set_d(six,4.0);break;
    default: arb_set_d(six,3.0);
    }

  uint64_t cn;

  while(true) // iterate over primitive characters
    {
      cn=dirichlet_char_exp(dg,chi);
      printf("Character number %lu\n",cn);

      // setup bottom contour = 0->0.5
      arb_zero(resb);
      arb_set_ui(dt,1);
      arb_mul_2exp_si(dt,dt,l2step); // 2^-l2step
      arb_mul_2exp_si(dt1,dt,1);
      arb_set(acb_realref(z),dt);
      arb_set_d(acb_imagref(z),0.0);
      arb_add_error(acb_realref(z),dt);
      while(true)
	{
	  //printf("z = ");acb_printd(z,20);printf("\n");
	  if(arb_ge(acb_realref(z),half)) break;
	  do_int1(tmp,z,dt,dt1,prec);
	  arb_add(resb,resb,tmp,prec);
	  arb_add(acb_realref(z),acb_realref(z),dt1,prec);
	}
      printf("Int over bottom = ");arb_printd(resb,20);printf("\n");


      // setup top contour = 0.5+6i->0+6i
      arb_zero(rest);
      arb_set_ui(dt,1);
      arb_mul_2exp_si(dt,dt,l2step); // 2^-l2step
      arb_mul_2exp_si(dt1,dt,1);
      arb_set(acb_realref(z),dt);
      arb_set(acb_imagref(z),six);
      arb_add_error(acb_realref(z),dt);
      while(true)
	{
	  //printf("z = ");acb_printd(z,20);printf("\n");
	  if(arb_ge(acb_realref(z),half)) break;
	  do_int1(tmp,z,dt,dt1,prec);
	  arb_sub(rest,rest,tmp,prec); // subtract as contour is r->l
	  arb_add(acb_realref(z),acb_realref(z),dt1,prec);
	}
      
      printf("Int over top = ");arb_printd(rest,20);printf("\n");

      // setup right contour = 0.5->0.5+6i
      arb_zero(resr);
      arb_set_ui(dt,1);
      arb_mul_2exp_si(dt,dt,l2step); // 2^-l2step
      arb_mul_2exp_si(dt1,dt,1);
      arb_set(acb_realref(z),half);
      arb_set(acb_imagref(z),dt);
      arb_add_error(acb_imagref(z),dt);
      while(true)
	{
	  //printf("z = ");acb_printd(z,20);printf("\n");
	  if(arb_ge(acb_imagref(z),six)) break;
	  do_int2(tmp,z,dt1,prec);
	  arb_add(resr,resr,tmp,prec);
	  arb_add(acb_imagref(z),acb_imagref(z),dt1,prec);
	}
      
      printf("Int over right = ");arb_printd(resr,20);printf("\n");

      // setup left contour = 6i->0
      arb_zero(resl);
      arb_set_ui(dt,1);
      arb_mul_2exp_si(dt,dt,l2step); // 2^-l2step
      arb_mul_2exp_si(dt1,dt,1);
      arb_zero(acb_realref(z));
      arb_set(acb_imagref(z),dt);
      arb_add_error(acb_imagref(z),dt);
      while(true)
	{
	  //printf("z = ");acb_printd(z,20);printf("\n");
	  if(arb_ge(acb_imagref(z),six)) break;
	  do_int2(tmp,z,dt1,prec);
	  arb_sub(resl,resl,tmp,prec);
	  arb_add(acb_imagref(z),acb_imagref(z),dt1,prec);
	}
      
      printf("Int over left = ");arb_printd(resl,20);printf("\n");


      arb_add(tmp,resl,resr,prec);
      arb_add(resl,rest,resb,prec);
      arb_add(rest,tmp,resl,prec);
      arb_div(resl,rest,two_pi,prec);
      printf("Re I = ");arb_printd(resl,20);printf("\n");
      if(arb_get_unique_fmpz(nz,resl))
	{
	  printf("result for q = %lu chi = %lu is N-P = ",q,cn);
	  fmpz_print(nz);
	}
      else
	{
	  arb_printd(resl,20);
	  printf(" does not bracket a unique integer.");
	}
      printf("\n");

      if(dirichlet_char_next_primitive(chi,dg)<0)
	break;
    }
  return 0;
}
