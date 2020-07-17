#include "inttypes.h"
#include "dirichlet.h"
#include "acb_dirichlet.h"
#include "quad.h"

#define max(a,b) (a>b ? a : b)

dirichlet_group_t dg;
dirichlet_char_t chi;
arb_t six,half;

// bottom contour from 0 to half
void arb_fb(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_ptr tmp;
  static acb_t z;
  if(!init)
    {
      init=true;
      tmp=(acb_ptr)malloc(sizeof(acb_t)*3);
      for(uint64_t n=0;n<3;n++)
	acb_init(tmp+n);
      acb_init(z);
    }
  arb_set(acb_realref(z),t);
  arb_zero(acb_imagref(z));
  acb_dirichlet_l_jet(tmp,z,dg,chi,0,3,prec);
  acb_mul_2exp_si(tmp+2,tmp+2,1);
  acb_div(tmp,tmp+2,tmp+1,prec);
  arb_set(res,acb_imagref(tmp));
  //printf("arb_fb at ");acb_printd(z,10);printf(" returning ");arb_printd(res,30);printf("\n");
}

// top contour from six i to half + six i
void arb_ft(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_ptr tmp;
  static acb_t z;
  if(!init)
    {
      init=true;
      tmp=(acb_ptr)malloc(sizeof(acb_t)*3);
      for(uint64_t n=0;n<3;n++)
	acb_init(tmp+n);
      acb_init(z);
    }
  arb_set(acb_realref(z),t);
  arb_set(acb_imagref(z),six);
  acb_dirichlet_l_jet(tmp,z,dg,chi,0,3,prec);
  acb_mul_2exp_si(tmp+2,tmp+2,1);
  acb_div(tmp,tmp+2,tmp+1,prec);
  arb_neg(res,acb_imagref(tmp));
  //printf("arb_ft at ");acb_printd(z,10);printf(" returning ");arb_printd(res,30);printf("\n");
}

// right from half to half + six i
void arb_fr(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_ptr tmp;
  static acb_t z;
  if(!init)
    {
      init=true;
      tmp=(acb_ptr)malloc(sizeof(acb_t)*3);
      for(uint64_t n=0;n<3;n++)
	acb_init(tmp+n);
      acb_init(z);
    }
  arb_set(acb_imagref(z),t);
  arb_set(acb_realref(z),half);
  acb_dirichlet_l_jet(tmp,z,dg,chi,0,3,prec);
  acb_mul_2exp_si(tmp+2,tmp+2,1);
  acb_div(tmp,tmp+2,tmp+1,prec);
  arb_set(res,acb_realref(tmp));
}

// left from 0 to six i
void arb_fl(arb_t res, const arb_t t, int64_t prec)
{
  static bool init=false;
  static acb_ptr tmp;
  static acb_t z;
  if(!init)
    {
      init=true;
      tmp=(acb_ptr)malloc(sizeof(acb_t)*3);
      for(uint64_t n=0;n<3;n++)
	acb_init(tmp+n);
      acb_init(z);
    }
  arb_set(acb_imagref(z),t);
  arb_zero(acb_realref(z));
  acb_dirichlet_l_jet(tmp,z,dg,chi,0,3,prec);
  acb_mul_2exp_si(tmp+2,tmp+2,1);
  acb_div(tmp,tmp+2,tmp+1,prec);
  arb_neg(res,acb_realref(tmp));
}

arb_t maxd;

void do_int1(arb_t res, acb_t z, arb_t dt1,
	     void (*f1)(arb_t,const arb_t,int64_t), 
	     int64_t prec)
{
  static bool init=false;
  static arb_t low,high,tmp;
  if(!init)
    {
      init=true;
      arb_init(low);arb_init(high);arb_init(tmp);
    }
  prec=max(prec,32);
  arb_get_mid_arb(tmp,acb_realref(z));
  //printf("Midpoint = ");arb_printd(tmp,20);printf("\n");
  //printf("dt1 = ");arb_printd(dt1,10);printf("\n");
  arb_sub(low,tmp,dt1,prec);
  arb_add(high,tmp,dt1,prec);
  //printf("integrating from ");arb_printd(low,10);printf(" to ");arb_printd(high,10);printf("\n");
  //printf("maxd = ");arb_printd(maxd,20);
  molin_int(res,20,f1,maxd,low,high,prec);
  if(!arb_is_finite(res))
    {
      printf("Integral returned ");arb_printd(res,20);printf(" on [");arb_printd(low,20);printf(",");arb_printd(high,20);printf("]\n");
      exit(0);
    }
  printf(" Int from ");arb_printd(low,10);printf(" to ");arb_printd(high,10);printf(" returning ");arb_printd(res,10);printf("\n");fflush(stdout);
}

void do_int2(arb_t res, acb_t z, arb_t dt1,
	     void (*f1)(arb_t,const arb_t,int64_t), 
	     int64_t prec)
{
  static bool init=false;
  static arb_t low,high,tmp;
  if(!init)
    {
      init=true;
      arb_init(low);arb_init(high);arb_init(tmp);
    }
  prec=max(prec,32);
  arb_get_mid_arb(tmp,acb_imagref(z));
  //printf("Midpoint = ");arb_printd(tmp,20);printf("\n");
  //printf("dt1 = ");arb_printd(dt1,10);printf("\n");
  arb_sub(low,tmp,dt1,prec);
  if(arb_ge(low,six))
    {
      arb_zero(res);
      return;
    }
  arb_add(high,tmp,dt1,prec);
  if(arb_gt(high,six))
    arb_set(high,six);
  //printf("integrating from ");arb_printd(low,10);printf(" to ");arb_printd(high,10);printf("\n");
  //printf("maxd = ");arb_printd(maxd,20);
  molin_int(res,20,f1,maxd,low,high,prec);
  if(!arb_is_finite(res))
    {
      printf("Integral returned ");arb_printd(res,20);printf(" on [");arb_printd(low,20);printf(",");arb_printd(high,20);printf("]\n");
      exit(0);
    }
  printf(" Int from ");arb_printd(low,10);printf(" to ");arb_printd(high,10);printf(" returning ");arb_printd(res,10);printf("\n");fflush(stdout);
}
  

// here we assume z is a (small) complex box
// and we want to check that L''/L' is analytic over that box
// i.e. check that L' does not vanish
// we take L' at centre and +/- max(abs(L''))*radius
//

#define MAX_MAXD ((double) 100000.0)
bool check_for_poles(acb_t z, int64_t prec, arb_t dt1, bool verbose)
{
  static bool init=false;
  static acb_ptr res;
  static arb_t tmp,max_maxd;
  static acb_t z1;
  if(!init)
    {
      init=true;
      res=(acb_ptr)malloc(sizeof(acb_t)*4);
      acb_init(res);acb_init(res+1);acb_init(res+2);acb_init(res+3);
      arb_init(tmp);
      acb_init(z1);
      arb_init(max_maxd);arb_set_d(max_maxd,MAX_MAXD);
    }
  if(verbose){printf("z=");acb_printd(z,20);printf("\n");}
  if(verbose){printf("dt = ");arb_printd(dt1,20);printf("\n");}  
  acb_dirichlet_l_jet(res,z,dg,chi,0,4,8); // we think 8 bits is the sweet spot
  acb_mul_ui(res+3,res+3,6,prec);
  acb_abs(tmp,res+3,prec);
  if(verbose){printf("|L'''| = ");arb_printd(tmp,20);printf("\n");}
  arb_mul(tmp,tmp,dt1,prec);
  if(verbose){printf("dt*|L'''| = ");arb_printd(tmp,20);printf("\n");}
  arb_get_mid_arb(acb_realref(z1),acb_realref(z));
  arb_get_mid_arb(acb_imagref(z1),acb_imagref(z));
  acb_dirichlet_l_jet(res,z1,dg,chi,0,3,(prec<32 ? 32 : prec));
  acb_mul_2exp_si(res+2,res+2,1);
  if(verbose){printf("L'' at centre = ");acb_printd(res+2,20);printf("\n");}
  arb_add_error(acb_realref(res+2),tmp);
  arb_add_error(acb_imagref(res+2),tmp);
  if(verbose){printf("L'' with error = ");acb_printd(res+2,20);printf("\n");}
  acb_abs(tmp,res+2,prec);
  arb_mul(tmp,tmp,dt1,prec);
  if(verbose){printf("dt*|L''| = ");arb_printd(tmp,20);printf("\n");}
  if(verbose){printf("L' at centre = ");acb_printd(res+1,20);printf("\n");}
  arb_add_error(acb_realref(res+1),tmp);
  arb_add_error(acb_imagref(res+1),tmp);
  if(verbose){printf("L' with error = ");acb_printd(res+1,20);printf("\n");}
  acb_div(res,res+2,res+1,prec);
  if(verbose){printf("L''/L' = ");acb_printd(res,20);printf("\n");}
  acb_abs(maxd,res,prec);
  if(verbose){printf("maxd = ");arb_printd(maxd,20);printf("\n");}
  if(arb_is_finite(maxd)&&arb_le(maxd,max_maxd))
    return true;
  else
    return false;
}


#define MAX_DEPTH (30)

// do top and bottom from z-w1+iT to z+w1+iT
// z is exact
// w2 is 2* w1
void int_piecet(arb_t res, acb_t z, arb_t w1, arb_t w2, void (*f1)(arb_t,const arb_t,int64_t), uint64_t depth, int64_t prec)
{
  if(depth>=MAX_DEPTH)
    {
      printf("Maximum depth reached. Aborting.\n");
      printf("Problem is near z = ");acb_printd(z,20);printf("\n");
      acb_t z1;
      acb_init(z1);
      acb_set(z1,z);
      arb_add_error(acb_realref(z1),w2);
      arb_add_error(acb_imagref(z1),w2);
      check_for_poles(z1,prec,w2,true);
      exit(0);
    }
  //printf("In int_piecet with depth = %lu z = ",depth);
  //acb_printd(z,20);printf(" w1 = ");arb_printd(w1,20);printf("\n");fflush(stdout);
  acb_t z1;
  acb_init(z1);
  acb_set(z1,z);
  arb_add_error(acb_realref(z1),w2);
  arb_add_error(acb_imagref(z1),w2);
  if(check_for_poles(z1,prec,w2,false))
    {
      do_int1(res,z,w1,f1,prec);
      //printf("int1 returned ");arb_printd(res,20);printf("\n");fflush(stdout);
      acb_clear(z1);
      return;
    }
  arb_set(acb_imagref(z1),acb_imagref(z));
  arb_t w0;
  arb_init(w0);
  arb_mul_2exp_si(w0,w1,-1);
  // do left half
  arb_sub(acb_realref(z1),acb_realref(z),w0,prec);
  arb_t res1;
  arb_init(res1);
  int_piecet(res1,z1,w0,w1,f1,depth+1,prec);
  // do right half
  arb_add(acb_realref(z1),acb_realref(z),w0,prec);
  arb_t res2;
  arb_init(res2);
  int_piecet(res2,z1,w0,w1,f1,depth+1,prec);
  arb_add(res,res1,res2,prec);
  acb_clear(z1);
  arb_clear(w0);
  arb_clear(res1);
  arb_clear(res2);
  return;
}
// do left and right
// z is exact
// w2 is 2* w1
void int_piecel(arb_t res, acb_t z, arb_t w1, arb_t w2, void (*f1)(arb_t,const arb_t,int64_t), uint64_t depth, int64_t prec)
{
  if(depth>=MAX_DEPTH)
    {
      printf("Maximum depth reached. Aborting.\n");
      printf("Problem is near z = ");acb_printd(z,20);printf("\n");
      acb_t z1;
      acb_init(z1);
      acb_set(z1,z);
      arb_add_error(acb_realref(z1),w2);
      arb_add_error(acb_imagref(z1),w2);
      check_for_poles(z1,prec,w2,true);
      exit(0);
    }
  //printf("In int_piecel with depth = %lu z = ",depth);
  //acb_printd(z,20);printf(" w1 = ");arb_printd(w1,20);printf("\n");fflush(stdout);
  acb_t z1;
  acb_init(z1);
  acb_set(z1,z);
  arb_add_error(acb_imagref(z1),w2);
  arb_add_error(acb_realref(z1),w2);
  if(check_for_poles(z1,prec,w2,false))
    {
      do_int2(res,z,w1,f1,prec);
      //printf("int1 returned ");arb_printd(res,20);printf("\n");fflush(stdout);
      acb_clear(z1);
      return;
    }
  arb_set(acb_realref(z1),acb_realref(z));
  arb_t w0;
  arb_init(w0);
  arb_mul_2exp_si(w0,w1,-1);
  // do left half
  arb_sub(acb_imagref(z1),acb_imagref(z),w0,prec);
  arb_t res1;
  arb_init(res1);
  int_piecel(res1,z1,w0,w1,f1,depth+1,prec);
  // do right half
  arb_add(acb_imagref(z1),acb_imagref(z),w0,prec);
  arb_t res2;
  arb_init(res2);
  int_piecel(res2,z1,w0,w1,f1,depth+1,prec);
  arb_add(res,res1,res2,prec);
  acb_clear(z1);
  arb_clear(w0);
  arb_clear(res1);
  arb_clear(res2);
  return;
}



int main(int argc, char **argv)
{
  uint64_t depth=0;

  printf("Command line:- ");
  for(uint64_t n=0;n<argc;n++)
    printf("%s ",argv[n]);
  printf("\n");

  if(argc!=4)
    {
      printf("Usage:- %s <q> <chi> <prec>\n",argv[0]);
      return 0;
    }

  uint64_t q=atol(argv[1]);
  if((q%4)==2)
    {
      printf("There are no primitive characters for q = 2 mod 4.\n");
      return 0;
    }
  int64_t prec=atol(argv[3]);
  uint64_t cn=atol(argv[2]);

  arb_init(maxd);
  fmpz_t nz;
  fmpz_init(nz);

  arb_t two_pi;
  arb_init(two_pi);
  arb_const_pi(two_pi,prec);
  arb_mul_2exp_si(two_pi,two_pi,1);

  dirichlet_group_init(dg,q);

  dirichlet_char_init(chi,dg);
  dirichlet_char_first_primitive(chi,dg);
  while(cn!=dirichlet_char_exp(dg,chi))
    {
      if(dirichlet_char_next_primitive(chi,dg)<0)
	{
	  printf("No such primitive character %lu mod %lu.\n",cn,q);
	  return 0;
	}
    }

  /*
  if(dirichlet_parity_char(dg,chi)==1)
    {
      printf("Character %lu mod %lu is odd.\n",cn,q);
	return 0;
    }
  */

  

  switch(q)
    {
    case 3: arb_set_d(six,6.0);break;
    case 4: arb_set_d(six,5.0);break;
    case 5: 
    case 7: arb_set_d(six,4.0);break;
    default: arb_set_d(six,3.0);
    }
  printf("t0 set to ");arb_printd(six,20);printf("\n");
  

  acb_t z;
  acb_init(z);
  arb_t w1,w2,bottom;
  arb_init(w1);
  arb_init(w2);

  arb_set_d(half,0.5);
  arb_init(bottom);
  arb_mul_2exp_si(w1,half,-1);
  arb_mul_2exp_si(w2,w1,1);
  arb_set_d(acb_imagref(z),0.0);
  arb_set(acb_realref(z),w1);
  int_piecet(bottom,z,w1,w2,arb_fb,0,prec);
  printf("Int over bottom returned ");arb_printd(bottom,20);printf("\n");fflush(stdout);

  arb_t top;
  arb_init(top);
  arb_set(acb_imagref(z),six);
  int_piecet(top,z,w1,w2,arb_ft,0,prec);
  printf("Int over top returned ");arb_printd(top,20);printf("\n");fflush(stdout);

  arb_t left;
  arb_init(left);
  arb_zero(acb_realref(z));
  arb_mul_2exp_si(acb_imagref(z),six,-1);
  arb_set(w1,acb_imagref(z));
  arb_mul_2exp_si(w2,w1,1);
  int_piecel(left,z,w1,w2,arb_fl,0,prec);
  printf("Int over left returned ");arb_printd(left,20);printf("\n");fflush(stdout);

  arb_t right;
  arb_init(right);
  arb_set(acb_realref(z),half);
  arb_mul_2exp_si(acb_imagref(z),six,-1);
  arb_set(w1,acb_imagref(z));
  arb_mul_2exp_si(w2,w1,1);
  int_piecel(right,z,w1,w2,arb_fr,0,prec);
  printf("Int over right returned ");arb_printd(right,20);printf("\n");fflush(stdout);
  arb_t tmp;
  arb_init(tmp);
  arb_add(tmp,left,right,prec);
  arb_add(left,top,bottom,prec);
  arb_add(top,tmp,left,prec);
  arb_div(left,top,two_pi,prec);
  //printf("I = ");arb_printd(resl,20);printf("\n");
  if(arb_get_unique_fmpz(nz,left))
    {
      printf("result for q = %lu chi = %lu parity = %s I = ",q,cn,(dirichlet_parity_char(dg,chi)==1 ? "odd" : "even"));arb_printd(left,20);
      printf(" is N-P = ");
      fmpz_print(nz);
    }
  else
    {
      printf("result for q = %lu chi = %lu parity = %s I = ",q,cn,(dirichlet_parity_char(dg,chi)==1 ? "odd" : "even"));arb_printd(left,20);
      printf(" is I does not bracket a unique integer.");
    }

  printf("\n");fflush(stdout);
  return 0;
}
