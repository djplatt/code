#include "acb.h"
#include "stdio.h"
#include "inttypes.h"
#include "stdbool.h"
#include "../quad/quad.h"

#define OP_ACC ((int64_t) 101)
arb_t *arb_zeros;
uint64_t zz;

// read our 13 byte structure representing a zero gap
// into an int_double
void in_bytes(arb_ptr res, FILE *infile, int64_t prec)
{
  uint64_t a;uint32_t b;uint8_t c;
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
  arb_set_ui(res,c);
  arb_mul_2exp_si(res,res,32);
  arb_add_ui(res,res,b,prec);
  arb_mul_2exp_si(res,res,64);
  arb_add_ui(res,res,a,prec);
  arb_mul_2exp_si(res,res,-OP_ACC);
  return;
}

void F_N(acb_t res, acb_t z, uint64_t N, arb_t *zeros, int64_t prec)
{
  static bool init=false;
  static acb_t tmp1,tmp2,tmp3;
  if(!init)
    {
      init=true;
      acb_init(tmp1);
      acb_init(tmp2);
      arb_set_d(acb_realref(tmp2),0.5);
      acb_init(tmp3);
    }
  acb_zero(res);
  for(uint64_t i=0;i<N;i++)
    {
      acb_mul_arb(tmp1,z,zeros[i],prec);
      acb_mul_onei(tmp3,tmp1);
      acb_exp(tmp1,tmp3,prec);
      arb_set(acb_imagref(tmp2),zeros[i]);
      acb_div(tmp3,tmp1,tmp2,prec);
      acb_add(res,res,tmp3,prec);
    }
  //printf("F_%lu(",N);acb_printd(z,20);printf(") returning ");acb_printd(res,20);printf("\n");
}

void Fd_N(acb_t res, acb_t z, uint64_t N, arb_t *zeros, int64_t prec)
{
  static bool init=false;
  static acb_t tmp1,tmp2,tmp3;
  if(!init)
    {
      init=true;
      acb_init(tmp1);
      acb_init(tmp2);
      arb_set_d(acb_realref(tmp2),0.5);
      acb_init(tmp3);
    }
  acb_zero(res);
  for(uint64_t i=0;i<N;i++)
    {
      acb_mul_arb(tmp1,z,zeros[i],prec);
      acb_mul_onei(tmp3,tmp1);
      acb_exp(tmp1,tmp3,prec);
      arb_set(acb_imagref(tmp2),zeros[i]);
      acb_div(tmp3,tmp1,tmp2,prec);
      acb_mul_arb(tmp1,tmp3,zeros[i],prec);
      acb_mul_onei(tmp3,tmp1);
      acb_add(res,res,tmp3,prec);
    }
  //printf("F_%lu(",N);acb_printd(z,20);printf(") returning ");acb_printd(res,20);printf("\n");
}

void newton(acb_t z1, acb_t z0, int64_t prec)
{
  static bool init=false;
  static acb_t tmp1,tmp2,tmp3;
  if(!init)
    {
      init=true;
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(tmp3);
    }
  F_N(tmp1,z0,zz,arb_zeros,prec);
  Fd_N(tmp2,z0,zz,arb_zeros,prec);
  acb_div(tmp3,tmp1,tmp2,prec);
  acb_sub(z1,z0,tmp3,prec);
}

void Fd_by_F_N(acb_t res, acb_t z, uint64_t N, arb_t *zeros, int64_t prec)
{
  static bool init=false;
  static acb_t tmp1,tmp2;
  if(!init)
    {
      init=true;
      acb_init(tmp1);
      acb_init(tmp2);
    }
  Fd_N(tmp1,z,N,zeros,prec);
  F_N(tmp2,z,N,zeros,prec);
  acb_div(res,tmp1,tmp2,prec);
  printf("F'/F(");acb_printd(z,20);printf(") rerurned ");acb_printd(res,20);printf("\n");
}

#define t_N (11)
#define y0 ((double) 6.69574675e-2)
#define y1 ((double) 1.21953870e-1)
#define u ((double) 1.468551615614841236e4)
#define v ((double) 7.983274721127996411e-2)
#define x0 ((double) u-0.05)
#define x1 ((double) u+0.05)

/*
#define C ((uint64_t) 10) // radius of contour around z

// F(t)=F'_N/F_N(z_0+e(t)/C)e(t)/C
void acb_quad_F(acb_t res, const acb_t t, int64_t prec) 
{
  static bool init=false;
  static acb_t z,tmp1,tmp2,tmp3;
  static arb_t two_pi;
  if(!init)
    {
      init=true;
      acb_init(z);
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(tmp3);
      arb_set_d(acb_realref(z),u);
      arb_set_d(acb_imagref(z),v);
      arb_init(two_pi);
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,1);
    }
  acb_mul_arb(tmp1,t,two_pi,prec);
  acb_mul_onei(tmp2,tmp1);
  acb_exp(tmp1,tmp2,prec);
  acb_div_ui(tmp2,tmp1,C,prec);
  acb_add(tmp3,tmp2,z,prec);
  Fd_by_F_N(tmp1,tmp3,zz,arb_zeros,prec);
  acb_mul(res,tmp2,tmp1,prec);
  printf("acb_quad_F(");acb_printd(t,20);printf(") returning ");acb_printd(res,20);printf("\n");
}

void arb_quad_F(arb_t res, const arb_t t, int64_t prec) 
{
  static bool init=false;
  static acb_t z,tmp1,tmp2,tmp3;
  static arb_t two_pi;
  if(!init)
    {
      init=true;
      acb_init(z);
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(tmp3);
      arb_set_d(acb_realref(z),u);
      arb_set_d(acb_imagref(z),v);
      arb_init(two_pi);
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,1);
    }
  arb_mul(acb_imagref(tmp2),two_pi,t,prec); // 2 pi t
  arb_zero(acb_realref(tmp2));
  acb_exp(tmp1,tmp2,prec); // e(t)
  acb_div_ui(tmp3,tmp1,C,prec); // e(t)/C
  acb_add(tmp2,tmp3,z,prec); // e(t)/C+z
  Fd_by_F_N(tmp1,tmp2,zz,arb_zeros,prec); // F'/F(e(t)+z)
  acb_mul(tmp2,tmp3,tmp1,prec); // F'/F(e(t)+z) * e(t)/C
  arb_set(res,acb_realref(tmp2));
}

void do_argument(int64_t prec)
{
  arb_t mx,lo,hi,res;
  arb_init(mx);arb_init(lo);arb_init(hi);arb_init(res);
  arb_zero(lo);
  arb_const_pi(hi,prec);
  arb_mul_2exp_si(hi,hi,1);
  arb_maxd(mx,acb_quad_F,lo,hi,prec,10);
  printf("maxd = ");arb_printd(mx,20);printf("\n");
  molin_int(res, 100, arb_quad_F, mx, lo,hi,prec);
  printf("int = ");arb_printd(res,20);printf("\n");
 
}
*/

void my_acb_arg(arb_t res, acb_t z, int64_t prec)
{
  static bool init=false;
  static acb_t tmp;
  static arb_t pi;
  if(!init)
    {
      init=true;
      acb_init(tmp);
      arb_init(pi);
      arb_const_pi(pi,prec);
    }
  if(arb_contains_zero(acb_imagref(z)))
    if(arb_is_negative(acb_realref(z)))
      {
	acb_neg(tmp,z);
	acb_arg(res,tmp,prec);
	arb_add(res,res,pi,prec);
	return;
      }
  acb_arg(res,z,prec);
}

#define F_ball ((int64_t) -20)
void do_argument1(arb_t min_F, acb_t z, int64_t prec)
{
  static bool init=false;
  static acb_t tmp1,tmp2,tmp3;
  static arb_t two_pi,t1,t2,t3;
  if(!init)
    {
      init=true;
      acb_init(tmp1);
      acb_init(tmp2);
      acb_init(tmp3);
      arb_init(t1);
      arb_init(t2);
      arb_init(t3);
      arb_set_d(acb_realref(z),u);
      arb_set_d(acb_imagref(z),v);
      arb_init(two_pi);
      arb_const_pi(two_pi,prec);
      arb_mul_2exp_si(two_pi,two_pi,1);
    }
  arb_set_ui(min_F,1000);
  double dt=1.0/1024;
  arb_set_d(t1,dt);
  for(uint64_t i=1;i<1024;i+=2)
    {
      arb_set_ui(t2,i);
      arb_mul_2exp_si(t2,t2,-10);
      arb_add_error(t2,t1);
      arb_zero(acb_realref(tmp1));
      arb_mul(acb_imagref(tmp1),t2,two_pi,prec);
      acb_exp(tmp2,tmp1,prec);
      acb_mul_2exp_si(tmp2,tmp2,F_ball);
      acb_add(tmp1,tmp2,z,prec);
      F_N(tmp2,tmp1,zz,arb_zeros,prec);
      acb_abs(t3,tmp2,prec);
      arb_min(min_F,min_F,t3,prec);
      my_acb_arg(t3,tmp2,prec);
      //printf("Arg(F_%lu(",zz);acb_printd(tmp1,20);printf(" = ");arb_printd(t3,20);printf("\n");
    }
}

void do_alpha(arb_t res, double dx, double dy, arb_t *zeros, int64_t prec)
{
  arb_t tmp,tmp1;arb_init(tmp);arb_init(tmp1);
  arb_set_ui(res,1);
  acb_t z,fn;
  acb_init(z);acb_init(fn);
  arb_set_d(acb_realref(z),x0);
  arb_set_d(acb_imagref(z),y0+dy/2.0);
  arb_t err;
  arb_init(err);
  arb_set_d(err,dy/2.0);
  arb_add_error(acb_imagref(z),err);
  arb_mul_2exp_si(err,err,1);
  arb_t ay1;arb_init(ay1);arb_set_d(ay1,y1);
  while(true)
    {
      arb_sub(tmp,ay1,acb_imagref(z),prec);
      if(arb_is_negative(tmp)) break;
      F_N(fn,z,t_N,zeros,prec);
      acb_abs(tmp,fn,prec);
      arb_min(tmp1,res,tmp,prec);
      arb_swap(tmp1,res);
      arb_add(acb_imagref(z),acb_imagref(z),err,prec);
    }
  arb_set_d(acb_realref(z),x1);
  arb_set_d(acb_imagref(z),y0+dy/2.0);
  arb_set_d(err,dy/2.0);
  arb_add_error(acb_imagref(z),err);
  arb_mul_2exp_si(err,err,1);
  while(true)
    {
      arb_sub(tmp,ay1,acb_imagref(z),prec);
      if(arb_is_negative(tmp)) break;
      F_N(fn,z,t_N,zeros,prec);
      acb_abs(tmp,fn,prec);
      arb_min(tmp1,res,tmp,prec);
      arb_swap(tmp1,res);
      arb_add(acb_imagref(z),acb_imagref(z),err,prec);
    }
  arb_set_d(acb_realref(z),x0+dx/2.0);
  arb_set_d(acb_imagref(z),y0);
  arb_set_d(err,dx/2.0);
  arb_add_error(acb_realref(z),err);
  arb_mul_2exp_si(err,err,1);
  arb_t ax1;arb_init(ax1);arb_set_d(ax1,x1);
  while(true)
    {
      arb_sub(tmp,ax1,acb_realref(z),prec);
      if(arb_is_negative(tmp)) break;
      F_N(fn,z,t_N,zeros,prec);
      acb_abs(tmp,fn,prec);
      arb_min(tmp1,res,tmp,prec);
      arb_swap(tmp1,res);
      arb_add(acb_realref(z),acb_realref(z),err,prec);
    }
  arb_set_d(acb_realref(z),x0+dx/2.0);
  arb_set_d(acb_imagref(z),y1);
  arb_set_d(err,dx/2.0);
  arb_add_error(acb_realref(z),err);
  arb_mul_2exp_si(err,err,1);
  while(true)
    {
      arb_sub(tmp,ax1,acb_realref(z),prec);
      if(arb_is_negative(tmp)) break;
      F_N(fn,z,t_N,zeros,prec);
      acb_abs(tmp,fn,prec);
      arb_min(tmp1,res,tmp,prec);
      arb_swap(tmp1,res);
      //printf("Max now ");arb_printd(res,20);printf("\n");
      arb_add(acb_realref(z),acb_realref(z),err,prec);
    }
}

void do_a1(arb_t res, arb_t gamma, arb_t y00, int64_t prec)
{
  static bool init=false;
  static arb_t t1,t2; 
  static acb_t tmp1,tmp2;
  if(!init)
    {
      init=true;
      acb_init(tmp1);
      arb_init(t1);
      arb_init(t2);
    }
  arb_mul(t1,y00,gamma,prec);
  arb_neg(t1,t1);
  arb_exp(t2,t1,prec);
  arb_set_d(acb_realref(tmp1),0.5);
  arb_set(acb_imagref(tmp1),gamma);
  acb_abs(t1,tmp1,prec);
  arb_div(res,t2,t1,prec);
}



void do_a(arb_t res, arb_t y00, uint64_t N, arb_t *zeros, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_zero(res);
  for(uint64_t i=0;i<N;i++)
    {
      do_a1(tmp,zeros[i],y00,prec);
      arb_add(res,res,tmp,prec);
    }
}

void do_b(arb_t res, arb_t y00, uint64_t N, uint64_t zz, arb_t *zeros, int64_t prec)
{
  static bool init=false;
  static arb_t tmp;
  if(!init)
    {
      init=true;
      arb_init(tmp);
    }
  arb_zero(res);
  for(uint64_t i=N;i<zz;i++)
    {
      do_a1(tmp,zeros[i],y00,prec);
      arb_add(res,res,tmp,prec);
    }
}

bool ok_q0(uint64_t q0, arb_t a, arb_t aw, arb_t b, arb_t bw, arb_t alpha, int64_t prec)
{
  static bool init=false;
  static arb_t lhs,rhs,t1,t2,t3,pi;
  if(!init)
    {
      init=true;
      arb_init(lhs);
      arb_init(rhs);
      arb_init(t1);
      arb_init(t2);
      arb_init(t3);
      arb_init(pi);
      arb_const_pi(pi,prec);
    }
  arb_set_ui(t1,q0);
  arb_const_pi(t2,prec);
  arb_div(t3,pi,t1,prec); // pi/q0
  arb_sin(t1,t3,prec);
  arb_mul_2exp_si(t1,t1,1);
  arb_mul(t3,t1,aw,prec); // t3 <- 2aw sin(pi/q0)
  arb_mul_2exp_si(t2,pi,1); // 2 pi
  arb_mul(t1,t2,a,prec); // 2 pi a
  arb_div_ui(t2,t1,q0,prec); // 2 pi a /q0
  arb_add(lhs,t3,t2,prec);
  arb_sub(t1,alpha,b,prec);
  arb_mul_2exp_si(t2,bw,1);
  arb_sub(rhs,t1,t2,prec);
  arb_sub(t1,lhs,rhs,prec);
  return arb_is_positive(t1);
}


void do_kacz(arb_t *zeros, uint64_t NN, double tmax, int64_t prec)
{
  printf("In do_kacz with %lu zeros up to height %f.\n",NN,tmax);
  acb_t z,res;acb_init(z);acb_init(res);
  arb_set_d(acb_realref(z),u);
  arb_set_d(acb_imagref(z),v);
  newton(z,z,prec); // only worth one iteration
  printf("z = ");acb_printd(z,30);printf("\n");

  arb_t alpha;
  arb_init(alpha);
  arb_t min_F;
  arb_init(min_F);
  do_argument1(min_F,z,prec);
  printf("Minimum F_N(z) around C_1 was ");arb_printd(min_F,20);printf("\n");
  do_alpha(alpha,0.000001,0.000001,zeros,prec);
  printf("alpha = ");arb_printd(alpha,20);printf("\n");
  arb_t a,aw,b,bw;
  arb_init(a);
  arb_init(aw);
  arb_init(b);
  arb_init(bw);
  arb_set_d(aw,y0);
  do_a(a,aw,t_N,arb_zeros,prec);
  do_b(b,aw,t_N,zz,arb_zeros,prec);
  printf("a = ");arb_printd(a,20);
  printf("\nb = ");arb_printd(b,20);printf("\n");
  arb_set_ui(aw,1);
  arb_mul_2exp_si(aw,aw,F_ball);
  arb_add_error(acb_realref(z),aw);
  arb_add_error(acb_imagref(z),aw);
  do_a(aw,acb_imagref(z),t_N,arb_zeros,prec);
  do_b(bw,acb_imagref(z),t_N,zz,arb_zeros,prec);
  printf("aw = ");arb_printd(aw,20);
  printf("\nbw = ");arb_printd(bw,20);printf("\n");
  uint64_t q0;
  for(q0=10;;q0++)
    if(!ok_q0(q0,a,aw,b,bw,alpha,prec)) break;
  printf("q0 = %lu\n",q0);
  arb_set_ui(a,q0);
  arb_log(b,a,prec);
  arb_mul_ui(a,b,t_N,prec);
  arb_neg(a,a);
  arb_exp(b,a,prec);
  arb_mul_2exp_si(b,b,1);
  printf("2 var kappa = ");arb_printd(b,20);printf("\n");


}

int main(int argc, char ** argv)
{
  printf("Command line: ");
  for(int64_t i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=3) exit(0);
  FILE*zfile=fopen(argv[1],"rb");
  if(!zfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  int64_t prec=atol(argv[2]);

  long int num_its;
  fread(&num_its,sizeof(long int),1,zfile);
  if(num_its>1)
    {
      printf("Hard wired to read only first block of zeros.\n");
      num_its=1;
    }
  
  double st[2];
  long int zs[2];
  int64_t ptr=0;
  arb_t del_t,t,rho,t_err;
  arb_init(del_t);arb_init(t);arb_init(rho);arb_init(t_err);
  arb_set_ui(t_err,1);
  arb_mul_2exp_si(t_err,t_err,-OP_ACC-1);
  for(uint64_t it=0;it<num_its;it++)
    {
      fread(st,sizeof(double),2,zfile);
      fread(&zs[0],sizeof(long int),1,zfile);
      if(st[0]==0.0)
	continue;
      arb_set_d(t,st[0]);
      fread(&zs[1],sizeof(long int),1,zfile);
      arb_zeros=(arb_t *)malloc(sizeof(arb_t)*(zs[1]-zs[0]));
      if(!arb_zeros)
	{
	  printf("Error allocating memory for arb_zeros. Exiting.\n");
	  exit(0);
	}
      zz=0;
      for(uint64_t z=zs[0]+1;z<=zs[1];z++)
	{
	  in_bytes(del_t,zfile,prec);
	  if(arb_is_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  arb_add(t,t,del_t,prec);
	  arb_set(rho,t);
	  arb_add_error(rho,t_err);
	  arb_init(arb_zeros[zz]);
	  arb_set(arb_zeros[zz],rho);
	  zz++;
	}
 
      do_kacz(arb_zeros,zz,st[1],prec);

      for(uint64_t i=0;i<zz;i++)
	arb_clear(arb_zeros[i]);
      free(arb_zeros);
    }
  return 0;
}
