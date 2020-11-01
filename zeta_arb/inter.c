#include "stdbool.h"
#include "acb.h"
#include "parameters.h"
#include "win_zeta.h"

// compute exp(-t^2/2H^2)
void inter_gaussian(arb_ptr res, arb_ptr t, int64_t prec)
{
  //arf_printf("inter_gaussian called with t=%.40Re\n",t);
  arb_div_d(res,t,H,prec);
  arb_mul(res,res,res,prec);
  arb_mul_2exp_si(res,res,-1);
  arb_neg(res,res);
  arb_exp(res,res,prec);
}

bool sincp;
// compute sinc(2*B*Pi*t) into ress, cos(2*B*Pi*t) into resc
// first time we call this (with sincp true) we compute sin and cos
// from then on (sincp false) just swap signs each time
void inter_sinc_cos(arb_ptr ress, arb_ptr resc, arb_ptr t, int64_t prec)
{
  static bool init=false;
  static arb_t sinc_tmp,msin,mcos,arb_two_pi_B;
  if(!init)
    {
      init=true;
      arb_init(sinc_tmp);
      arb_init(msin);
      arb_init(mcos);
      arb_init(arb_two_pi_B);
      arb_const_pi(sinc_tmp,prec);
      arb_div_d(arb_two_pi_B,sinc_tmp,INTER_A,prec);
    }
  arb_mul(sinc_tmp,t,arb_two_pi_B,prec);
  if(sincp)
    {
      arb_sin_cos(msin,mcos,sinc_tmp,prec);
      sincp=false;
    }
  else
    {
      arb_neg(msin,msin);
      arb_neg(mcos,mcos);
    }
  arb_set(resc,mcos);
  arb_div(ress,msin,sinc_tmp,prec);
}

// t is distance from t0 to t implied by f_ptr
// if fd_res is NULL, don't calc differential
void inter_point(arb_ptr f_res, arb_ptr fd_res, acb_t *f_vec, arb_ptr t, int f_ptr, int64_t prec)
{
  static bool init=false;
  static arb_t ip_tmp1,ip_tmp2,ip_tmp3,ip_tmp4;
  if(!init)
    {
      init=true;
      arb_init(ip_tmp1);
      arb_init(ip_tmp2);
      arb_init(ip_tmp3);
      arb_init(ip_tmp4);
    }
  arb_set(ip_tmp1,acb_realref(f_vec[f_ptr]));
  inter_gaussian(ip_tmp2,t,prec);
  inter_sinc_cos(ip_tmp3,ip_tmp4,t,prec);
  // compute f_res
  arb_mul(ip_tmp2,ip_tmp1,ip_tmp2,prec); //f*exp
  arb_mul(f_res,ip_tmp2,ip_tmp3,prec); //f*exp*sinc
  // compute fd_res
  if(!fd_res)
    return;
  arb_div(fd_res,ip_tmp4,t,prec); // cos/t
  arb_set_ui(ip_tmp4,1);
  arb_div(ip_tmp4,ip_tmp4,t,prec); // 1/t
  arb_div_d(ip_tmp1,t,H*H,prec); // t/H^2
  arb_add(ip_tmp4,ip_tmp4,ip_tmp1,prec); // 1/t+t/H^2
  arb_mul(ip_tmp4,ip_tmp4,ip_tmp3,prec); // sinc*(1/t+t/H^2)
  arb_sub(fd_res,fd_res,ip_tmp4,prec);   // cos/t-sinc*(1/t+t/H^2)
  arb_mul(fd_res,fd_res,ip_tmp2,prec);   // f*exp*(cos/t-sinc....)
  //
  // WARNING. COMPLETELY ARBITARY NEGATION TO MAKE f'(t) CORRECT SIGN
  //
  arb_neg(fd_res,fd_res);
}

// t0 points into f_vec
// let t=t0+(t_ptr-N/2)*one_over_A
// returns f(t) in f_res, f'(t) in fd_res
// return f_vec[t0]->re if t_ptr is integral 
void arb_inter_t(arb_ptr f_res, arb_ptr fd_res, acb_t *f_vec, arb_ptr t_ptr, int64_t prec)
{
  static bool init=false;
  static arb_t inter_tmp,inter_tmp1,inter_tmp2;
  if(!init)
    {
      init=true;
      arb_init(inter_tmp);
      arb_init(inter_tmp1);
      arb_init(inter_tmp2);
    }
  int i,j,nearest_t0;
  //printf("In inter_t\n");
  double t0_d=arb_get_d(t_ptr);
  arb_set_si(inter_tmp,t0_d);
  if(arb_equal(t_ptr,inter_tmp)) // t_ptr is integral so = a lattice pt
    nearest_t0=t0_d+1.5;
  else
    nearest_t0=t0_d+0.5;
  //printf("Nearest t0=%ld\n",nearest_t0);
  arb_zero(f_res);
  if (fd_res) arb_zero(fd_res);
  sincp=true;
  if(nearest_t0>t0_d) // start from right of t0
    {
      for(i=nearest_t0,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
      sincp=true;
      for(i=nearest_t0-INTER_SPACING,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
    }
  else // start from left
    {
      for(i=nearest_t0+INTER_SPACING,j=0;j<Ns;i+=INTER_SPACING,j++)
	{
	  //printf("Doing right tail with i=%d\n",i);
	  //
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  //printf("Calling ip\n");
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  //arf_printf("ip returned %.40Re\n",inter_tmp1);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  //
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	  //
	}
      sincp=true;
      for(i=nearest_t0,j=0;j<Ns;i-=INTER_SPACING,j++)
	{
	  arb_sub_ui(inter_tmp,t_ptr,i,prec);
	  arb_mul_d(inter_tmp,inter_tmp,-one_over_A,prec);
	  inter_point(inter_tmp1,inter_tmp2,f_vec,inter_tmp,i,prec);
	  arb_add(f_res,f_res,inter_tmp1,prec);
	  if (fd_res) arb_add(fd_res,fd_res,inter_tmp2,prec);
	}
    }
}

void acb_inter_t(arb_ptr f_res, acb_t *f_vec, arb_ptr t_ptr, int64_t prec)
{
  arb_inter_t(f_res, NULL, f_vec, t_ptr, prec);
}
