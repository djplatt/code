/*
Code to implement Turing's method.
See Booker: Artin's Conjecture, Turing's Method, and the Riemann Hypothesis - Exp. Math.
    Platt: Isolating some non-trivial zeros of zeta - Math. Comp.
    Platt & Trudgian: The Riemann hypothesis is true up to 3\cdot 10^{12} - To appear 
*/

#include "stdbool.h"
#include "acb.h"
#include "parameters.h"
#include "win_zeta.h"
#include "inter.h"

/*
On RH, there are no local minima of completed zeta above the real line
and no local maxima below it. Thus, if we find something that looks
like on,e then we should look harder as it probably hides two zeros.
*/
 
// we have found a stationary pt between left, left+1 and left+2
// we don't assume RH, so go find the zeros
bool resolve_stat_point(arb_ptr tl, arb_ptr ftl, arb_ptr tm, arb_ptr ftm, arb_ptr tr, arb_ptr ftr, int left, sign_t this_sign, acb_t *f_vec, int64_t prec)
{
  static bool init=false;
  static arb_t sp_t,sp_ft,sp_t1,sp_ft1;  
  if(!init)
    {
      init=true;
      arb_init(sp_t);
      arb_init(sp_ft);
      arb_init(sp_t1);
      arb_init(sp_ft1);
    }
  dir_t dir1,dir2,dir3,dir4;
  arb_set_ui(tl,left);
  arb_set_ui(tm,left+1);
  arb_set_ui(tr,left+2);
  arb_set(ftl,acb_realref(f_vec[left]));
  arb_set(ftm,acb_realref(f_vec[left+1]));
  arb_set(ftr,acb_realref(f_vec[left+2]));

  while(true)
    {
      arb_add(sp_t,tl,tm,prec);
      arb_mul_2exp_si(sp_t,sp_t,-1); // sp_t <- (tl+tm)/2
      arb_inter_t(sp_ft,NULL,f_vec,sp_t,prec); // sp_ft <- f(sp_t)
      if((arb_sign(sp_ft)&this_sign)==0) // different signs, found them
	{
	  arb_set(tr,tm);
	  arb_set(ftr,ftm);
	  arb_set(tm,sp_t);
	  arb_set(ftm,sp_ft);
	  return true;
	}
      dir1=arb_dir(ftl,sp_ft);
      dir2=arb_dir(sp_ft,ftm);
      if(((this_sign==POS)&&(dir1==DOWN)&&(dir2==UP))||((this_sign==NEG)&&(dir1==UP)&&(dir2==DOWN)))
	{
	  arb_set(tr,tm);
	  arb_set(ftr,ftm);
	  arb_set(tm,sp_t);
	  arb_set(ftm,sp_ft);
	  continue;
	}

      arb_add(sp_t1,tm,tr,prec);
      arb_mul_2exp_si(sp_t1,sp_t1,-1);
      arb_inter_t(sp_ft1,NULL,f_vec,sp_t1,prec);
      if((arb_sign(sp_ft1)&this_sign)==0)
	{
	  arb_set(tl,tm);
	  arb_set(ftl,ftm);
	  arb_set(tm,sp_t1);
	  arb_set(ftm,sp_ft1);
	  return true;
	}
      dir3=arb_dir(ftm,sp_ft1);
      dir4=arb_dir(sp_ft1,ftr);
      if(((this_sign==POS)&&(dir2==DOWN)&&(dir3==UP))||((this_sign==NEG)&&(dir2==UP)&&(dir3==DOWN)))
	{
	  arb_set(tl,sp_t);
	  arb_set(ftl,sp_ft);
	  arb_set(tr,sp_t1);
	  arb_set(ftr,sp_ft1);
	  continue;
	}
      if(((this_sign==POS)&&(dir3==DOWN)&&(dir4==UP))||((this_sign==NEG)&&(dir3==UP)&&(dir4==DOWN)))
	{
	  arb_set(tl,tm);
	  arb_set(ftl,ftm);
	  arb_set(tm,sp_t1);
	  arb_set(ftm,sp_ft1);
	  continue;
	}
      printf("Stat point resolver failed to converge.\n");
      return false; // failed to converge
    }
}

// true if there is a stationary point l-m-r
bool stat_pt(acb_ptr l, acb_ptr m, acb_ptr r)
{
  sign_t sm=acb_sign(m);
  if(sm==UNK) // can't tell
    return false;
  if(sm==POS) // look for l>m and r>m
    return arb_gt(acb_realref(l),acb_realref(m))&&arb_gt(acb_realref(r),acb_realref(m));
  else // look for m>l and m>r
    return arb_gt(acb_realref(m),acb_realref(l))&&arb_gt(acb_realref(m),acb_realref(r));
}


// count zeros (using stationary points) between start and end
// this could be more efficient as we end up taking sign of
// entries multiple times, but it's insignificant anyway
int zeros_st(acb_t *f_vec, int start, int end, int64_t prec)
{
  static bool init=false;
  static arb_t t1,t2,t3,t4,t5,t6;
  if(!init)
    {
      init=true;
      arb_init(t1);
      arb_init(t2);
      arb_init(t3);
      arb_init(t4);
      arb_init(t5);
      arb_init(t6);
    }
  // start by looking for sign changes
  int ptr=start,count=0;
  sign_t last_sign=acb_sign(f_vec[ptr++]);
  while(last_sign==UNK)
    if(ptr>end)
      return 0; // the vector was all UNK, a problem
    else
      last_sign=acb_sign(f_vec[ptr++]);
  while(ptr<=end)
    {
      sign_t this_sign=acb_sign(f_vec[ptr++]);
      if(this_sign==UNK)
	continue;
      if(this_sign==last_sign)
	continue;
      count++;
      last_sign=this_sign;
    }
  // now look for stationary points
  ptr=start+2;
  while(ptr<=end)
    {
      if(stat_pt(f_vec[ptr-2],f_vec[ptr-1],f_vec[ptr])) // stat pt
	if(resolve_stat_point(t1,t2,t3,t4,t5,t6,ptr-2,arb_sign(acb_realref(f_vec[ptr-2])),f_vec,prec)) // found two zeros
	  count+=2;
      ptr++;
    }
  return count;
}

// used by im_int
void im_int1(arb_t res, arb_t t, int64_t prec)
{
  static bool init=false;
  static arb_t im_t,im_1;
  if(!init)
    {
      init=true;
      arb_init(im_t);
      arb_init(im_1);
    }
  arb_mul_2exp_si(im_t,t,2); // 4t
  arb_atan(res,im_t,prec);
  //printf("atan(4t)=");mpfi_printn(res,30);
  arb_mul(res,res,t,prec);
  arb_mul_d(res,res,-0.25,prec); // -t/4*atan(4t)
  //printf("-t/4atan(4t)=");mpfi_printn(res,30);
  arb_mul(im_t,t,t,prec); // t^2
  arb_mul_d(im_1,im_t,-3.0/4.0,prec);
  arb_add(res,res,im_1,prec); //-3/4*t^2
  arb_mul_ui(im_1,im_t,16,prec); // 16t^2
  arb_add_ui(im_1,im_1,1,prec);
  arb_log(im_1,im_1,prec); // log(1+16t^2)
  arb_mul_d(im_1,im_1,1.0/32.0,prec);
  arb_add(res,res,im_1,prec);
  arb_add_d(res,res,-1.0/64.0,prec);
  arb_add_d(im_t,im_t,1.0/16.0,prec); // t^2+1/16
  arb_log(im_1,im_t,prec);
  arb_mul(im_t,im_t,im_1,prec); // (t^2+1/16)log(t^2+1/16)
  arb_mul_d(im_t,im_t,0.25,prec);
  arb_add(res,res,im_t,prec);
}

//int_t0^t1 Im log gamma(1/4+it/2) dt
void im_int(arb_t res, arb_t t0, arb_t t1, int64_t prec)
{
  static bool init=false;
  static arb_t im_err,im_t0,im_t1,im_1,im_2;
  if(!init)
    {
      init=true;
      arb_init(im_err);
      arb_init(im_t0);
      arb_init(im_t1);
      arb_init(im_1);
      arb_init(im_2);
    }
  // |error| is < int 1/(6|z|) < (t1-t0)/(3 t0) 
  arb_sub(im_1,t1,t0,prec);
  arb_div_ui(im_2,im_1,3,prec);
  arb_div(im_err,im_2,t0,prec);
  arb_mul_2exp_si(im_t0,t0,-1);
  arb_mul_2exp_si(im_t1,t1,-1);
  im_int1(res,im_t1,prec);
  im_int1(im_2,im_t0,prec);
  arb_sub(res,res,im_2,prec);
  arb_mul_ui(res,res,2,prec);
  arb_add_error(res,im_err);
}


// assume that if a zero lies between t0 and t1 then it is at t0
double Nleft_int(long int t0_ptr, long int t1_ptr, acb_t *f_vec, double delta, int64_t prec)
{
  static bool init=false;
  static arb_t t1,t2,t3,t4,t5,t6;
  if(!init)
    {
      init=true;
      arb_init(t1);
      arb_init(t2);
      arb_init(t3);
      arb_init(t4);
      arb_init(t5);
      arb_init(t6);
    }
  long int res=0;
  long int ptr=t0_ptr,last_ptr;
  sign_t last_sign=acb_sign(f_vec[ptr++]),this_sign;
  while(last_sign==UNK)
    last_sign=acb_sign(f_vec[ptr++]);
  last_ptr=ptr-1;
  for(;ptr<=t1_ptr;)
    {
      this_sign=acb_sign(f_vec[ptr++]);
      if(this_sign==UNK) // might be a sign change coming
	continue;
      if(this_sign!=last_sign) // a sign change after last_ptr 
	{
	  res-=(last_ptr-t0_ptr);
	  last_ptr=ptr-1;
	  last_sign=this_sign;
	  continue;
	}
      if(ptr-3>=t0_ptr) // far enough in to find stat points
	if(stat_pt(f_vec[ptr-3],f_vec[ptr-2],f_vec[ptr-1])) // suspected stat pt
	  if(resolve_stat_point(t1,t2,t3,t4,t5,t6,ptr-3,arb_sign(acb_realref(f_vec[ptr-2])),f_vec,prec)) // found two zeros after ptr-3
	    res-=2*(ptr-3-t0_ptr);
      last_ptr=ptr-1;
    }
  //printf("Nleft_int returning %f\n",res*delta);
  return res*delta;
}

// assume that if a zero lies between t0 and t1 then it is at t1
double Nright_int(long int t0_ptr, long int t1_ptr, acb_t *f_vec, double delta, int64_t prec)
{
  static bool init=false;
  static arb_t t1,t2,t3,t4,t5,t6;
  if(!init)
    {
      init=true;
      arb_init(t1);
      arb_init(t2);
      arb_init(t3);
      arb_init(t4);
      arb_init(t5);
      arb_init(t6);
    }
  long int res=0;
  long int ptr=t0_ptr;
  sign_t last_sign=acb_sign(f_vec[ptr++]),this_sign;
  while(last_sign==UNK)
    last_sign=acb_sign(f_vec[ptr++]);
  for(;ptr<=t1_ptr;)
    {
      this_sign=acb_sign(f_vec[ptr++]);
      if(this_sign==UNK) // no use if we can't tell sign
	continue;
      if(this_sign!=last_sign) // a sign change before ptr-1
	{
	  res+=(t1_ptr-ptr+1);
	  last_sign=this_sign;
	  continue;
	}

      if(ptr-3>=t0_ptr) // check for stat pt
	if(stat_pt(f_vec[ptr-3],f_vec[ptr-2],f_vec[ptr-1]))
	  if(resolve_stat_point(t1,t2,t3,t4,t5,t6,ptr-3,arb_sign(acb_realref(f_vec[ptr-2])),f_vec,prec)) // 2 zeros before ptr-1
	    res+=2*(t1_ptr-ptr+1); 
    }  
  return res*delta;
}

// returns bound on int_t0^{t0+h} S(t) dt
// t0=t0_ptr*delta
// t0>168*Pi
// Trudgian Thesis Theorem 5.2.2
// < 0.059 log(t0+h) + 2.067
void St_int(arb_t res, arb_t t, int64_t prec) 
{
  static bool init=false;
  static arb_t c1,c2;
  if(!init)
    {
      init=true;
      arb_init(c1);
      arb_init(c2);
      arb_set_ui(c1,59);
      arb_div_ui(c1,c1,1000,prec);
      arb_set_ui(c2,2067);
      arb_div_ui(c2,c2,1000,prec);
    }
  arb_log(res,t,prec);
  arb_mul(res,res,c1,prec);
  arb_add(res,res,c2,prec);
}

// the log term integrated
void ln_term(arb_t res, arb_t t0, arb_t t1, arb_t h1, int64_t prec)
{
  static bool init=false;
  static arb_t ln_pi;
  if(!init)
    {
      init=true;
      arb_init(ln_pi);
      arb_const_pi(res,prec);
      arb_log(ln_pi,res,prec);
    }
  arb_add(res,t0,t1,prec);
  //arb_mul_d(res,res,-0.25,prec);
  arb_mul_2exp_si(res,res,-2);
  arb_neg(res,res);
  arb_mul(res,res,ln_pi,prec);
  arb_mul(res,res,h1,prec);
}

// the maximum number of zeros <=a based on Turing region [a,b]
long int turing_max(acb_t *f_vec, long int a_ptr, double a, long int b_ptr, double b, int64_t prec)
{
  static bool init=false;
  static fmpz_t fz;
  static arb_t pi,tm_t0,tm_t1,tm_h,tm_tmp,tm_tmp1,tm_res;
  if(!init)
    {
      init=true;
      arb_init(pi);
      arb_const_pi(pi,prec);
      fmpz_init(fz);
      arb_init(tm_t0);
      arb_init(tm_t1);
      arb_init(tm_h);
      arb_init(tm_tmp);
      arb_init(tm_tmp1);
      arb_init(tm_res);
    }
  arb_set_d(tm_t0,a);
  arb_set_d(tm_t1,b);
  arb_sub(tm_h,tm_t1,tm_t0,prec);  
  St_int(tm_res,tm_t1,prec);
  arb_sub_d(tm_res,tm_res,Nright_int(a_ptr,b_ptr,f_vec,one_over_A,prec),prec);
  ln_term(tm_tmp,tm_t0,tm_t1,tm_h,prec);
  im_int(tm_tmp1,tm_t0,tm_t1,prec);
  arb_add(tm_tmp,tm_tmp,tm_tmp1,prec);
  arb_div(tm_tmp,tm_tmp,pi,prec);
  arb_add(tm_res,tm_res,tm_tmp,prec);
  arb_div(tm_tmp,tm_res,tm_h,prec);
  arb_floor(tm_res,tm_tmp,prec);
  if(!arb_get_unique_fmpz(fz,tm_res)) // should never happen
    {
      printf("Turing max did not bracket an integer. Exiting.\n");
      exit(0);
    }
  return fmpz_get_si(fz)+1; // +1 because of sqrt(omega) =+/-1
}

// the minimum number of zeros below b based on Turing region [a,b]
long int turing_min(acb_t *f_vec, long int a_ptr, double a, long int b_ptr, double b, int64_t prec)
{
  static bool init=false;
  static fmpz_t fz;
  static arb_t pi,tm_t0,tm_t1,tm_h,tm_tmp,tm_tmp1,tm_res;
  if(!init)
    {
      init=true;
      arb_init(pi);
      arb_const_pi(pi,prec);
      fmpz_init(fz);
      arb_init(tm_t0);
      arb_init(tm_t1);
      arb_init(tm_h);
      arb_init(tm_tmp);
      arb_init(tm_tmp1);
      arb_init(tm_res);
    }
  arb_set_d(tm_t0,a);
  arb_set_d(tm_t1,b);;
  arb_sub(tm_h,tm_t1,tm_t0,prec);
  St_int(tm_res,tm_t1,prec);
  arb_neg(tm_res,tm_res);
  arb_sub_d(tm_res,tm_res,Nleft_int(a_ptr,b_ptr,f_vec,one_over_A,prec),prec);
  ln_term(tm_tmp,tm_t0,tm_t1,tm_h,prec);
  im_int(tm_tmp1,tm_t0,tm_t1,prec);
  arb_add(tm_tmp,tm_tmp,tm_tmp1,prec);
  arb_div(tm_tmp,tm_tmp,pi,prec);
  arb_add(tm_res,tm_res,tm_tmp,prec);
  arb_div(tm_res,tm_res,tm_h,prec);
  arb_ceil(tm_res,tm_res,prec);
  if(!arb_get_unique_fmpz(fz,tm_res))
    {
      printf("Turing min did not bracket an integer. Exiting.\n");
      exit(0);
    }
  return(fmpz_get_si(fz)+1); // +1 because of sqrt(omega) =+/-1 
}

// Use Turing's method to estimate and then check number of zeros in [a,b]
// last_max is zeros to a (if known, 0 otherwise)
// on success returns the number of zeros of zeta up to b
// on failure return 0
// a_ptr,b_ptr are pointers to a,b respectively in f_vec
long int turing(acb_t *f_vec, long int a_ptr, double a, long int b_ptr, double b, long int last_max, int64_t prec)
{
  // estimate maximum for N(b)
  long int min_lo,max_hi=turing_max(f_vec,b_ptr,b,b_ptr+TURING_LEN,b+TURING_WIDTH,prec);
  if(last_max==0) // No previous run or previous run failed, so estimate min N(a)
    min_lo=turing_min(f_vec,a_ptr-TURING_LEN,a-TURING_WIDTH,a_ptr,a,prec);
  else // If previous run succeeded, zeros to height N(a) is known
    min_lo=last_max;
  printf("looking for %ld-%ld=%ld zeros\n",max_hi,min_lo,max_hi-min_lo);
  long int num_found,num_exp=max_hi-min_lo;
  sign_t ls,rs;
  ls=acb_sign(f_vec[a_ptr]);
  if(ls==UNK)
    {
      printf("Unknown at start of data.\n");
      printf("Missed All/All zeros in region %f to %f.\n",a,b);
      return 0;
    }
  rs=acb_sign(f_vec[b_ptr]);
  if(rs==UNK)
    {
      printf("Unknown at end of data.\n");
      printf("Missed All/All zeros in region %f to %f.\n",a,b);
      return 0;
    }

  if(ls==rs) // same sign at both ends, count should be even
    {
      if(num_exp&1)
	{
	  printf("Problem in Turing Zone. Num zeros should be even.\n");
	  printf("Missed All/All zeros in region %f to %f.\n",a,b);
	  return 0;
	}
    }
  else // count should be odd
    {
      if(!(num_exp&1)) 
	{
	  printf("Problem in Turing Zone. Num zeros should be odd.\n");
	  printf("Missed All/All zeros in region %f to %f.\n",a,b);
	  return 0;
	}
    }

  num_found=zeros_st(f_vec,a_ptr,b_ptr,prec); // go find zeros, using stat pts as well

  if(num_exp==num_found) // excellent
    {
      printf("All %ld zeros found in region %f to %f using stat points.\n",num_exp,a,b);
      return max_hi; // use this for next iteration
    }
  // bugger. Turing's method failed to account for all the zeros.
  // could be. 1) we missed some in [a,b]. look harder next time
  //           2) we missed some in [a-tz,a] and or [b,b+tz]
  //              where tz is the Turing Zone. Either look harder
  //              or shift a,b 
  printf("Missed %ld/%ld zeros (running with stat points) in region %f to %f.\n",num_exp-num_found,num_exp,a,b);
  return 0;
}
