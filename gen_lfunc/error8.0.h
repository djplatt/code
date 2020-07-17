bool M_error(arb_t res, arb_t x, uint64_t M, L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  static bool init=false;
  static arb_t C_bit,arb_M,tmp1,tmp2,tmp3,tmp4,u,X,half_log_N,pi,XM2r,Mc;
  if(!init)
    {
      init=true;
      arb_init(arb_M);
      arb_init(C_bit);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(u);
      arb_init(X);
      arb_init(half_log_N);
      arb_init(pi);
      arb_init(XM2r);
      arb_init(Mc);
      arb_init(tmp4);
      // pi = pi
      arb_const_pi(pi,prec);
      // C_bit = C2^(r/2)
      arb_sqrt_ui(tmp1,3,prec);
      arb_set_ui(tmp2,2);
      arb_set_d(tmp3,Lf->r/2.0);
      arb_pow(tmp4,tmp2,tmp3,prec);
      arb_mul(C_bit,tmp1,tmp4,prec);
      // half_log_N = 1/2 log(N)
      arb_log_ui(half_log_N,Lfu->N,prec);
      arb_mul_2exp_si(half_log_N,half_log_N,-1);
    }
  // u = x-1/2 log N
  arb_sub(u,x,half_log_N,prec);
  //printf("u=");arb_printd(u,10);printf("\n");
  // X = Pi r exp(2u/r)
  arb_div_ui(tmp1,u,Lf->r,prec);
  arb_mul_2exp_si(tmp1,tmp1,1);
  arb_exp(tmp2,tmp1,prec);
  arb_mul_ui(tmp1,tmp2,Lf->r,prec);
  arb_mul(X,tmp1,pi,prec);
  //printf("X=");arb_printd(X,10);printf("\n");
  // XM2r=XM^(2/r)
  arb_set_d(tmp1,2.0);
  arb_div_ui(tmp2,tmp1,Lf->r,prec);
  arb_set_ui(arb_M,M);
  arb_pow(tmp1,arb_M,tmp2,prec); // M^(2/r)
  arb_mul(XM2r,tmp1,X,prec); // XM^(2/r)
  arb_set_d(tmp1,Lf->r/2.0);
  arb_sub(tmp2,XM2r,tmp1,prec);
  if(!arb_is_positive(tmp2))
    {
      printf("Failed XM^(r/2)>r/2 test in Merror. Exiting.\n");
      return(false);
    }
  arb_mul(tmp2,tmp1,Lf->c,prec);
  arb_sub_ui(tmp1,tmp2,1,prec);
  arb_sub(tmp2,XM2r,tmp1,prec);
  if(!arb_is_positive(tmp1))
    {
      printf("Failed XM^(r/2)>cr/2-1 test in Merror. Exiting.\n");
      return(false);
    }
  
  //printf("X M^(2/r)=");arb_printd(XM2r,10);printf("\n");
  // Mc=M^c-1
  arb_sub_ui(tmp1,Lf->c,1,prec);
  arb_pow(Mc,arb_M,tmp1,prec);
  //printf("M^{c-1}=");arb_printd(Mc,10);printf("\n");
  arb_mul(tmp1,Lf->nu,u,prec);
  arb_sub(tmp3,tmp1,XM2r,prec); // nu*u-XM^(2/r)
  arb_exp(tmp1,tmp3,prec);
  arb_mul(tmp3,tmp1,Mc,prec); // M^{c-1} exp(nu*u-XM^(2/r))
  arb_add_ui(tmp1,XM2r,1,prec);
  arb_mul_2exp_si(tmp1,tmp1,1); // 2XM^(r/2)+2
  arb_mul_ui(tmp2,Lf->c,Lf->r,prec); // cr
  arb_sub(tmp4,tmp1,tmp2,prec);
  arb_mul_ui(tmp1,arb_M,Lf->r,prec);
  arb_div(tmp2,tmp1,tmp4,prec);
  arb_add_ui(tmp1,tmp2,1,prec);
  arb_mul(tmp2,tmp1,C_bit,prec);
  arb_mul(res,tmp2,tmp3,prec);
  //printf("before product=");arb_printd(res,10);printf("\n");
  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_mul_ui(tmp2,Lf->nus[j],Lf->r,prec);
      arb_div(tmp1,tmp2,XM2r,prec);
      arb_add_ui(tmp2,tmp1,1,prec);
      arb_pow(tmp1,tmp2,Lf->nus[j],prec);
      arb_mul(res,res,tmp1,prec);
    }
  return(true);
}

bool F_hat_twiddle_error(arb_t res, arb_t x, L_family_t *Lf, L_func_t *Lfu, L_comp_t *Lc, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3,u,X,half_log_N,pi,twoXbyr;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(u);
      arb_init(X);
      arb_init(half_log_N);
      arb_init(pi);
      arb_init(twoXbyr);
      // pi = pi
      arb_const_pi(pi,prec);
      // C_bit = Cr2^(r/2-1)
      // half_log_N = 1/2 log(N)
      arb_log_ui(half_log_N,Lfu->N,prec);
      arb_mul_2exp_si(half_log_N,half_log_N,-1);
    }
  // u = x-1/2 log N
  arb_sub(u,x,half_log_N,prec);
  //printf("u=");arb_printd(u,10);printf("\n");
  // X = Pi r exp(2u/r)
  arb_div_ui(tmp1,u,Lf->r,prec);
  arb_mul_2exp_si(tmp1,tmp1,1);
  arb_exp(tmp2,tmp1,prec);
  arb_mul_ui(tmp1,tmp2,Lf->r,prec);
  arb_mul(X,tmp1,pi,prec);
  // check X>r/2
  arb_set_d(tmp1,Lf->r/2.0);
  arb_sub(tmp2,X,tmp1,prec);
  if(!arb_is_positive(tmp2))
    {
      printf("Fhat twiddle error failed X>r/2 test. Exiting.\n");
      return(false);
    }
  //printf("X=");arb_printd(X,10);printf("\n");
  // 2X/r
  arb_div_ui(twoXbyr,X,Lf->r,prec);
  arb_mul_2exp_si(twoXbyr,twoXbyr,1);
  arb_zeta(tmp2,twoXbyr,prec);
  arb_pow_ui(tmp1,tmp2,Lf->r,prec);
  arb_set_d(tmp2,0.5);
  arb_sub(tmp3,tmp2,twoXbyr,prec);
  arb_mul(tmp2,tmp3,Lc->arb_A,prec);
  arb_mul(tmp3,tmp2,pi,prec);
  arb_mul_2exp_si(tmp3,tmp3,1); // 2piA(1/2-2X/r)
  arb_exp(tmp2,tmp3,prec);
  arb_sub_ui(tmp3,tmp2,1,prec);
  arb_neg(tmp3,tmp3);
  arb_div(tmp2,tmp1,tmp3,prec);
  arb_mul(tmp1,u,Lf->nu,prec);
  arb_sub(tmp3,tmp1,X,prec);
  arb_exp(tmp1,tmp3,prec);
  arb_mul(tmp3,tmp1,tmp2,prec);
  arb_sqrt_ui(tmp1,2,prec);
  arb_pow_ui(tmp2,tmp1,Lf->r,prec);
  arb_mul(res,tmp2,tmp3,prec);
  for(uint64_t j=0;j<Lf->r;j++)
    {
      arb_mul_ui(tmp1,Lf->nus[j],Lf->r,prec);
      arb_div(tmp2,tmp1,X,prec);
      arb_add_ui(tmp1,tmp2,1,prec);
      arb_pow(tmp2,tmp1,Lf->nus[j],prec);
      arb_mul(res,res,tmp2,prec);
    }
  return(true);
}

// thanks Bober
void acb_reasonable_sqrt(acb_t out, acb_t in, slong prec) {
    if(arb_is_negative(acb_realref(in)) && arb_contains_zero(acb_imagref(in))) {
        acb_neg(in, in);
        acb_sqrt(out, in, prec);
        acb_mul_onei(out, out);
        acb_neg(in, in);
    }
    else {
        acb_sqrt(out, in, prec);
    }
}

void fix_epsilon(acb_t res, acb_t x, acb_t y,uint64_t prec)
{
  static acb_t th1,th2;
  arb_t ep;
  static bool init=false;
  if(!init)
    {
      init=true;
      acb_init(th1);
      acb_init(th2);
      arb_init(ep);
    }
  acb_mul(th1,x,y,prec);
  acb_reasonable_sqrt(th2,th1,prec);
  //printf("sqrt=");acb_printd(th2,10);printf("\n");
  acb_arg(ep,th2,prec);
  arb_sin_cos(acb_imagref(res),acb_realref(res),ep,prec);
}


// when this is called, we have convolved with all
// a_m from M0 to M and then done m=1->M0 simplistically
// on all n from 0 to hi_i
// thus n=0 has all m<=M
// n=n has all m<=max(floor(exp((hi_i-n+0.5)*2*Pi/B)*sqrt(N)),M0)

// given n, what is the last a_m have we summed into res[n]
// we need to add the error for sum m>=m+1
// returns 0 if no coefficients have been added in
// then use F_hat_twiddle bound
uint64_t inv_m(uint64_t hi_i, uint64_t n, uint64_t M0, double one_over_B, double dc)
{
  return(exp(((double) hi_i-(double) n+0.5)*2.0*M_PI*one_over_B)*dc);
}

bool do_pre_iFFT_errors(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, L_error_t *Le, uint64_t prec)
{
  static arb_t err,x,fhattwiddle,pi,two_pi_A,tmp,th,th1,th2;
  static acb_t ctmp1,ctmp2;
  static bool init=false;
  if(!init)
    {
      arb_init(th);
      arb_init(th2);
      arb_init(th2);
      arb_init(fhattwiddle);
      arb_init(err);arb_init(x);arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(two_pi_A);
      arb_mul(two_pi_A,pi,Lc->arb_A,prec);
      arb_mul_2exp_si(two_pi_A,two_pi_A,1);
      arb_init(tmp);
      acb_init(ctmp1);
      acb_init(ctmp2);

      init=true;
    }
  uint64_t i=0;
  arb_zero(x);
  for(i=0;i<=Lc->N/2;i++)
    {
      uint64_t M=inv_m(Lf->hi_i,i,Lc->M0,Lf->one_over_B,Lfu->dc);
      //printf("i=%lu => M=%lu ",i,M);
      if(M==0)
	break;
      if(!M_error(err,x,M+1,Lf,Lfu,prec))
	return(false);
      //printf("M Error for n=%lu is ",i);
      //arb_printd(err,10);
      //printf("\n");
      arb_add_error(acb_realref(Lc->res[i]),err);
      arb_add_error(acb_imagref(Lc->res[i]),err);
      arb_add(x,x,Lf->two_pi_by_B,prec);
    }

  // error sum k\neq 0 F_hat(x+2\pi k A)
  arb_sub(err,two_pi_A,x,prec);
  if(!F_hat_twiddle_error(fhattwiddle,err,Lf,Lfu,Lc,prec))
    return(false);
  arb_mul_2exp_si(fhattwiddle,fhattwiddle,1);
  //printf("F_hat_twiddle error at n=%lu = ",i);
  //arb_printd(fhattwiddle,10);
  //printf("\n");
  arb_mul(err,Le->eq59,Lfu->sum_ans,prec);
  printf("Adding eq 5-9 error = ");arb_printd(err,10);printf("\n");

  for(uint64_t j=0;j<i;j++)
    {
      arb_add_error(acb_realref(Lc->res[i]),fhattwiddle);
      arb_add_error(acb_imagref(Lc->res[i]),fhattwiddle);
      arb_add_error(acb_realref(Lc->res[i]),err);
      arb_add_error(acb_imagref(Lc->res[i]),err);
    }

  // figure out epsilon
  
  if(acb_contains_zero(Lc->res[0]))
    {
      //printf("F_hat[1]=");acb_printd(Lc->res[1],10);printf("\n");
      //printf("F_hat[-1]=");acb_printd(Lc->res[Lc->N-1],10);printf("\n");
      //printf("Using F_hat[1] and F_hat[-1] to fix epsilon.\n");
      fix_epsilon(Lfu->epsilon,Lc->res[1],Lc->res[Lc->N-1],prec);
    }
  else
    {  
      //printf("F_hat[0]=");acb_printd(Lc->res[1],10);printf("\n");
      if((arb_is_negative(acb_realref(Lc->res[0])))&&(arb_contains_zero(acb_imagref(Lc->res[0]))))
	{
	  acb_neg(ctmp1,Lc->res[0]);
	  acb_arg(th,ctmp1,prec);
	  arb_sin_cos(th1,th2,th,prec);
	  arb_set(acb_realref(Lfu->epsilon),th2);
	  arb_set(acb_imagref(Lfu->epsilon),th1);
	}
      else
	{
	  acb_arg(th,Lc->res[0],prec);
	  //acb_arg(th2,Lc->res[Lc->N/2],prec);
	  //arb_set(th,th1); //arb_intersect(th,th1,th2,prec);
	  arb_neg(th,th);
	  arb_sin_cos(th1,th2,th,prec);
	  arb_set(acb_realref(Lfu->epsilon),th2);
	  arb_set(acb_imagref(Lfu->epsilon),th1);
	}
    }
  printf("epsilon set to ");
  acb_printd(Lfu->epsilon,10);
  printf("\n");
  arb_sub_ui(th,acb_realref(Lfu->epsilon),1,prec);
  Lfu->self_dual_p=arb_contains_zero(th);
  if(Lfu->self_dual_p)
    printf("Looks like this L-function is self dual. Assuming that.\n");
  else
    printf("Looks like this L-function is not self dual. Assuming that.\n");
  
  for(uint64_t n=0;n<i;n++)
    acb_mul(Lc->res[n],Lc->res[n],Lfu->epsilon,prec);
  
  // the rest of the vector we just approximate with F_hat_twiddle error

  // error sum F_hat(x+2\pi k A)
  if(!F_hat_twiddle_error(fhattwiddle,x,Lf,Lfu,Lc,prec))
    return(false);
  // we have sum k\geq 0 F_hat(x+2\pi k A)
  // since x<\pi A we can just double this
  arb_mul_2exp_si(fhattwiddle,fhattwiddle,1);
  //printf("F_hat_twiddle error beyond n=%lu = ",i);
  //arb_printd(fhattwiddle,10);
  //printf("\n");
  for(;i<=Lc->NN/2;i++)
    {
      acb_zero(Lc->res[i]);
      arb_add_error(acb_realref(Lc->res[i]),fhattwiddle);
      arb_add_error(acb_imagref(Lc->res[i]),fhattwiddle);
    }
  for(uint64_t n=Lc->NN/2+1;n<Lc->NN;n++)
    acb_conj(Lc->res[n],Lc->res[Lc->NN-n]);
  /*
  for(uint64_t n=0;n<Lc->N/2;n++)
    {
      printf("Pre iFFT res[%lu] = ",n);
      acb_printd(Lc->res[n],10);
      printf("\n");
    }
  */
  return(true);
}

