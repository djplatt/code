//#define BEST_POSS

#ifdef BEST_POSS
arb_set_ui(tmp1,927436);
arb_div_ui(A6,tmp1,100000,prec);
#else
arb_set_ui(tmp1,872606);
arb_div_ui(A6,tmp1,100000,prec);
#endif
arb_set_si(tmp1,-813199);
arb_div_ui(A7,tmp1,100000,prec);

// alpha =2/5
arb_set_ui(tmp1,2);
arb_div_ui(tmp1,tmp1,5,prec);
arb_sub(tmp1,tmp1,alpha,prec);
if(!arb_contains_zero(tmp1))
  {
    printf("Looks like you've compiled with wrong constatnts for this alpha. Exiting.\n");
    exit(0);
  }

//c(alpha)=0.6877 (at x->6)
arb_set_ui(tmp1,1458001753);
arb_div_ui(A8,tmp1,100000,prec);

arb_set_ui(tmp1,2763359);
arb_div_ui(A9,tmp1,100000,prec);


#ifdef BEST_POSS
arb_zero(A7);
arb_zero(A8);
arb_zero(A9);
#endif

