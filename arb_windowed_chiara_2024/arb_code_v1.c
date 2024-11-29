/*
Look for extremal values of S(t)
This is for zeros beyond the database of zeros to ~ 3e10
so we use acb_dirichlet_platt_hardy_z_zeros to find them

Copyright DJ Platt 2024

*/
#include <inttypes.h>
#include <stdlib.h>
#include <stdbool.h>
#include <flint/fmpz.h>
#include <flint/acb_dirichlet.h>

arb_t big_st,little_st; // largest and smallest S(t) found
int64_t big_zero,little_zero; // the respective zero numbers

// Estimate S(t) just before and just after the z'th zero at 1/2+i gamma
// S(t)=N(T)-1/pi Phi(t)-1
// Phi(t)=-t log pi/2+Im log_gamma(1/4+it/2)
// Im log_gamma(z)=Im[(z-1/2)log(z/e)]-Theta(1/(8Im z)) Booker 2006 
void St(arb_t gamma, int64_t z, int64_t prec)
{
  static acb_t ctmp1,ctmp2,ctmp3;
  static arb_t pi,log_pi,half,rtmp1,rtmp2,rtmp3,t_2_log_pi,err;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(pi);arb_init(log_pi);arb_init(half);
      arb_init(rtmp1);arb_init(rtmp2);arb_init(rtmp3);
      arb_init(t_2_log_pi);arb_init(err);
      arb_set_d(half,0.5);
      arb_const_pi(pi,prec);
      arb_log(log_pi,pi,prec); // log(pi)
      arb_inv(pi,pi,prec); // 1/pi
      acb_init(ctmp1);acb_init(ctmp2);acb_init(ctmp3);
    }
  //printf("gamma = ");arb_printd(gamma,20);printf("\n");
  arb_inv(err,gamma,prec);
  arb_mul_2exp_si(err,err,-2); //1/(4 gamma)
  arb_mul_2exp_si(acb_imagref(ctmp1),gamma,-1); // i gamma/2
  arb_set_d(acb_realref(ctmp1),-0.25); // -1/4+i gamma/2
  arb_mul(t_2_log_pi,log_pi,acb_imagref(ctmp1),prec); // t/2 log pi
  arb_set(acb_imagref(ctmp2),acb_imagref(ctmp1)); // it/2
  arb_set_d(acb_realref(ctmp2),0.25); // 1/4+it/2
  acb_log(ctmp3,ctmp2,prec); // log(1/4+it/2)
  arb_sub_ui(acb_realref(ctmp3),acb_realref(ctmp3),1,prec); // log(1/4+it/2)-1
  acb_mul(ctmp2,ctmp1,ctmp3,prec); // (-1/4+it/2)(log(1/4+it/2)-1)
  //printf("()*()=");acb_printd(ctmp2,20);printf("\n");
  arb_sub(rtmp1,acb_imagref(ctmp2),t_2_log_pi,prec); // t/2 log pi -Im[()()]
  arb_add_error(rtmp1,err);
  arb_mul(rtmp2,rtmp1,pi,prec); // -1/pi[Im[()()-t/2 log pi +/-1/4t]]

  //printf("Argument bit = ");arb_printd(rtmp2,20);printf("\n");
  
  arb_set_ui(rtmp1,z-1); // N(t)-1
  arb_sub(rtmp3,rtmp1,rtmp2,prec); // S(T) just after the zero
  arb_sub(rtmp1,big_st,rtmp3,prec);
  if(!arb_is_positive(rtmp1))
    {
      if(arb_is_negative(rtmp1)) // S(T) is bigger than before
	{
	  arb_set(big_st,rtmp3);
	  big_zero=z;      
	  printf("New largest at %ld ",z);arb_printd(big_st,20);printf("\n");
	  fflush(stdout);
	}
      else
	{
	  printf("Indeterminate comparison just after zero = %ld. Skipping it.\n",z);
	  fflush(stdout);
	}
    }
  
  arb_sub_ui(rtmp1,rtmp3,1,prec); // S(T) just before the zero
  arb_sub(rtmp2,little_st,rtmp1,prec);
  if(!arb_is_negative(rtmp2))
    {
      if(arb_is_positive(rtmp2))
	{
	  arb_set(little_st,rtmp1);
	  little_zero=z;      
	  printf("New smallest at %ld ",z);arb_printd(little_st,20);printf("\n");
	  fflush(stdout);
	}
      else
	{
	  printf("Indeterminate comparison just before zero = %ld. Skipping it.\n",z);
	  fflush(stdout);
	}
    }
}


int main(int argc, char**argv)
{
  printf("Command line:-");
  for(uint64_t i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=5)
    {
      printf("Usage:- %s <zeros per iteration> <starting zero> <# iterations> <prec>\n",argv[0]);
      return 0;
    }
  
  int64_t n_zeros=atol(argv[1]);
  int64_t istart_zero=atol(argv[2]);
  fmpz_t start_zero,this_start;
  fmpz_init(start_zero);
  int64_t num_its=atol(argv[3]);
  int64_t prec=atol(argv[4]);

  arb_init(big_st);arb_init(little_st);
  
  arb_t *zeros;
  zeros=(arb_t *)malloc(sizeof(arb_t)*n_zeros);
  for(uint64_t i=0;i<n_zeros;i++)
    arb_init(zeros[i]);

  //acb_dirichlet_hardy_z_zeros((arb_ptr) zeros,start_zero,n_zeros,prec);

  int64_t st=istart_zero;
  printf("First zero processed was %ld\n",st);
  for(int64_t it=0;it<num_its;it++)
    {
      fmpz_set_si(start_zero,st);
      int64_t this_prec=prec;
      while(acb_dirichlet_platt_hardy_z_zeros((arb_ptr) zeros,start_zero,n_zeros,this_prec)!=n_zeros) // some zeros missed, so increase precision
	{
	  this_prec+=40; // 40 more bits (~ 12 dp)
	  if(this_prec>200) // that's enough, let's quit
	    {
	      printf("Failed at zero %ld even at 200 bits. Exiting.\n",st);
	      return 0;
	    }
	  printf("Increasing precision at zero %ld. Now %ld\n",st,this_prec);
	}
      for(int64_t z=0;z<n_zeros;z++) // check S(t) at each zero
	St(zeros[z],z+st,prec);
      st+=n_zeros;
    }

  printf("Last zero processed was %ld\n",st-1);
  
  return 0;
}
