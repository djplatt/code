#include <inttypes.h>
#include <stdlib.h>
#include <stdbool.h>
#include <flint/fmpz.h>
#include <flint/acb_dirichlet.h>


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
      arb_log(log_pi,pi,prec);
      arb_inv(pi,pi,prec);
      acb_init(ctmp1);acb_init(ctmp2);acb_init(ctmp3);
    }
  //printf("gamma = ");arb_printd(gamma,20);printf("\n");
  arb_inv(err,gamma,prec);
  arb_mul_2exp_si(err,err,-2); //1/4t
  arb_mul_2exp_si(acb_imagref(ctmp1),gamma,-1); // it/2
  arb_set_d(acb_realref(ctmp1),-0.25); // -1/4+it/2
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
  
  arb_set_ui(rtmp1,z-1);
  arb_sub(rtmp3,rtmp1,rtmp2,prec); // S(T) just after the zero
  printf("S(t) just after zero %ld is ",z);arb_printd(rtmp3,20);printf("\n");
  arb_sub_ui(rtmp1,rtmp3,1,prec); // S(T) just before the zero
  printf("S(t) just before zero %ld is ",z);arb_printd(rtmp1,20);printf("\n");
  printf("T = ");arb_printd(gamma,20);printf("\n");
}


int main(int argc, char**argv)
{

  printf("Command line:-");
  for(uint64_t i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=3)
    {
      printf("Usage:- %s <list of zeros> <prec>\n",argv[0]);
      return 0;
    }

  FILE *infile=fopen(argv[1],"r");
  if(!infile)
    {
      perror("Error opening file:-");
      return 0;
    }
  int64_t prec=atol(argv[2]);

  arb_t t;
  arb_init(t);
  int64_t z;
  fmpz_t fz;
  fmpz_init(fz);
  while(fscanf(infile,"%ld",&z)==1)
    {
      fmpz_set_si(fz,z);
      acb_dirichlet_hardy_z_zero(t,fz,prec);
      St(t,z,prec);
    }
  
  return 0;
}
