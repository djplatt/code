#include "stdlib.h"
#include "stdio.h"
#include "gmp.h"
#include "mpfr.h"

#define T (100) // powers of t in expansion of f
#define K (10) // number of error terms in RS

mpfr_t qd_tmp;

void init_out_qd()
{
  mpfr_init(qd_tmp);
}

void out_qd(mpfr_ptr x, FILE *outfile)
{
  long unsigned int e;
  double d;

  d=mpfr_get_d(x,GMP_RNDZ);
  fwrite(&d,sizeof(double),1,outfile);
  //printf("d=%e\n",d);
  mpfr_sub_d(qd_tmp,x,d,GMP_RNDN);

  d=mpfr_get_d(qd_tmp,GMP_RNDZ);
  fwrite(&d,sizeof(double),1,outfile);
  //printf("d=%e\n",d);
  mpfr_sub_d(qd_tmp,qd_tmp,d,GMP_RNDN);

  d=mpfr_get_d(qd_tmp,GMP_RNDZ);
  fwrite(&d,sizeof(double),1,outfile);
  //printf("d=%e\n",d);
  mpfr_sub_d(qd_tmp,qd_tmp,d,GMP_RNDN);

  d=mpfr_get_d(qd_tmp,GMP_RNDZ);
  fwrite(&d,sizeof(double),1,outfile);
  //printf("d=%e\n",d);
  //mpfr_sub_d(x,x,d,GMP_RNDN);

}

int main(int argc,char **argv)
{
  // coefficients of the Taylor expansion
  // of f(t)= cos(pi/2(t^2+3/4))/cos(pi t)
  // and first 3K of its derivatives
  // ts[a][b] = coeff of t^b in f^(a)(t) 
  mpfr_t ts[3*K+1][T+1]; 
  
  // coefficients of RS terms
  // cs[a][b]= coefficient of f^(b)(2pi/t)^(a/2) in a'th RS term
  // so cs[0][0]=1, cs[0][*]=0
  // cs[1][3]=-1/96/pi^2, cs[1][*]=0
  // cs[2][6]=1/18432/pi^4, cs[2][2]=1/64/pi^2, cs[2][*]=0
  // ...
  mpfr_t cs[K+1][3*K+1];

  mpfr_t pi_pows[2*K+1];

  unsigned long int t,k;
  if(argc!=5)
    {
      printf("Usage:- rs_coeff <prec> <ts_file> <cs_file> <out_file>\n");
      exit(0);
    }
  mpfr_set_default_prec(atol(argv[1]));
  //printf("Initialising ts.\n");
  for(k=0;k<=3*K;k++)
    for(t=0;t<=T;t++)
      mpfr_init(ts[k][t]);
  for(k=0;k<=2*K;k+=2)
    mpfr_init(pi_pows[k]);
  mpfr_set_ui(pi_pows[0],1,GMP_RNDN);
  mpfr_const_pi(pi_pows[2],GMP_RNDN);
  mpfr_mul(pi_pows[2],pi_pows[2],pi_pows[2],GMP_RNDN);
  for(k=4;k<=2*K;k+=2)
    mpfr_mul(pi_pows[k],pi_pows[k-2],pi_pows[2],GMP_RNDN);
  FILE *tfile=fopen(argv[2],"r");
  FILE *cfile=fopen(argv[3],"r");
  FILE *outfile=fopen(argv[4],"wb");
  if((!tfile)||(!cfile)||(!outfile))
    {
      if(!tfile)
	printf("main: error opening %s. Exiting.\n",argv[2]);
      if(!cfile)
	printf("main: error opening %s. Exiting.\n",argv[3]);
      if(!outfile)
	printf("main: error opening %s. Exiting.\n",argv[4]);
      exit(0);
    }
  //read the Taylor coefficients of f(t)=cos(pi/2*(t^2+3/4))/code(pi*t)
  // at t=0
  // we want g(t)=cos(2pi(t^2-t-1/16))/cos(2pi t)
  // so multiply coeff of t^n by 2^n
  //odd ones are zero
  mpfr_t pow4;
  mpfr_init(pow4);
  mpfr_set_ui(pow4,1,GMP_RNDN);
  for(t=0;t<=T;t+=2)
    {
      if(!mpfr_inp_str(ts[0][t],tfile,10,GMP_RNDN))
	{
	  printf("main: error reading coefficient %lu from tfile. Exiting.\n",t);
	  exit(0);
	}
      mpfr_mul(ts[0][t],ts[0][t],pow4,GMP_RNDN);
      mpfr_mul_ui(pow4,pow4,4,GMP_RNDN);
    }
  mpfr_clear(pow4);
  for(t=1;t<=T;t+=2)
    mpfr_set_ui(ts[0][t],0,GMP_RNDN);
  // now differentiate 30 times
  for(k=1;k<=K*3;k++)
    {
      mpfr_set_ui(ts[k][T],0,GMP_RNDN);
      for(t=0;t<T;t++)
	mpfr_mul_ui(ts[k][t],ts[k-1][t+1],t+1,GMP_RNDN);
    }
  fclose(tfile);

  long int d;
  long int den;
  unsigned long int num,pis,diff;
  for(k=0;k<=K;k++)
    {
      d=k*3; // maximum differential we will see
      while(d>=0)
	{
	  fscanf(cfile,"%ld %lu %lu %lu",&den,&num,&pis,&diff);
	  //printf("Read %ld %lu %lu %lu\n",den,num,pis,diff);
	  mpfr_init(cs[k][diff]);
	  if(den<0)
	    {
	      mpfr_set_ui(cs[k][diff],-den,GMP_RNDN);
	      mpfr_neg(cs[k][diff],cs[k][diff],GMP_RNDN);
	    }
	  else
	    mpfr_set_ui(cs[k][diff],den,GMP_RNDN);
	  mpfr_div_ui(cs[k][diff],cs[k][diff],num,GMP_RNDN);
	  mpfr_div(cs[k][diff],cs[k][diff],pi_pows[pis],GMP_RNDN);
	  //printf("cs[%lu][%lu] set to ",k,diff);mpfr_out_str(stdout,10,0,cs[k][diff],GMP_RNDN);printf("\n");
	  d-=4;
	}
    }

  mpfr_t out_ts[T+1];
  mpfr_t tmp;
  mpfr_init(tmp);
  for(t=0;t<=T;t++)
    mpfr_init(out_ts[t]);
  init_out_qd();
  k=K;
  fwrite(&k,sizeof(unsigned long int),1,outfile);
  t=T;
  fwrite(&t,sizeof(unsigned long int),1,outfile);
  long int m1=-1;
  for(k=0;k<=K;k++)
    {
      fwrite(&k,sizeof(unsigned long int),1,outfile);
      for(t=0;t<=T;t++)
	mpfr_set_ui(out_ts[t],0,GMP_RNDN);
      for(d=3*k;d>=0;d-=4)
	for(t=0;t<=T;t++)
	  {
	    mpfr_mul(tmp,cs[k][d],ts[d][t],GMP_RNDN);
	    mpfr_add(out_ts[t],out_ts[t],tmp,GMP_RNDN);
	  }
      for(t=0;t<=T;t++)
	{
	  if(!mpfr_zero_p(out_ts[t]))
	    {
	      //printf("Coeff for t^%lu =",t);
	      //mpfr_out_str(stdout,10,0,out_ts[t],GMP_RNDN);
	      //printf("\n");
	      fwrite(&t,sizeof(unsigned long int),1,outfile);
	      out_qd(out_ts[t],outfile);
	    }
	}
      fwrite(&m1,sizeof(long int),1,outfile); // end of t's for this k
    }
  fclose(outfile);
  return(0);
}
