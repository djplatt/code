//
// hurwitz.c
//
// compute zeta(1/2+it_n+k,M+j/N_ROWS)
// where j goes from N_ROWS to 1
// k goes from 0 to N_COLUMNS-1
//
// t_0=t_start
// t_{n+1}=t_n+DEL_T
//
//
// Taylor error reconstructing zeta(1/2+it,alpha) will be < 1.521e-96
// Lemma 7.2.3 My Thesis
//
#include "inttypes.h"
#include "fmpr.h"
#include "fmprb.h"
#include "fmpcb.h"
#define N_COLUMNS (50) // should give us 316 bits of absolute precision
#define N_ROWS (512)
#define DEL_T ((double) 5.0/128.0)
#define PREC (300)
#define M (2)

int main(int argc, char **argv)
{
  if(argc!=4)
    {
      printf("Usage:- %s <t start> <t end> <outfile>.\n",argv[0]);
      exit(0);
    }
  double t_start=atof(argv[1]);
  double t_end=atof(argv[2]);
  FILE *outfile=fopen(argv[3],"wb");
  if(!outfile)
    {
      printf("Failed to open file %s for binary output. Exiting.\n",argv[3]);
      exit(0);
    }
  uint64_t n_COLUMNS=N_COLUMNS;
  uint64_t n_ROWS=N_ROWS;
  uint64_t mM=M;
  double dEL_T=DEL_T;

  if(!fwrite(&n_COLUMNS,sizeof(uint64_t),1,outfile)||
     !fwrite(&n_ROWS,sizeof(uint64_t),1,outfile)||
     !fwrite(&mM,sizeof(uint64_t),1,outfile)||
     !fwrite(&t_start,sizeof(double),1,outfile)||
     !fwrite(&t_end,sizeof(double),1,outfile)||
     !fwrite(&dEL_T,sizeof(double),1,outfile))
    {
      printf("Error writing to %s, Exiting.\n",argv[3]);
      exit(0);
    }
  uint64_t i,row,column,s;
  fmprb_t r[N_COLUMNS],as[N_ROWS+1],temp;
  fmprb_init(temp);
  for(i=0;i<N_COLUMNS;i++)
    {
      fmprb_init(r[i]);
      fmprb_set_ui(temp,1+(i<<1));
      fmprb_div_ui(r[i],temp,2,32);
      //printf("r[%lu]=",i);fmprb_printd(r[i],10);printf("\n");
    }
  for(i=0;i<=N_ROWS;i++)
    {
      fmprb_init(as[i]);
      fmprb_set_ui(as[i],i);
      fmprb_div_ui(temp,as[i],N_ROWS,32);
      fmprb_add_ui(as[i],temp,M,32);
    }

  fmpr_t t;
  double td;
  fmpr_init(t);

  fmpz_t left,right,exp;
  fmpz_init(left);
  fmpz_init(right);
  fmpz_init(exp);

  mpz_t m;
  mpz_init(m);

  fmprb_t pi4,tpi4;
  fmprb_init(pi4);
  fmprb_init(tpi4);
  fmprb_const_pi(tpi4,PREC);
  fmprb_mul_2exp_si(pi4,tpi4,-2);

  fmpcb_t a,z,w,w1;
  fmpcb_init(a);
  fmpcb_init(z);
  fmpcb_init(w);
  fmpcb_init(w1);
  fmprb_set_ui(fmpcb_imagref(a),0);

  for(td=t_start;td<t_end;td+=DEL_T)
    {
      fmpr_set_d(t,td);
      fmprb_mul_fmpr(tpi4,pi4,t,PREC);
      fmprb_set_fmpr(fmpcb_imagref(z),t);
      fmprb_set(fmpcb_realref(z),r[0]);
      fmpcb_mul_2exp_si(w,z,-1);
      //fmpcb_printd(w,20);printf("\n");
      fmpcb_lgamma(w1,w,PREC);
      fmpcb_add_fmprb(w,w1,tpi4,PREC);
      //fmpcb_printd(w,20);printf("\n");
      fmprb_get_interval_fmpz_2exp(left,right,exp,fmpcb_realref(w));
      fmpz_get_mpz(m,left);
      mpz_out_raw(outfile,m);
      fmpz_get_mpz(m,right);
      mpz_out_raw(outfile,m);
      fmpz_get_mpz(m,exp);
      mpz_out_raw(outfile,m);
      fmprb_get_interval_fmpz_2exp(left,right,exp,fmpcb_imagref(w));
      fmpz_get_mpz(m,left);
      mpz_out_raw(outfile,m);
      fmpz_get_mpz(m,right);
      mpz_out_raw(outfile,m);
      fmpz_get_mpz(m,exp);
      mpz_out_raw(outfile,m);
      fmpcb_add_ui(w1,z,1,PREC);
      fmpcb_mul_2exp_si(w,w1,-1);
      fmpcb_lgamma(w1,w,PREC);
      fmpcb_add_fmprb(w,w1,tpi4,PREC);
      fmprb_get_interval_fmpz_2exp(left,right,exp,fmpcb_realref(w));
      fmpz_get_mpz(m,left);
      mpz_out_raw(outfile,m);
      fmpz_get_mpz(m,right);
      mpz_out_raw(outfile,m);
      fmpz_get_mpz(m,exp);
      mpz_out_raw(outfile,m);
      fmprb_get_interval_fmpz_2exp(left,right,exp,fmpcb_imagref(w));
      fmpz_get_mpz(m,left);
      mpz_out_raw(outfile,m);
      fmpz_get_mpz(m,right);
      mpz_out_raw(outfile,m);
      fmpz_get_mpz(m,exp);
      mpz_out_raw(outfile,m);
      

      for(row=0;row<=N_ROWS;row++)
	{
	  fmprb_set(fmpcb_realref(a),as[row]);
	  //printf("r[0]=");fmprb_printd(r[0],10);printf("\n");
	  for(column=0;column<N_COLUMNS;column++)
	    {
	      fmprb_set(fmpcb_realref(z),r[column]);
	      //printf("Doing Hurwitz on ");fmpcb_printd(z,20);
	      //printf("\nand ");fmpcb_printd(a,20);	      
	      fmpcb_hurwitz_zeta(w,z,a,PREC);
	      //printf("\ngot ");fmpcb_printd(w,20);printf("\n");	      
	      fmprb_get_interval_fmpz_2exp(left,right,exp,fmpcb_realref(w));
	      fmpz_get_mpz(m,left);
	      mpz_out_raw(outfile,m);
	      fmpz_get_mpz(m,right);
	      mpz_out_raw(outfile,m);
	      fmpz_get_mpz(m,exp);
	      mpz_out_raw(outfile,m);
	      fmprb_get_interval_fmpz_2exp(left,right,exp,fmpcb_imagref(w));
	      fmpz_get_mpz(m,left);
	      mpz_out_raw(outfile,m);
	      fmpz_get_mpz(m,right);
	      mpz_out_raw(outfile,m);
	      fmpz_get_mpz(m,exp);
	      mpz_out_raw(outfile,m);
	    }
	}
    }
}
