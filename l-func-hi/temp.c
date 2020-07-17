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
#define N_COLUMNS (40) // should give us 318 bits of precision
#define N_ROWS (4096)
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
  uint64_t i,row,column,s;
  fmprb_t r[N_COLUMNS],as[N_ROWS],temp;
  fmprb_init(temp);
  for(i=0;i<N_COLUMNS;i++)
    {
      fmprb_init(r[i]);
      fmprb_set_ui(temp,1+(i<<1));
      fmprb_div_ui(r[i],temp,2,32);
    }
  for(i=0;i<N_ROWS;i++)
    {
      fmprb_init(as[i]);
      fmprb_set_ui(as[i],N_ROWS-i);
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

  fmpcb_t a,z,w;
  fmpcb_init(a);
  fmpcb_init(z);
  fmpcb_init(w);
  fmprb_set_ui(fmpcb_imagref(a),0);

  for(td=t_start;td<t_end;td+=DEL_T)
    {
      fmpr_set_d(t,td);
      fmprb_set_fmpr(fmpcb_imagref(z),t);
      for(row=0;row<N_ROWS;row++)
	{
	  fmprb_set(fmpcb_realref(a),as[row]);
	  for(column=0;column<N_COLUMNS;column++)
	    {
	      fmprb_set(fmpcb_realref(z),r[column]);
	      fmpcb_hurwitz_zeta(w,z,a,PREC);
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
