//
// File: pintz_no_grh.c
//
// DJ Platt, Heilbronn Institute, 9 October 2013
//
// Implement Pintz and Ruzsa approach to Linnik's constants
//
// As deciphered by T Trudgian, ANU
//
// gcc pintz.c -I /home/dave/flint2-trunk/ -lgmp -lflint -larb -lpthread
//
// Uses the arb library for multi precision interval arithmetic
//
#include "inttypes.h"
#include "math.h"
#include "fmpr.h"
#include "fmprb.h"
#include "fmprb_mat.h"
//#include "fmpcb.h"

uint64_t L,k;

double lambdad,cd;

#define PREC ((uint64_t) 200)
uint64_t cnum,cden;

fmprb_t *sigmas,*rhos,*bs,*twopims,*costwopims,*sintwopims,
  *lamcostwopims,*explamcostwopims,pi,temp,temp1,temp2,c;
fmprb_t lambda,log_2;

fmprb_mat_t U,U1;

// set up memory for vectors
void calc_twopims()
{
  sigmas=(fmprb_t *)malloc(sizeof(fmprb_t)*(k+1));
  rhos=(fmprb_t *)malloc(sizeof(fmprb_t)*(k+1));
  bs=(fmprb_t *)malloc(sizeof(fmprb_t)*(k+1));
  twopims=(fmprb_t *)malloc(sizeof(fmprb_t)*(k+1));
  costwopims=(fmprb_t *)malloc(sizeof(fmprb_t)*(k+1));
  sintwopims=(fmprb_t *)malloc(sizeof(fmprb_t)*(k+1));
  lamcostwopims=(fmprb_t *)malloc(sizeof(fmprb_t)*(k+1));
  explamcostwopims=(fmprb_t *)malloc(sizeof(fmprb_t)*(k+1));
  
  if(!sigmas||!rhos||!bs||!twopims||!costwopims||!sintwopims||!lamcostwopims||!explamcostwopims)
    {
      printf("Error allocating memory. Exiting.\n");
      exit(0);
    }
  fmprb_mat_init(U,k,k);
  fmprb_mat_init(U1,k,k);
  fmprb_init(temp);
  fmprb_init(temp1);
  fmprb_init(temp2);
  fmprb_init(pi);
  fmprb_const_pi(pi,PREC);
  fmprb_init(log_2);
  fmprb_const_log2(log_2,PREC);

  fmprb_init(c);
  fmprb_set_ui(temp,cnum);
  fmprb_div_ui(temp1,temp,cden,PREC);
  fmprb_mul(c,temp1,log_2,PREC);

  uint64_t m;
  for(m=0;m<=k;m++)
    {
      fmprb_mul_ui(temp,pi,m<<1,PREC); // 2 Pi m
      fmprb_init(twopims[m]);
      fmprb_div_ui(twopims[m],temp,k+1,PREC); // (2 Pi m)/(k+1)
      fmprb_init(costwopims[m]);
      fmprb_init(sintwopims[m]);
      fmprb_sin_cos(sintwopims[m],costwopims[m],twopims[m],PREC);
      fmprb_init(explamcostwopims[m]);
      fmprb_init(lamcostwopims[m]);
    }
}

// compute sigma_j=1/(k+1)sum_{m=0}^k exp(lam cos(2pi m/(k+1))cos(2pi mj/(k+1))
void make_sigmas()
{
  calc_twopims();
  uint64_t j,m,mj;
  for(m=0;m<=k;m++)
    {
      fmprb_mul(lamcostwopims[m],costwopims[m],lambda,PREC);
      fmprb_exp(explamcostwopims[m],lamcostwopims[m],PREC);
    }
  for(j=0;j<=k;j++)
    {
      fmprb_set_ui(sigmas[j],0);
      for(m=0;m<=k;m++)
	{
	  mj=(m*j)%(k+1);
	  fmprb_mul(temp,explamcostwopims[m],costwopims[mj],PREC);
	  fmprb_add(temp1,sigmas[j],temp,PREC);
	  fmprb_set(sigmas[j],temp1);
	}
      fmprb_div_ui(temp,sigmas[j],k+1,PREC);
      fmprb_set(sigmas[j],temp);
    }
  /*
  for(m=0;m<=k;m++)
    {
      printf("sigma[%lu]=",m);
      fmprb_printd(sigmas[m],10);
      printf("\n");
    }
  */
}

void make_new_sigmas()
{
  uint64_t j,m,mj;
  for(m=0;m<=k;m++)
    {
      fmprb_mul(lamcostwopims[m],costwopims[m],lambda,PREC);
      fmprb_exp(explamcostwopims[m],lamcostwopims[m],PREC);
    }
  for(j=0;j<=k;j++)
    {
      fmprb_set_ui(sigmas[j],0);
      for(m=0;m<=k;m++)
	{
	  mj=(m*j)%(k+1);
	  fmprb_mul(temp,explamcostwopims[m],costwopims[mj],PREC);
	  fmprb_add(temp1,sigmas[j],temp,PREC);
	  fmprb_set(sigmas[j],temp1);
	}
      fmprb_div_ui(temp,sigmas[j],k+1,PREC);
      fmprb_set(sigmas[j],temp);
    }
  /*
  for(m=0;m<=k;m++)
    {
      printf("sigma[%lu]=",m);
      fmprb_printd(sigmas[m],10);
      printf("\n");
    }
  */
}


void make_rhos()
{
  uint64_t j,m,mj;
  for(j=0;j<=k;j++)
    {
      fmprb_set_ui(rhos[j],0);
      for(m=0;m<=k;m++)
	{
	  mj=(m*j)%(k+1);
	  fmprb_mul(temp1,sintwopims[m],sintwopims[mj],PREC);
	  fmprb_mul(temp,temp1,explamcostwopims[m],PREC);
	  fmprb_add(temp1,rhos[j],temp,PREC);
	  fmprb_set(rhos[j],temp1);
	}
      fmprb_div_ui(temp,rhos[j],k+1,PREC);
      fmprb_neg(temp1,temp);
      fmprb_mul(rhos[j],temp1,lambda,PREC);
    }
  /*
  for(m=0;m<=k;m++)
    {
      printf("rhos[%lu]=",m);
      fmprb_printd(rhos[m],10);
      printf("\n");
    }
  */
}

void make_bs()
{
  fmprb_init(bs[0]);
  fmprb_set(bs[0],sigmas[0]);
  uint64_t j;
  for(j=1;j<=k;j++)
    {
      fmprb_mul_ui(temp,sigmas[j],k+1-j,PREC);
      fmprb_sub(temp1,temp,rhos[j],PREC);
      fmprb_mul_2exp_si(temp,temp1,1);
      fmprb_init(bs[j]);
      fmprb_div_ui(bs[j],temp,k+1,PREC);
    }
  /*
  for(j=0;j<=k;j++)
    {
      printf("b[%lu]=",j);
      fmprb_printd(bs[j],10);
      printf("\n");
    }
  */
}

void make_U()
{
  uint64_t n,l,i;
  int64_t j,n0,j1,j2;
  for(n=0;n<k;n++)
    for(l=0;l<k;l++)
      fmprb_set_ui(fmprb_mat_entry(U1,n,l),0);

  for(j=1;j<=k;j++)
    {
      n0=(j-1)%2;
      for(n=n0;n<k;n+=2)
	{
	  j1=(j-1+n)>>1;
	  j2=(j-1-n);
	  if(j2<0)j2=-j2;
	  j2>>=1;
	  fmprb_mul_2exp_si(temp,bs[j-1],-1);
	  fmprb_add(temp1,fmprb_mat_entry(U1,n,j1),temp,PREC);
	  fmprb_set(fmprb_mat_entry(U1,n,j1),temp1);
	  fmprb_add(temp1,fmprb_mat_entry(U1,n,j2),temp,PREC);
	  fmprb_set(fmprb_mat_entry(U1,n,j2),temp1);
	}
    }
  /*

  for(n=0;n<k;n++)
    for(l=0;l<k;l++)
      {
	fmprb_set_ui(temp,0);
	j=2*l;
	j-=n;
	if((j>=0)&&(j<=k))
	  fmprb_set(temp,bs[j]);
	j=2*l+n;
	if((j>=0)&&(j<=k))
	  fmprb_add(temp1,temp,bs[j],PREC);
	else
	  fmprb_set(temp1,temp);
	j=n;
	j-=2*l;
	if((j>=0)&&(j<=k))
	  fmprb_add(temp,temp1,bs[j],PREC);
	else
	  fmprb_set(temp,temp1);
	fmprb_mul_2exp_si(fmprb_mat_entry(U1,n,l),temp,-1);
      }
  */
  //printf("\nU before exponentiation:\n");fmprb_mat_printd(U1,10);
  fmprb_mat_pow_ui(U,U1,L,PREC);
  //printf("\nU after exponentiation:\n");fmprb_mat_printd(U,10);

  fmprb_set(temp,fmprb_mat_entry(U,0,0));
  for(l=1;l<k;l++)
    {
      fmprb_add(temp1,temp,fmprb_mat_entry(U,0,l),PREC);
      fmprb_set(temp,temp1);
    }
    
  // here temp=exp(L psi(lambda))

  //printf("sum U_{0,n}=");fmprb_printd(temp1,30);printf("\n");


  fmprb_log(temp,temp1,PREC); // temp1=L psi(lambda)
  fmprb_div_ui(temp1,temp,L,PREC);
  //printf("/2^L=");fmprb_printd(temp1,10);printf("\n");

  fmprb_add(temp,temp1,c,PREC);
  fmprb_div(temp1,temp,lambda,PREC);
    
  fmprb_printd(temp1,30);printf("\n");
  

}

int main(int argc, char**argv)
{

  if(argc!=8)
    {
      printf("Usage:-%s <lambda num start> <lambda num end> <lambda den> <L> <k> <cnum> <cden>.\n",argv[0]);
      printf("compute (psi(lambda)+<cnum>/<cden> log 2)/lambda for lambda = <num start>/<den>..<num end>/<den>.\n");
      exit(0);
    }

  L=atol(argv[4]);
  k=atol(argv[5]);
  cnum=atol(argv[6]);
  cden=atol(argv[7]);

  fmprb_init(lambda);
  uint64_t num=atol(argv[1]),fin=atol(argv[2]),den=atol(argv[3]);
  lambdad=num/(double) den;
  printf("lambda=%lu/%lu: ",num,den);
  fmprb_set_ui(temp,num);
  fmprb_div_ui(lambda,temp,den,PREC);
  make_sigmas();
  make_rhos();
  make_bs();
  make_U();
  num++;
  for(;num<=fin;num++)
    {
      printf("lambda=%lu/%lu: ",num,den);
      lambdad=num/(double) den;
      fmprb_set_ui(temp,num);
      fmprb_div_ui(lambda,temp,den,PREC);
      make_new_sigmas();
      make_rhos();
      make_bs();
      make_U();
    }
  return(0);
}
