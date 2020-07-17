
// gcc pintz.c -I /home/dave/flint2-trunk/ -lgmp -lflint -larb -lpthread

#include "inttypes.h"
#include "fmpr.h"
#include "fmprb.h"
#include "fmprb_mat.h"
//#include "fmpcb.h"

uint64_t L,k;

double lambdad,cd;

#define PREC ((uint64_t) 200)

fmprb_t *sigmas,*rhos,*bs,*twopims,*costwopims,*sintwopims,
  *lamcostwopims,*explamcostwopims,pi,temp,temp1;
fmpr_t lambda,c;

fmprb_mat_t U,U1;

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
  fmprb_init(pi);
  fmprb_const_pi(pi,PREC);

  uint64_t m;
  for(m=0;m<=k;m++)
    {
      fmprb_mul_ui(temp,pi,m<<1,PREC);
      fmprb_init(twopims[m]);
      fmprb_div_ui(twopims[m],temp,k+1,PREC);
      fmprb_init(costwopims[m]);
      fmprb_init(sintwopims[m]);
      fmprb_sin_cos(sintwopims[m],costwopims[m],twopims[m],PREC);
      fmprb_init(explamcostwopims[m]);
      fmprb_init(lamcostwopims[m]);
    }
}

void make_sigmas()
{
  calc_twopims();
  uint64_t j,m,mj;
  for(m=0;m<=k;m++)
    {
      fmprb_mul_fmpr(lamcostwopims[m],costwopims[m],lambda,PREC);
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
      fmprb_mul_fmpr(lamcostwopims[m],costwopims[m],lambda,PREC);
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
      fmprb_div_si(temp,rhos[j],-k-1,PREC);
      fmprb_mul_fmpr(rhos[j],temp,lambda,PREC);
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
  for(j=0;j<=k;j++)
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
  int64_t j;
  for(n=0;n<k;n++)
    for(l=0;l<k;l++)
      {
	fmprb_set_ui(temp,0);
	j=2*l-n;
	if((j>0)&&(j<k))
	  fmprb_set(temp,bs[j]);
	j=2*l+n;
	if((j>0)&&(j<k))
	  fmprb_add(temp1,temp,bs[j],PREC);
	else
	  fmprb_set(temp1,temp);
	j=n-2*l;
	if((j>0)&&(j<k))
	  fmprb_add(temp,temp1,bs[j],PREC);
	else
	  fmprb_set(temp,temp1);
	fmprb_mul_2exp_si(fmprb_mat_entry(U1,n,l),temp,-1);
      }

  //fmprb_mat_printd(U1,10);

  fmprb_mat_pow_ui(U,U1,L,PREC);

  fmprb_set(temp,fmprb_mat_entry(U,0,0));
  for(l=1;l<k;l++)
    {
      fmprb_add(temp1,temp,fmprb_mat_entry(U,0,l),PREC);
      fmprb_set(temp,temp1);
    }
    
  //printf("exp(%lu psi(%f))<=",L,lambdad);
  //fmprb_printd(temp,10);
  //printf("\n");
  
  fmprb_log(temp1,temp,PREC);
  fmprb_div_ui(temp,temp1,L,PREC);
  //printf("psi(%f)<=",lambdad);
  //fmprb_printd(temp,10);
  //printf("\n");
  //return;

  fmprb_add_fmpr(temp1,temp,c,PREC);
  fmprb_div_fmpr(temp,temp1,lambda,PREC);

  printf("(psi(%f)+%f)/%f<=",lambdad,cd,lambdad);
  fmprb_printd(temp,10);
  printf("\n");
  

}

int main(int argc, char**argv)
{

  if(argc!=7)
    {
      printf("Usage:-%s <lambda num start> <lambda num end> <lambda den> <L> <k> <c>.\n",argv[0]);
      exit(0);
    }

  L=atol(argv[4]);
  k=atol(argv[5]);
  cd=atof(argv[6]);
  fmpr_init(lambda);
  fmpr_init(c);
  fmpr_set_d(c,cd);
  uint64_t num=atol(argv[1]),fin=atol(argv[2]),den=atol(argv[3]);
  lambdad=num/(double) den;
  fmpr_set_d(lambda,lambdad);
  make_sigmas();
  make_rhos();
  make_bs();
  make_U();
  num++;
  for(;num<=fin;num++)
    {
      lambdad=num/(double) den;
      fmpr_set_d(lambda,lambdad);
      make_new_sigmas();
      make_rhos();
      make_bs();
      make_U();
    }
  return(0);
}
