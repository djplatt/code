#include "structs.h"
#include "acb_fft.h"

// read the parameters for this computation
bool read_comp(L_comp_t *L, FILE *cfile, uint64_t prec)
{
  arb_t tmp;
  arb_init(tmp);

  L->N=1<<10; // size of FFT
  L->M=100; // Number of terms in sum over G
  L->K=40; // Highest differential of G to use
  L->eta=0.0;
  arb_init(L->delta);
  arb_const_pi(L->delta,prec);
  arb_mul_2exp_si(L->delta,L->delta,-1);
  arb_set_d(tmp,1.0-L->eta);
  arb_mul(L->delta,L->delta,tmp,prec);

  L->w=(acb_t *)malloc(sizeof(acb_t)*L->N/2);
  if(!L->w)
    {
      printf("Error allocating memory for w. Exiting.\n");
      return(0);
    }

  for(uint64_t i=0;i<L->N/2;i++)
    acb_init(L->w[i]);
  acb_initfft(L->w,L->N,prec);


  arb_clear(tmp);
  return(true);
}

// read conductor from infile
uint64_t read_cond(FILE *infile)
{
  return(11);
}

uint64_t read_u64(FILE *infile)
{
  uint64_t r;
  if(fscanf(infile,"%lu\n",&r)!=1)
    return(0);
  else
    return(r);
}

// read r from file
uint64_t read_r(FILE *infile)
{
  return(read_u64(infile));
}

  

bool read_arb(arb_ptr res, FILE *infile)
{
  static bool init=false;
  static fmpz_t a,b;
  static mpz_t x,e;
  static arb_t radius;
  if(!init)
    {
      init=true;
      fmpz_init(a);
      fmpz_init(b);
      mpz_init(x);
      mpz_init(e);
    }

 
  if(!mpz_inp_str(x,infile,10))
    return(false);
  if(!mpz_inp_str(e,infile,10))
    return(false);
  fmpz_set_mpz(a,x);
  fmpz_set_mpz(b,e);
  arb_set_fmpz_2exp(res,a,b);

  if(!mpz_inp_str(x,infile,10))
    return(false);
  if(!mpz_inp_str(e,infile,10))
    return(false);
  fmpz_set_mpz(a,x);
  fmpz_set_mpz(b,e);
  arb_set_fmpz_2exp(radius,a,b);

  arb_add_error(res,radius);

  return(true);
}

bool read_acb(acb_ptr res, FILE *infile)
{
  if(!read_arb(acb_realref(res),infile))
    return(false);
  return(read_arb(acb_imagref(res),infile));
}

bool read_mu(acb_ptr res, FILE *infile)
{
  return(read_acb(res,infile));
}

int64_t read_m(FILE *infile)
{
  return(read_u64(infile));
}

bool read_lambda(arb_ptr res, FILE *infile)
{
  return(read_arb(res,infile));
}

bool read_residue(acb_ptr res, FILE *infile)
{
  return(read_acb(res,infile));
}


bool read_C_alpha(double *C, double *alpha, FILE *infile)
{
  return(fscanf(infile,"%f %f\n",C,alpha)==2);
}

uint64_t read_max_K(FILE *infile)
{
  return(40);
}
double read_max_B(FILE *infile)
{
  return(128.0);
}
int64_t read_low_i(FILE *infile)
{
  return(-24);
}
uint64_t read_num_i(FILE *infile)
{
  return(1024);
}

double read_gap(FILE *infile)
{
  double res;
  if(fscanf(infile,"%e",&res)!=1)
    return(0.0);
  else
    return(res);
}

// hardwired for L function of quadratic character mod 5
bool read_family(L_family_t *L_family, L_comp_t *Lc, FILE *ffile, uint64_t prec)
{
  //printf("Defaulting to family of odd dirichlet characters mod 11.\n");
  //fflush(stdout);
  arb_t tmp;
  arb_init(tmp);

  // should read N and r from file
  L_family->r=read_r(ffile);
  if(!L_family->r)
    return(false);
  printf("r=%lu\n",L_family->r);
  L_family->mus=(acb_t *) malloc(sizeof(acb_t)*L_family->r);
  if(!L_family->mus)
    {
      printf("Error allocating memory for mus in read_family.\n");
      return(false);
    }
  L_family->nus=(arb_t *) malloc(sizeof(arb_t)*L_family->r);
  if(!L_family->nus)
    {
      printf("Error allocating memory for nus in read_family.\n");
      return(false);
    }
  for(uint64_t i=0;i<L_family->r;i++)
    {
      acb_init(L_family->mus[i]);
      arb_init(L_family->nus[i]);
    }

  for(uint64_t i=0;i<L_family->r;i++)
    if(!read_mu(L_family->mus[i],ffile))
      {
	printf("Error reading mu from family file. Exiting.\n");
	//printf("mu[%lu] was ",i);acb_printd(L_family->mus[i],10);printf("\n");
	return(false);
      }
    else
      {
	printf("mu[%lu]=",i);
	acb_printd(L_family->mus[i],10);
	printf("\n");
      }
  // Lemma 5.2 defines mu
  acb_init(L_family->mu);
  acb_zero(L_family->mu);
  for(uint64_t i=0;i<L_family->r;i++)
    acb_add(L_family->mu,L_family->mu,L_family->mus[i],prec);
  acb_div_ui(L_family->mu,L_family->mu,L_family->r,prec);
  arb_set_d(tmp,0.5);
  arb_sub(acb_realref(L_family->mu),acb_realref(L_family->mu),tmp,prec);

  // Lemma 5.2 defines nu_j
  arb_set_ui(tmp,L_family->r*2);
  arb_inv(tmp,tmp,prec);
  for(uint64_t i=0;i<L_family->r;i++)
    {
      arb_sub_ui(L_family->nus[i],acb_realref(L_family->mus[i]),1,prec);
      arb_mul_2exp_si(L_family->nus[i],L_family->nus[i],-1);
      arb_add(L_family->nus[i],L_family->nus[i],tmp,prec);
    }

  // should read m from file
  L_family->m=read_m(ffile); // No. of poles
  if(L_family->m<0)
    return(false);
  printf("m=%lu\n",L_family->m);
  if(L_family->m>0)
    {
      printf("Need to implement sum over residues in 5-2. Exiting.\n");
      exit(0);
    }
  if(L_family->m>0) // there are some poles
    {
      L_family->lambdas=(arb_t *) malloc(sizeof(arb_t)*L_family->r);
      if(!L_family->lambdas)
	{
	  printf("Error allocating memory for lambdas in read_family.\n");
	  return(false);
	}
      L_family->residues=(acb_t *) malloc(sizeof(acb_t)*L_family->r);
      if(!L_family->residues)
	{
	  printf("Error allocating memory for residues in read_family.\n");
	  return(false);
	}
      for(uint64_t i=0;i<L_family->m;i++)
	{
	  arb_init(L_family->lambdas[i]);
	  if(!read_lambda(L_family->lambdas[i],ffile))
	    return(false);
	  acb_init(L_family->residues[i]);
	  if(!read_residue(L_family->residues[i],ffile))
	    return(false);
	}
    }

  if(!read_C_alpha(&L_family->C,&L_family->alpha,ffile))
    return(false);

  // c,c_dash defined in Lemma 5.4
  arb_init(L_family->c);
  arb_init(L_family->c_dash);
  arb_set_d(tmp,0.5+L_family->alpha);
  arb_add(L_family->c,acb_realref(L_family->mu),tmp,prec);
  arb_mul_ui(tmp,L_family->c,L_family->r,prec);
  arb_mul_2exp_si(tmp,tmp,-1);
  arb_sub_ui(tmp,tmp,1,prec);
  if(arb_contains_negative(tmp))
    {
      if(arb_contains_positive(tmp))
	{
	  printf("cr/2-1 straddles zero. Use more precision. Exiting.\n");
	  return(false);
	}
      arb_zero(L_family->c_dash);
    }
  else
    arb_set(L_family->c_dash,tmp);

  L_family->file_two_pi_by_B=read_gap(ffile);
  if(L_family->file_two_pi_by_B==0.0)
    return(false);
  printf("File has 2pi/B=%e\n",L_family->two_pi_by_B);

  L_family->max_K=read_max_K(ffile);
  if(!L_family->max_K)
    return(false);
  L_family->max_B=read_max_B(ffile);
  if(!L_family->max_B)
    return(false);
  L_family->low_i=read_low_i(ffile);
  if(L_family->low_i>0)
    return(false);
  L_family->num_i=read_num_i(ffile);
  if(!L_family->num_i)
    return(false);
  arb_clear(tmp);
  return(true);
}

bool read_ans(acb_t *ans, uint64_t M, FILE *infile)
{
  printf("Defaulting to the quadratic character mod 11.\n");
  for(uint64_t i=0;i<M;i+=11)
    acb_set_ui(ans[i],1);
  for(uint64_t i=1;i<M;i+=11)
    acb_set_si(ans[i],-1);
  for(uint64_t i=3;i<M;i+=11)
    acb_set_ui(ans[i],1);
  for(uint64_t i=7;i<M;i+=11)
    acb_set_si(ans[i],-1);
  for(uint64_t i=4;i<M;i+=11)
    acb_set_ui(ans[i],1);
  for(uint64_t i=9;i<M;i+=11)
    acb_set_si(ans[i],-1);
  for(uint64_t i=8;i<M;i+=11)
    acb_set_ui(ans[i],1);
  for(uint64_t i=6;i<M;i+=11)
    acb_set_si(ans[i],-1);
  for(uint64_t i=2;i<M;i+=11)
    acb_set_ui(ans[i],1);
  for(uint64_t i=5;i<M;i+=11)
    acb_set_si(ans[i],-1);
  for(uint64_t i=10;i<M;i+=11)
    acb_zero(ans[i]);
  return(true);
}

bool init_func(L_func_t *L, uint64_t M)
{
  for(uint64_t side=0;side<2;side++)
    {
      acb_init(L->epsilon[side]);
      L->ans[side]=(acb_t *)malloc(sizeof(acb_t)*M);
      if(!L->ans[side])
	{
	  printf("Error allocating memory for a_n's. Exiting.\n");
	  return(false);
	}
      for(uint64_t i=0;i<M;i++)
	acb_init(L->ans[side][i]);
    }
  return(true);
}


// M values of a_n
bool read_func(L_func_t *L, FILE *fufile, uint64_t prec, uint64_t M)
{
  // check that the function file contains at least M coefficients!

  printf("In read func.\n");fflush(stdout);
  L->N=read_cond(fufile);
  if(!L->N)
    {
      printf("Error reading conductor for function.\n");
      return(false);
    }
  L->dc=sqrt((double) L->N);

  L->self_dual_p=true;

  // read the M values a_n for n=1..M into ans[0]..ans[M-1]
  // hardwired for quadratic character mod 11
  if(!read_ans(L->ans[0],M,fufile))
    return(false);



  arb_t tmp;
  arb_init(tmp);

  if(L->self_dual_p)
    for(uint64_t i=1;i<M;i++)
      {
	arb_sqrt_ui(tmp,i+1,prec);
	acb_div_arb(L->ans[0][i],L->ans[0][i],tmp,prec);
      }
  else
    {
      if(!read_ans(L->ans[1],M,fufile))
	return(false);
      for(uint64_t i=1;i<M;i++)
	{
	  arb_sqrt_ui(tmp,i+1,prec);
	  acb_div_arb(L->ans[0][i],L->ans[0][i],tmp,prec);
	  acb_div_arb(L->ans[1][i],L->ans[1][i],tmp,prec);
	}
    }

  arb_clear(tmp);

  return(true);
}

int64_t calc_m(uint64_t i, double B, double dc)
{
  double x=log((double) i/dc)*B/(2*M_PI);
  return(round(x));
}

bool read_Gs(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, FILE *infile, uint64_t prec)
{
  char line[1024];
  for(uint64_t k=0;k<Lc->K;k++)
    for(uint64_t n=0;n<Lc->N;n++)
      {
	if(fscanf(infile,"%s\n",line)!=1)
	  {
	    printf("Error reading n=%lu k=%lu.\n",k,n);
	    return(false);
	  }
	if(arb_set_str(acb_realref(Lc->Gs[k][n]),line,prec))
	  {
	    printf("Error converting %s to arb_t.",line);
	    return(false);
	  }
      }
  return(true);
}

bool finalise_comp(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, FILE *ffile, uint64_t prec)
{
  if(Lc->K>Lf->max_K)
    {
      printf("Insufficient differentials (%lu) to support K=%lu\n",
	     Lf->max_K,Lc->K);
      return(false);
    }
  Lc->B=Lf->max_B;
  printf("working B set to %10.8e\n",Lc->B);
  Lc->A=Lc->N/Lc->B;
  printf("working A set to %10.8e\n",Lc->A);

  Lc->offset=calc_m(1,Lc->B,Lfu->dc);
  if(Lc->offset!=Lf->low_i)
    {
      printf("G values need to go down to %ld. We only have down to %ld.\n",
	     Lc->offset,Lf->low_i);
      return(false);
    }
  if(Lc->N!=Lf->num_i)
    {
      printf("G values need to go up to %ld. We only have up to %ld.\n",
	     Lc->N,Lf->num_i);
      return(false);
    }


  // allocate vector of vectors to hold G^(k)(i/B)
  Lc->Gs=(acb_t **)malloc(sizeof(acb_t *)*(Lc->K+1));
  for(uint64_t k=0;k<=Lc->K;k++)
    {
      Lc->Gs[k]=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
      for(uint64_t n=0;n<Lc->N;n++)
	acb_init(Lc->Gs[k][n]);
    }

  // now we read the G^(k)(2 pi i/B) into the right slots
  // G(offset/B) into [0], ... G(N/2B)into [Lc->G_len-1]
  if(!read_Gs(Lc,Lf,Lfu,ffile,prec))
    return(false);

  // pre-FFT the Gs (1/3 of a convolve)
  for(uint64_t k=0;k<=Lc->K;k++)
    acb_fft(Lc->Gs[k],Lc->N,Lc->w,prec);

  Lc->ms=(int64_t *)malloc(sizeof(int64_t)*Lc->M);
  for(uint64_t m=1;m<=Lc->M;m++)
    Lc->ms[m-1]=calc_m(m,Lc->B,Lfu->dc);
  
  for(uint64_t m=0;m<10;m++) printf("ms[%lu]=%ld\n",m,Lc->ms[m]);

  arb_init(Lc->arb_B);
  arb_set_d(Lc->arb_B,Lc->B);
  arb_inv(Lc->arb_B,Lc->arb_B,prec); // 1/B
  /*
  arb_init(Lc->two_pi_by_B);
  arb_const_pi(Lc->two_pi_by_B,prec);
  arb_mul_2exp_si(Lc->two_pi_by_B,Lc->two_pi_by_B,1);
  arb_mul(Lc->two_pi_by_B,Lc->two_pi_by_B,Lc->arb_B,prec);
  */

  Lc->ums=(arb_t *)malloc(sizeof(arb_t)*Lc->M);
  for(uint64_t m=0;m<Lc->M;m++)
    {
      arb_init(Lc->ums[m]);
      arb_mul_si(Lc->ums[m],Lc->two_pi_by_B,Lc->ms[m],prec);
    }
  
  for(uint64_t m=0;m<10;m++)
    {
      printf("ums[%lu]=",m);
      arb_printd(Lc->ums[m],10);
      printf("\n");
    }
  

  Lc->sks=(arb_t **)malloc(sizeof(arb_t *)*Lc->K);
  for(uint64_t k=1;k<=Lc->K;k++)
    {
      Lc->sks[k-1]=(arb_t *)malloc(sizeof(arb_t)*Lc->M);
      for(uint64_t m=0;m<Lc->M;m++)
	arb_init(Lc->sks[k-1][m]);
    }
  arb_t tmp,tmp1,arb_sc;
  arb_init(tmp);arb_init(tmp1);arb_init(arb_sc);
  arb_sqrt_ui(tmp,Lfu->N,prec);
  arb_inv(arb_sc,tmp,prec);
  for(uint64_t m=0;m<Lc->M;m++)
    {
      arb_mul_ui(tmp,arb_sc,m+1,prec);
      arb_log(tmp1,tmp,prec);
      arb_sub(Lc->sks[0][m],tmp1,Lc->ums[m],prec);
    }

  arb_clear(tmp);arb_clear(tmp1);arb_clear(arb_sc);

  for(uint64_t k=2;k<=Lc->K;k++)
    for(uint64_t m=0;m<Lc->M;m++)
      arb_pow_ui(Lc->sks[k-1][m],Lc->sks[0][m],k,prec);

  Lc->kres=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
  Lc->skm=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
  Lc->res[0]=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
  Lc->res[1]=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
  for(uint64_t n=0;n<Lc->N;n++)
    {
      acb_init(Lc->kres[n]);
      acb_init(Lc->skm[n]);
      acb_init(Lc->res[0][n]);
      acb_init(Lc->res[1][n]);
    }

  return(true);
}  

bool set_errors(L_error_t *errs)
{
  arb_init(errs->lem54);
  arb_init(errs->lem56);
  arb_init(errs->lem57);
  arb_init(errs->eq59);
  printf("Setting all error terms to zero!.\n");
  arb_zero(errs->lem54);
  arb_zero(errs->lem56);
  arb_zero(errs->lem57);
  arb_zero(errs->eq59);
  return(true);
}

// where does the entry for m+1 go?
uint64_t bucket(uint64_t m, L_comp_t *Lc)
{
  int64_t res=(Lc->N-(Lc->ms[m]-Lc->offset))%Lc->N;
  return(res);
}


bool do_convolves(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  acb_t tmp;
  acb_init(tmp);
  for(uint64_t side=0;side<(Lfu->self_dual_p ? 1 : 2);side++)
    {
      for(uint64_t n=0;n<Lc->N;n++)
	acb_zero(Lc->res[side][n]);
      // let's do k=0
      for(uint64_t n=0;n<Lc->N;n++)
	acb_zero(Lc->skm[n]);
      for(uint64_t m=0;m<Lc->M;m++)
	{
	  uint64_t b=bucket(m,Lc);			
	  acb_add(Lc->skm[b],Lc->skm[b],Lfu->ans[side][m],prec);
	  printf("Put ");acb_printd(Lfu->ans[side][m],20);printf(" into bucket %lu\n",b);
	}
      acb_convolve1(Lc->res[side],Lc->skm,Lc->Gs[0],Lc->N,Lc->w,prec);
      
      for(uint64_t k=1;k<=Lc->K;k++)
	{
	  for(uint64_t n=0;n<Lc->N;n++)
	    acb_zero(Lc->skm[n]);
	  for(uint64_t m=0;m<Lc->M;m++)
	    {
	      uint64_t b=bucket(m,Lc);			
	      acb_mul_arb(tmp,Lfu->ans[side][m],Lc->sks[k-1][m],prec);
	      acb_add(Lc->skm[b],Lc->skm[b],tmp,prec);
	      printf("Put ");acb_printd(tmp,20);printf(" into bucket %lu\n",b);
	    }
	  acb_convolve1(Lc->kres,Lc->skm,Lc->Gs[k],Lc->N,Lc->w,prec);
	  for(uint64_t n=0;n<=Lc->N/2;n++)
	    acb_add(Lc->res[side][n],Lc->res[side][n],Lc->kres[n],prec);
	}
    }		    
  acb_clear(tmp);
  return(true);
}

void arb_intersect(arb_ptr res, arb_ptr a, arb_ptr b, uint64_t prec)
{
  arf_t al,ah,bl,bh;
  arf_init(al);arf_init(ah);arf_init(bl);arf_init(bh);
  arb_get_interval_arf(al,ah,a,prec);
  arb_get_interval_arf(bl,bh,b,prec);

  if(arf_cmp(al,bl)<0)
    arf_set(al,bl);
  if(arf_cmp(ah,bh)>0)
    arf_set(ah,bh);
  if(arf_cmp(al,ah)>=0) // enpty intersection
    {
      printf("arb_intersect produced empty intersection. Exiting.\n");
      printf("intersected ");arb_printd(a,10);
      printf("\n with ");arb_printd(b,10);
      printf("\n");
      exit(0);
    }
  arb_set_interval_arf(res,al,ah,prec);
  arf_clear(al);arf_clear(ah);arf_clear(bl);arf_clear(bh);
}

bool final_ifft(L_comp_t *Lc, L_func_t *Lfu, uint64_t prec)
{
  // first use res[side][0] and res[side][N/2] to figure our epsilon
  arb_t th1,th2,th;
  arb_init(th1);arb_init(th2);arb_init(th);

  for(uint64_t side=0;side<(Lfu->self_dual_p ? 1 : 2);side++)
    {
      acb_arg(th1,Lc->res[side][0],prec);
      acb_arg(th2,Lc->res[side][Lc->N/2],prec);
      arb_set(th,th1); //arb_intersect(th,th1,th2,prec);
      arb_neg(th,th);
      arb_sin_cos(th1,th2,th,prec);
      arb_set(acb_realref(Lfu->epsilon[side]),th2);
      arb_set(acb_imagref(Lfu->epsilon[side]),th1);
      for(uint64_t n=0;n<Lc->N/2;n++)
	acb_mul(Lc->res[side][n],Lc->res[side][n],Lfu->epsilon[side],prec);
      for(uint64_t n=Lc->N/2+1;n<Lc->N;n++)
	acb_conj(Lc->res[side][n],Lc->res[side][Lc->N-n]);
      acb_ifft(Lc->res[side],Lc->N,Lc->w,prec);
      for(uint64_t n=0;n<Lc->N;n++)
	acb_mul_arb(Lc->res[side][n],Lc->res[side][n],Lc->two_pi_by_B,prec);
    }
      
  arb_clear(th1),arb_clear(th2);arb_clear(th);
  return(true);
}

int main(int argc, char **argv)
{
  printf("Command line:- %s",argv[0]);
  for(uint64_t i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=6)
    {
      printf("Usage:- %s <prec> <comp file> <family file> <func file>  <out file>.\n",argv[0]);
      return(0);
    }

  uint64_t arb_prec=atol(argv[1]);
  L_family_t L_family;
  FILE *ffile=fopen(argv[3],"r");
  if(!ffile)
    {
      printf("Failed to open family file %s. Exiting.\n",argv[3]);
      return(0);
    }
  L_comp_t L_comp;
  FILE *cfile;
  L_func_t L_func;
  FILE *fufile;
  L_error_t L_errors;
  
  if(!read_comp(&L_comp,cfile,arb_prec))
    {
      printf("Error reading computation file. Exiting.\n");
      return(0);
    }
  
  if(!read_family(&L_family,&L_comp,ffile,arb_prec))
    {
      printf("Error reading L function family file. Exiting.\n");
      return(0);
    }
  
  print_family(&L_family);


  if(!set_errors(&L_errors))
    {
      printf("Failed to set error terms. Exiting.\n");
      return(0);
    }


  if(!init_func(&L_func,L_comp.M))
    {
      printf("Error initialising L-func structure. Exiting.\n");
      return(0);
    }

  // plan would be to iterate here over functions in same family
  if(!read_func(&L_func,fufile,arb_prec,L_comp.M))
    {
      printf("Error reading function file. Exiting.\n");
      return(0);
    }

  if(!finalise_comp(&L_comp, &L_family, &L_func, ffile, arb_prec))
    {
      printf("Error finalising computation parameters. Exiting.\n");
      return(0);
    }

  printf("Doing convolutions...\n");
  if(!do_convolves(&L_comp,&L_family,&L_func,arb_prec))
    {
      printf("Error doing computation. Exiting.\n");
      return(0);
    }
  printf("Finished convolutions...\n");

  for(uint64_t side=0;side<(L_func.self_dual_p ? 1 : 2);side++)
    for(uint64_t n=0;n<L_comp.N;n++)
      {
	printf("L[%lu][%lu]=",side,n);
	arb_printd(acb_realref(L_comp.res[side][n]),20);
	printf("\n");
      }   


  if(!final_ifft(&L_comp,&L_func,arb_prec))
    {
      printf("Error doing final IFFT. Exiting.\n");
      return(0);
    }

  for(uint64_t side=0;side<(L_func.self_dual_p ? 1 : 2);side++)
    for(uint64_t n=0;n<L_comp.N;n++)
      {
	printf("L[%lu][%lu]=",side,n);
	arb_printd(acb_realref(L_comp.res[side][n]),20);
	printf("\n");
      }   

  return(0);
}
