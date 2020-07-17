// 3.0 only copy those G values we actually need
// means we can't pre-FFT the G data
// only store (log(m/sqrt(N))-um) not all powers

// 4.0 compute Taylor series the obvious way (no convolution)
// see do_convolves_no_fft
// also remove the two sides (one is just the conjugate of the other anyway)

// 6.0 add error term for truncating sum over coefficients
// and F_hat twiddle error

// added cheap upsample at iFFT stage

#include "structs7.0.h"
#include "acb_fft.h"
#include "error7.0.h"

#define OUTPUT_RATIO (8) // we will analyse this portion of B
#define TURING_RATIO (64) // we will use this portion for Turing

#define MAX_M (1LL<<18) // maximum no. of Dirichlet coeffs we can handle

// read the parameters for this computation
// all we do is set FFT length N and eta=0
bool read_comp(L_comp_t *L, FILE *cfile, uint64_t prec)
{
  arb_t tmp;
  arb_init(tmp);

  L->N=1<<10; // size of FFT
  L->NN=1<<16;
  printf("Set FFT length to %lu. Not checking whether this is sufficient. Maybe need FFT Length > 2*(# G values used).\n",L->N);

  L->G=(acb_t *)malloc(sizeof(acb_t)*L->N);
  for(uint64_t n=0;n<L->N;n++)
    acb_init(L->G[n]);

  L->eta=0.0;
  arb_init(L->delta);
  arb_const_pi(L->delta,prec);
  arb_mul_2exp_si(L->delta,L->delta,-1); // pi/2
  arb_set_d(tmp,1.0-L->eta);
  arb_mul(L->delta,L->delta,tmp,prec); // (1-eta)pi/2

  L->w=(acb_t *)malloc(sizeof(acb_t)*L->N/2); // twiddles for FFT
  if(!L->w)
    {
      printf("Error allocating memory for w. Exiting.\n");
      return(0);
    }

  for(uint64_t i=0;i<L->N/2;i++)
    acb_init(L->w[i]);
  acb_initfft(L->w,L->N,prec); // set twiddles

  L->ww=(acb_t *)malloc(sizeof(acb_t)*L->NN/2); // twiddles for FFT
  if(!L->ww)
    {
      printf("Error allocating memory for ww. Exiting.\n");
      return(0);
    }

  for(uint64_t i=0;i<L->NN/2;i++)
    acb_init(L->ww[i]);
  acb_initfft(L->ww,L->NN,prec); // set twiddles

  L->ms=(int64_t *)malloc(sizeof(int64_t)*MAX_M);

  L->ums=(arb_t *)malloc(sizeof(arb_t)*MAX_M);
  for(uint64_t m=0;m<MAX_M;m++)
    arb_init(L->ums[m]);

  arb_clear(tmp);
  return(true);
}


uint64_t read_u64(FILE *infile)
{
  uint64_t r;
  if(fscanf(infile,"%lu\n",&r)!=1)
    return(0);
  else
    return(r);
}

#define BAD_64 (1LL<<62)
int64_t read_64(FILE *infile)
{
  int64_t r;
  if(fscanf(infile,"%ld\n",&r)!=1)
    return(BAD_64);
  else
    return(r);
}


// read r from file
uint64_t read_r(FILE *infile)
{
  return(read_u64(infile));
}

bool read_err(arb_ptr res, FILE *infile)
{  
  static bool init=false;
  static fmpz_t a,b;
  static mpz_t x,e;
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
  return(true);
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
      arb_init(radius);
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
  double mu;
  if(!fscanf(infile,"%f",&mu)==1)
    {
      printf("Error reading mu from file.\n");
      return(false);
    }
  arb_set_d(acb_realref(res),mu);
  arb_zero(acb_imagref(res));
  return(true);
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


bool read_C_alpha_line(FILE *infile)
{
  double x[3];
  return(fscanf(infile,"%lf %lf %lf\n",x,x+1,x+2)==3);
}

int64_t read_low_i(FILE *infile)
{
  return(read_64(infile));
}

int64_t read_hi_i(FILE *infile)
{
  return(read_64(infile));
}

uint64_t read_max_K(FILE *infile)
{
  return(read_u64(infile));
}

double read_gap(FILE *infile)
{
  double res;
  if(fscanf(infile,"%lf",&res)!=1)
    return(0.0);
  else
    return(res);
}

bool read_eq59(arb_ptr res, FILE *infile)
{
  return(read_err(res,infile));
}

bool read_Gs(L_comp_t *Lc, L_family_t *Lf, FILE *infile, uint64_t prec)
{

  int64_t fi,i,j;uint64_t fk;
  arb_t G;
  arb_init(G);

  for(i=Lf->low_i,j=0;i<=Lf->hi_i;i++,j++)
    for(uint64_t k=0;k<Lf->max_K;k++)
      {
	if(fscanf(infile,"%ld %lu",&fi,&fk)!=2)
	  return(false);
	if(fi!=i)
	  return(false);
	if(fk!=k)
	  return(false);
	if(!read_arb(G,infile))
	  return(false);
	arb_set(Lf->Gs[k][j],G);
      }
  arb_clear(G);
  return(true);
}


// hardwired for L function of quadratic character mod 5
bool read_family(L_family_t *L_family, L_comp_t *Lc, L_error_t *Le, FILE *ffile, uint64_t prec)
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
  L_family->mus=(double *) malloc(sizeof(double)*L_family->r);
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
      arb_init(L_family->nus[i]);
    }

  for(uint64_t i=0;i<L_family->r;i++)
    if(fscanf(ffile,"%lf",L_family->mus+i)!=1)
      {
	printf("Error reading mu from family file. Exiting.\n");
	//printf("mu[%lu] was ",i);acb_printd(L_family->mus[i],10);printf("\n");
	return(false);
      }
    else
      printf("mu[%lu]=%10.8e\n",i,L_family->mus[i]);

  // Lemma 5.2 defines mu
  acb_init(L_family->mu);
  acb_zero(L_family->mu);
  for(uint64_t i=0;i<L_family->r;i++)
    {
      arb_set_d(tmp,L_family->mus[i]);
      acb_add_arb(L_family->mu,L_family->mu,tmp,prec);
    }
  acb_div_ui(L_family->mu,L_family->mu,L_family->r,prec);
  arb_set_d(tmp,0.5);
  arb_sub(acb_realref(L_family->mu),acb_realref(L_family->mu),tmp,prec);

  // lemma 3 defines nu_j and nu
  for(uint64_t i=0;i<2;i++)
    {
      arb_set_d(tmp,L_family->mus[i]-0.5);
      arb_mul_2exp_si(L_family->nus[i],tmp,-1);
    }
  for(uint64_t i=2;i<L_family->r;i++)
    {
      arb_set_d(tmp,L_family->mus[i]-1.0);
      arb_mul_2exp_si(L_family->nus[i],tmp,-1);
    }
  arb_init(L_family->nu);
  arb_set(L_family->nu,L_family->nus[0]);
  for(uint64_t i=1;i<L_family->r;i++)
    arb_add(L_family->nu,L_family->nu,L_family->nus[i],prec);
  arb_div_ui(L_family->nu,L_family->nu,L_family->r,prec);
  arb_mul_2exp_si(L_family->nu,L_family->nu,1);
  arb_set_d(tmp,0.5);
  arb_add(L_family->nu,L_family->nu,tmp,prec);

  if(!read_C_alpha_line(ffile))
    return(false);
  printf("Assuming Ramanujan bound.\n");
  //printf("C=%10.8e alpha=%10.8e\n",L_family->C,L_family->alpha);

  // alpha must be >=1/r so use 1/r
  arb_init(L_family->alpha);
  arb_set_d(L_family->alpha,1.0);
  arb_div_ui(L_family->alpha,L_family->alpha,L_family->r,prec);

  // c defined in Lemma 4  
  arb_init(L_family->c);
  arb_add(L_family->c,L_family->nu,L_family->alpha,prec);
  arb_set_d(tmp,0.5);
  arb_add(L_family->c,L_family->c,tmp,prec);

  // C defined in lemma 4
  if(L_family->r!=2)
    {
      printf("Value for C in lemma 4 hardwired for r=2. Exiting.\n");
      return(0);
    }
  arb_init(L_family->C);
  arb_sqrt_ui(L_family->C,3,prec);
  /*
  // c,c_dash defined in Lemma 5.4
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
  */
  L_family->one_over_B=read_gap(ffile);
  if(L_family->one_over_B==0.0)
    return(false);
  printf("B=%10.8e\n",1.0/L_family->one_over_B);
  arb_init(L_family->two_pi_by_B);
  arb_set_d(L_family->two_pi_by_B,L_family->one_over_B*2.0);
  arb_const_pi(tmp,prec);
  arb_mul(L_family->two_pi_by_B,L_family->two_pi_by_B,tmp,prec);

  printf("File has 2pi/B=");arb_printd(L_family->two_pi_by_B,10);printf("\n");

  L_family->low_i=read_low_i(ffile);
  if(L_family->low_i==BAD_64)
    return(false);
  L_family->hi_i=read_hi_i(ffile);
  if(L_family->hi_i==BAD_64)
    return(false);
  L_family->G_len=L_family->hi_i-L_family->low_i+1;

  printf("We have Gs from %ld to %lu inclusive.\n",L_family->low_i,L_family->hi_i);

  L_family->max_K=read_max_K(ffile);
  if(!L_family->max_K)
    return(false);
  Lc->K=L_family->max_K;
  printf("Number of differentials in file K = %lu\n",L_family->max_K);

  Lc->sks=(arb_t *)malloc(sizeof(arb_t)*MAX_M);
  for(uint64_t m=0;m<MAX_M;m++)
    arb_init(Lc->sks[m]);

  arb_init(Le->eq59);
  if(!read_eq59(Le->eq59,ffile))
    return(false);
  printf("Error for equation 5-9 (Taylor truncation) = ");
  arb_printd(Le->eq59,10);
  printf("\n");


  L_family_t *Lf=L_family;



  Lc->A=Lc->NN*Lf->one_over_B;
  printf("A=%10.8e\n",Lc->A);
  arb_init(Lc->arb_A);
  arb_set_d(Lc->arb_A,Lc->A);
  arb_init(Lc->one_over_A);
  arb_inv(Lc->one_over_A,Lc->arb_A,prec);
  // allocate vector of vectors to hold G^(k)(i/B)
  Lf->Gs=(arb_t **)malloc(sizeof(arb_t *)*Lc->K);
  for(uint64_t k=0;k<Lc->K;k++)
    {
      Lf->Gs[k]=(arb_t *)malloc(sizeof(arb_t)*Lf->G_len);
      for(uint64_t n=0;n<Lf->G_len;n++)
	arb_init(Lf->Gs[k][n]);
    }

  // now we read the G^(k)(2 pi i/B)
  // read them all for now, we'll select later
  if(!read_Gs(Lc,Lf,ffile,prec))
    return(false);

  Lc->kres=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
  Lc->skm=(acb_t *)malloc(sizeof(acb_t)*Lc->N);
  Lc->res=(acb_t *)malloc(sizeof(acb_t)*Lc->NN);
  for(uint64_t n=0;n<Lc->N;n++)
    {
      acb_init(Lc->kres[n]);
      acb_init(Lc->skm[n]);
    }
  for(uint64_t n=0;n<Lc->NN;n++)
    acb_init(Lc->res[n]);


  arb_clear(tmp);
  return(true);
}

bool init_func(L_func_t *L)
{
  arb_init(L->sum_ans);
  acb_init(L->epsilon);
  L->ans=(acb_t *)malloc(sizeof(acb_t)*MAX_M);
  if(!L->ans)
    {
      printf("Error allocating memory for a_n's. Exiting.\n");
      return(false);
    }
  for(uint64_t i=0;i<MAX_M;i++)
    acb_init(L->ans[i]);
  return(true);
}


// M values of a_n
bool read_func(L_func_t *L, L_comp_t *Lc, L_family_t *Lf, FILE *infile, double normalisation, uint64_t prec)
{
  printf("Conductor=%lu\n",L->N); 
  L->dc=sqrt((double) L->N);

  L->M=L->dc*exp(2*M_PI*(Lf->hi_i+0.5)*Lf->one_over_B);
  printf("M computed from hi_i = %lu\n",L->M);

  if(L->M>MAX_M)
    {
      printf("Can't handle M>%lu. Exiting.\n",MAX_M);
      return(false);
    }


  // read the M values a_n for n=1..M into ans[0]..ans[M-1]
  for(uint64_t m=0;m<L->M;m++)
    {
      if(!read_acb(L->ans[m],infile))
	{
	  printf("Error reading coefficient number %lu from file. Exiting.\n",m);
	  return(false);
	}
      //printf("Coeff[%lu] = ",m);acb_printd(L->ans[m],10);printf("\n");
    }
  fclose(infile);

  arb_t tmp,tmp1,n;
  arb_init(tmp);arb_init(tmp1);arb_init(n);

  // now divide each coefficient by the normalisation
  // we need n^1/2 for our algorithm
  // and n^normalisation for the rest
  arb_set_d(tmp,-normalisation-0.5);
  arb_zero(L->sum_ans);
  for(uint64_t i=1;i<L->M;i++)
    {
      arb_set_ui(n,i+1);
      arb_pow(tmp1,n,tmp,prec); // n^(-(k-1)/2-1/2)
      acb_mul_arb(L->ans[i],L->ans[i],tmp1,prec);
      acb_abs(tmp1,L->ans[i],prec);
      arb_add(L->sum_ans,L->sum_ans,tmp1,prec);
    }
  printf("sum |an/sqrt(n)|=");arb_printd(L->sum_ans,10);printf("\n");
  
  arb_clear(tmp);
  return(true);
}

int64_t calc_m(uint64_t i, double two_pi_by_B, double dc)
{
  double x=log((double) i/dc)/two_pi_by_B;
  return(round(x));
}


bool finalise_comp(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, FILE *ffile, uint64_t prec)
{
  double two_pi_by_B=Lf->one_over_B*2*M_PI;
  Lc->offset=calc_m(1,two_pi_by_B,Lfu->dc);
  if(Lc->offset<Lf->low_i)
    {
      printf("G values need to go down to %ld. We only have down to %ld.\n",
	     Lc->offset,Lf->low_i);
      return(false);
    }

  for(uint64_t m=1;m<=Lfu->M;m++)
    Lc->ms[m-1]=calc_m(m,two_pi_by_B,Lfu->dc);
  
  //for(uint64_t m=0;m<10;m++) printf("ms[%lu]=%ld\n",m,Lc->ms[m]);

  for(uint64_t m=0;m<Lfu->M;m++)
    arb_mul_si(Lc->ums[m],Lf->two_pi_by_B,Lc->ms[m],prec);
  
  /*
  for(uint64_t m=0;m<10;m++)
    {
      printf("ums[%lu]=",m);
      arb_printd(Lc->ums[m],10);
      printf("\n");
    }
  */
  arb_t tmp,tmp1,arb_sc;
  arb_init(tmp);arb_init(tmp1);arb_init(arb_sc);
  arb_sqrt_ui(tmp,Lfu->N,prec);
  arb_inv(arb_sc,tmp,prec);
  for(uint64_t m=0;m<Lfu->M;m++)
    {
      arb_mul_ui(tmp,arb_sc,m+1,prec);
      arb_log(tmp1,tmp,prec);
      arb_sub(Lc->sks[m],tmp1,Lc->ums[m],prec);
    }
  /*
  for(uint64_t m=0;m<10;m++)
    {
      printf("sks[0][%lu]=",m);
      arb_printd(Lc->sks[0][m],10);
      printf("\n");
    }
  */
  arb_clear(tmp);arb_clear(tmp1);arb_clear(arb_sc);
  /*
  for(uint64_t k=2;k<Lc->K;k++)
    for(uint64_t m=0;m<Lfu->M;m++)
      arb_pow_ui(Lc->sks[k-1][m],Lc->sks[0][m],k,prec);
  */

  return(true);
}  

bool set_errors(L_error_t *errs)
{
  arb_init(errs->lem54);
  arb_init(errs->lem56);
  arb_init(errs->lem57);
  //arb_init(errs->eq59);
  printf("Setting all other error terms to zero!.\n");

  arb_zero(errs->lem54);
  arb_zero(errs->lem56);
  arb_zero(errs->lem57);
  //arb_zero(errs->eq59);
  return(true);
}

// where does the entry for m+1 go?
uint64_t bucket(uint64_t m, L_comp_t *Lc, L_family_t *Lf)
{
  int64_t res=(-Lc->ms[m])%Lc->N;
  return(res);
}


bool do_convolves(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  acb_t tmp;arb_t tmp1;
  acb_init(tmp);
  arb_init(tmp1);
  int64_t n0=round(log(Lc->M0/Lfu->dc)/(Lf->one_over_B*2*M_PI));
  printf("Taking G values from %ld to %ld\n",n0,Lf->hi_i);
  /*
  for(uint64_t n=0;n<Lc->N;n++)
    acb_zero(Lc->res[n]);
  */
  // let's do k=0
  for(uint64_t n=0;n<Lc->N;n++)
    acb_zero(Lc->skm[n]);
  for(uint64_t m=Lc->M0;m<Lfu->M;m++)
    {
      uint64_t b=bucket(m,Lc,Lf);			
      acb_add(Lc->skm[b],Lc->skm[b],Lfu->ans[m],prec);
    }

  // just copy those G we actually need
  // i.e. from hi_i down to u_m=round(log(1/sqrt{conductor})*B/2/Pi)
  for(int64_t n=0;n<n0;n++)
    acb_zero(Lc->G[n]); 
  for(int64_t n=n0,n2=n0-Lf->low_i;n<=Lf->hi_i;n++,n2++)
    {
      int64_t n1=n%Lc->N;
      //printf("Putting G[%ld] in slot %ld.\n",n,n1);
      arb_set(acb_realref(Lc->G[n1]),Lf->Gs[0][n2]);
      arb_zero(acb_imagref(Lc->G[n1]));
    }
  for(int64_t n=Lf->hi_i+1;n<Lc->N;n++)
    {
      //printf("Zeroed G[%ld].\n",n);
      acb_zero(Lc->G[n]);
    }
  /*
    for(int64_t i=0;i<Lc->N;i++)
    {
    printf("G^[%lu]=",i);
    acb_printd(Lc->G[i],10);
    printf("\n");
    }
  */
  acb_convolve(Lc->res,Lc->skm,Lc->G,Lc->N,Lc->w,prec);
  //printf("G[0] after FFT = ");acb_printd(Lc->G[0],10);printf("\n");
  //printf("res[0]=");acb_printd(Lc->res[0],10);printf("\n");
  
  for(uint64_t k=1;k<Lc->K;k++)
    {
      for(uint64_t n=0;n<Lc->N;n++)
	acb_zero(Lc->skm[n]);
      for(uint64_t m=Lc->M0;m<Lfu->M;m++)
	{
	  uint64_t b=bucket(m,Lc,Lf);	
	  arb_pow_ui(tmp1,Lc->sks[m],k,prec);		
	  acb_mul_arb(tmp,Lfu->ans[m],tmp1,prec);
	  //acb_mul_arb(tmp,Lfu->ans[side][m],Lc->sks[k-1][m],prec);
	  acb_add(Lc->skm[b],Lc->skm[b],tmp,prec);
	  //printf("Put ");acb_printd(tmp,20);printf(" into bucket %lu\n",b);
	}
      
      for(int64_t n=0;n<n0;n++)
	acb_zero(Lc->G[n]); 
      for(int64_t n=n0,n2=n0-Lf->low_i;n<=Lf->hi_i;n++,n2++)
	{
	  int64_t n1=n%Lc->N;
	  //printf("Putting G[%ld] in slot %ld.\n",n,n1);
	  arb_set(acb_realref(Lc->G[n1]),Lf->Gs[k][n2]);
	  arb_zero(acb_imagref(Lc->G[n1]));
	}
      for(int64_t n=Lf->hi_i+1;n<Lc->N;n++)
	{
	  //printf("Zeroed G[%ld].\n",n);
	  acb_zero(Lc->G[n]);
	}

      /*
      for(int64_t n=n0,n2=n0-Lf->low_i;n<=Lf->hi_i;n++,n2++)
	{
	  int64_t n1=n%Lc->N;
	  //printf("Putting G[%ld] in slot %ld.\n",n,n1);
	  arb_set(acb_realref(Lc->G[n1]),Lf->Gs[k][n2]);
	  arb_zero(acb_imagref(Lc->G[n1]));
	}
      for(int64_t n=Lf->hi_i+1;n<Lc->N+n0;n++)
	{
	  //printf("Zeroed G[%ld].\n",n);
	  acb_zero(Lc->G[n]);
	}
      */

      acb_convolve(Lc->kres,Lc->skm,Lc->G,Lc->N,Lc->w,prec);
      //printf("G[%lu][0] after FFT = ",k);acb_printd(Lc->G[0],10);printf("\n");
      //printf("kres[0]=");acb_printd(Lc->kres[0],10);printf("\n");
      
      for(int64_t n=0;n<=Lc->N/2;n++)
	acb_add(Lc->res[n],Lc->res[n],Lc->kres[n],prec);      
    }
		    
  arb_clear(tmp1);
  acb_clear(tmp);
  return(true);
}

bool finish_convolves(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  acb_t tmp,tmp2;arb_t tmp1;
  acb_init(tmp);acb_init(tmp2);
  arb_init(tmp1);
  for(uint64_t k=0;k<Lc->K;k++)
    {
      for(uint64_t m=0;m<Lc->M0;m++)
	{
	  arb_pow_ui(tmp1,Lc->sks[m],k,prec);
	  acb_mul_arb(tmp,Lfu->ans[m],tmp1,prec); // an/sqrt(n)(log(m/sqrt(N))-um)^k
	  for(uint64_t n=0;;n++)
	    {
	      int64_t nn=Lc->ms[m]+n;
	      if(nn>Lf->hi_i) // run out of G values
		break;
	      acb_mul_arb(tmp2,tmp,Lf->Gs[k][nn-Lf->low_i],prec);
	      acb_add(Lc->res[n],Lc->res[n],tmp2,prec);
	    }
	}
    }
  /*
  for(uint64_t n=0;n<=Lc->N/2;n++)
    {
      printf("F_hat[%lu]=",n);
      acb_printd(Lc->res[n],10);
      printf("\n");
    }
  */
  acb_clear(tmp);acb_clear(tmp2);
  arb_clear(tmp1);
  return(true);
}



bool do_convolves_no_fft(L_comp_t *Lc, L_family_t *Lf, L_func_t *Lfu, uint64_t prec)
{
  acb_t tmp,tmp2;arb_t tmp1;
  acb_init(tmp);acb_init(tmp2);
  arb_init(tmp1);
  for(uint64_t n=0;n<Lc->N/2;n++)
    acb_zero(Lc->res[n]);

  for(uint64_t k=0;k<Lc->K;k++)
    {
      for(uint64_t m=0;m<Lfu->M;m++)
	{
	  arb_pow_ui(tmp1,Lc->sks[m],k,prec);
	  acb_mul_arb(tmp,Lfu->ans[m],tmp1,prec); // an/sqrt(n)(log(m/sqrt(N))-um)^k
	  for(uint64_t n=0;;n++)
	    {
	      int64_t nn=Lc->ms[m]+n;
	      if(nn>Lf->hi_i)
		break;
	      acb_mul_arb(tmp2,tmp,Lf->Gs[k][nn-Lf->low_i],prec);
	      acb_add(Lc->res[n],Lc->res[n],tmp2,prec);
	    }
	}
    }
  /*
  for(uint64_t n=0;n<=Lc->N/2;n++)
    {
      printf("F_hat[%lu]=",n);
      acb_printd(Lc->res[n],10);
      printf("\n");
    }
  */
  acb_clear(tmp);acb_clear(tmp2);
  arb_clear(tmp1);
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

bool final_ifft(L_comp_t *Lc, L_func_t *Lfu, L_family_t *Lf, L_error_t *Le, uint64_t prec)
{
  // first use res[side][0] and res[side][N/2] to figure our epsilon

  printf("Going to iFFT.\n");
  acb_ifft(Lc->res,Lc->NN,Lc->ww,prec);
  for(uint64_t n=0;n<Lc->NN;n++)
    acb_mul_arb(Lc->res[n],Lc->res[n],Lf->two_pi_by_B,prec);
  
  return(true);
}

void normalise(L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Lu, uint64_t prec)
{
  /*
  arb_t res,term;
  arb_t pi,one_over_A,t;
  acb_t s_plus_mu;

  acb_init(s_plus_mu);
  arb_init(res);arb_init(term);
  arb_init(pi);
  arb_const_pi(pi,prec);


  arb_init(t);
  arb_zero(t);


  for(uint64_t n=0;n<=Lc->N/2;n++,arb_add(t,t,one_over_A,prec))
    {
      arb_set_ui(res,1);
      for(uint64_t r=0;r<Lf->r;r++)
	{
	  arb_set_d(acb_realref(s_plus_mu),0.5+Lf->mus[r]);
	  arb_set(acb_imagref(s_plus_mu),t); // s=1/2+mu+it
	  //printf("s+mu=");acb_printd(s_plus_mu,10);printf("\n");
	  acb_mul_2exp_si(s_plus_mu,s_plus_mu,-1); // s/2
	  arb_pow(term,pi,acb_realref(s_plus_mu),prec); // |Pi^(s+mu)/2|
	  //printf("Pi^((s+mu)/2)=");arb_printd(term,10);printf("\n");
	  acb_mul_arb(Lc->res[n],Lc->res[n],term,prec);
	  acb_gamma(s_plus_mu,s_plus_mu,prec);
	  acb_abs(term,s_plus_mu,prec); // |GAMMA((s+mu)/2)|
	  //printf("|GAMMA((s+mu)/2)|=");arb_printd(term,10);printf("\n");
	  acb_div_arb(Lc->res[n],Lc->res[n],term,prec);
	}
    }
  */
  // Conrey upsampling assumes normalisation by exp(Pi t r/4)
  static arb_t pir4,one_over_A,t,tmp,tmp1;
  static bool init=false;
  if(!init)
    {
      arb_init(pir4);
      arb_const_pi(pir4,prec);
      arb_mul_ui(pir4,pir4,Lf->r,prec); // Pi r
      arb_mul_2exp_si(pir4,pir4,-2); // Pi r/4
      arb_init(one_over_A);
      arb_set_d(one_over_A,Lc->A);
      arb_inv(one_over_A,one_over_A,prec);
      arb_init(t);
      arb_init(tmp);arb_init(tmp1);
      init=true;
    }
  arb_zero(t);
  for(uint64_t n=0;n<=Lc->NN/2;n++,arb_add(t,t,one_over_A,prec))
    {
      arb_mul(tmp,pir4,t,prec);
      arb_exp(tmp1,tmp,prec);
      arb_mul(acb_realref(Lc->res[n]),acb_realref(Lc->res[n]),tmp1,prec);
      //arb_mul(acb_realref(Lc->res[Lc->N-n-1]),acb_realref(Lc->res[Lc->N-n-1]),tmp1,prec);
    }
  
  arb_set(t,one_over_A);
  for(uint64_t n=Lc->NN-1;n>Lc->NN/2;n--,arb_add(t,t,one_over_A,prec))
    {
      arb_mul(tmp,pir4,t,prec);
      arb_neg(tmp,tmp);
      arb_exp(tmp1,tmp,prec);
      arb_mul(acb_realref(Lc->res[n]),acb_realref(Lc->res[n]),tmp1,prec);
    }
  

}

bool setup_upsampling(L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, uint64_t prec)
{
  arb_init(Lu->H);
  arb_set_d(Lu->H,1.0/8.0);
  arb_init(Lu->inv_2H2); // -1/2H^2
  arb_set_d(Lu->inv_2H2,-32.0);
  Lu->N=70;
  Lu->stride=7; // use every 7th sample
  Lu->no_values=Lc->NN/OUTPUT_RATIO+Lc->NN/TURING_RATIO+Lu->N*2*Lu->stride;
  Lu->values[0]=(arb_t *)malloc(sizeof(arb_t)*Lu->no_values);
  Lu->values[1]=(arb_t *)malloc(sizeof(arb_t)*Lu->no_values);
  for(uint64_t n=0;n<Lu->no_values;n++)
    {
      arb_init(Lu->values[0][n]);
      arb_init(Lu->values[1][n]);
    }
  return(true);
}

// sin(pi x)/(pi x)
void my_sinc(arb_t res, arb_t x, uint64_t prec)
{
  static arb_t pi,tmp,tmp1;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(tmp);
      arb_init(tmp1);
    }
  arb_sin_pi(tmp,x,prec);
  arb_mul(tmp1,x,pi,prec);
  arb_div(res,tmp,tmp1,prec);
}



int64_t nearest_n(arb_ptr diff, arb_ptr t0, arb_t A, uint64_t prec)
{
  static arb_t tmp;
  static fmpz_t fmpz_tmp;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      fmpz_init(fmpz_tmp);
    }

  arb_mul(tmp,t0,A,prec);
  arb_floor(tmp,tmp,prec);
  if(arb_get_unique_fmpz(fmpz_tmp,tmp)==0)
    {
      printf("Error rounding to int in upsample routines. Exiting.\n");
      return(0);
    }
  int64_t res=fmpz_get_si(fmpz_tmp); // this is going to be 1st point in left tail
  arb_set_si(tmp,res);
  arb_div(tmp,tmp,A,prec);
  arb_sub(diff,tmp,t0,prec); // will be -ve n pts to left of t0
  return(res);
}

void do_point(arb_ptr res, L_upsample_t *Lu, L_comp_t *Lc, int64_t n, arb_t delta, arb_t A, uint64_t prec)
{
  static arb_t tmp,tmp1,tmp2;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
    }
  n=n%Lc->NN;
  arb_mul(tmp,delta,delta,prec);
  arb_mul(tmp1,tmp,Lu->inv_2H2,prec);
  arb_exp(tmp,tmp1,prec);
  //printf("exp term = ");arb_printd(tmp,10);printf("\n");
  arb_mul(tmp2,delta,A,prec);
  my_sinc(tmp1,tmp2,prec);
  //printf("sinc term = ");arb_printd(tmp1,10);printf("\n");
  arb_mul(tmp2,tmp1,tmp,prec);
  //printf("W term = ");arb_printd(acb_realref(Lc->res[n]),10);printf("\n");
  arb_mul(res,tmp2,acb_realref(Lc->res[n]),prec);
  //printf("Upsample contribution from %lu (delta = ",n);arb_printd(delta,10);
  //printf(" was ");arb_printd(res,10);printf("\n");
}

bool upsample_stride(arb_ptr res, arb_ptr t0, L_upsample_t *Lu, L_comp_t *Lc, L_family_t *Lf, uint64_t prec)
{
  static arb_t step,diff,this_diff,term,A;
  static bool init=false;
  if(!init)
    {
      init=true;
      arb_init(step);
      arb_mul_ui(step,Lc->one_over_A,Lu->stride,prec);
      arb_init(A);
      arb_div_ui(A,Lc->arb_A,Lu->stride,prec);
      arb_init(diff);
      arb_init(this_diff);
      arb_init(term);
    }
  //printf("upsampling with t0=");arb_printd(t0,100);printf("\n");
  int64_t n=nearest_n(diff,t0,Lc->arb_A,prec);
  //printf("nearest point is n=%lu\n",n);
  //printf("delta to t0 is ");arb_printd(diff,20);printf("\n");
  // do nearest point
  do_point(res,Lu,Lc,n,diff,A,prec);
  //printf("nearest point contributed ");arb_printd(res,20);printf("\n");
  // do points to left of nearest
  arb_sub(this_diff,diff,step,prec);
  for(int64_t nn=Lu->stride;nn<=Lu->N*Lu->stride;nn+=Lu->stride,arb_sub(this_diff,this_diff,step,prec))
    {
      //printf("doing point %lu with diff = ",n-nn);arb_printd(this_diff,20);printf("\n");
      do_point(term,Lu,Lc,n-nn,this_diff,A,prec);
      //printf("this point contributes ");arb_printd(term,20);printf("\n");
      arb_add(res,res,term,prec);
    }
  // do points to right of nearest
  arb_add(this_diff,diff,step,prec);
  for(int64_t nn=Lu->stride;nn<=Lu->N*Lu->stride;nn+=Lu->stride,arb_add(this_diff,this_diff,step,prec))
    {
      //printf("doing point %lu with diff = ",n+nn);arb_printd(this_diff,20);printf("\n");
      do_point(term,Lu,Lc,n+nn,this_diff,A,prec);
      //printf("this point contributes ");arb_printd(term,20);printf("\n");
      arb_add(res,res,term,prec);
    }
  //printf("upsampled value is ");arb_printd(res,100);printf("\n");
  //exit(0);
  return(true);
}

#define POS (1)
#define NEG (2)
#define UNK (0)  

uint8_t sign(arb_t x)
{
  if(arb_is_positive(x))
    return(POS);
  if(arb_is_negative(x))
    return(NEG);
  printf("Indeterminate value when upsampling. Exiting.\n");
  exit(0);
  return(UNK);
}

#define OP_ACC (101)


// use simple binary chop
bool isolate_zero(arb_t res, arb_t tt0, arb_t tt1, arb_t ff0, arb_t ff1, int8_t s0, int8_t s1, L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Ls, uint64_t prec)
{
  static bool init=false;
  static arb_t tmp1,tmp2,tmp3,t0,t1,f0,f1;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
      arb_init(t0);
      arb_init(t1);
      arb_init(f0);
      arb_init(f1);

    }
  arb_set(t0,tt0);
  arb_set(t1,tt1);
  arb_set(f0,ff0);
  arb_set(f1,ff1);
  while(true)
    {
      /*
      printf("t0 = ");arb_printd(t0,40);
      printf("\nt1 = ");arb_printd(t1,40);
      printf("\nf0 = ");arb_printd(f0,40);
      printf("\nf1 = ");arb_printd(f1,40);
      */
      arb_add(tmp1,t0,t1,prec);
      arb_mul_2exp_si(tmp1,tmp1,-1);
      printf("\nupsampling at ");arb_printd(tmp1,40);
      upsample_stride(tmp2,tmp1,Ls,Lc,Lf,prec);
      printf("\ngot ");arb_printd(tmp2,40);printf("\n");
      //exit(0);
      uint8_t new_sign=sign(tmp2);
      if(new_sign!=s0) // change is between s0 and tmp1
	{
	  arb_sub(tmp3,tmp1,t0,prec);
	  arb_mul_2exp_si(tmp3,tmp3,OP_ACC);
	  arb_sub_ui(tmp3,tmp3,1,prec);
	  if(arb_is_negative(tmp3))
	    {
	      arb_set(res,tmp1);
	      return(true);
	    }
	  arb_set(t1,tmp1);
	  arb_set(f1,tmp2);
	}
      else // change is between tmp1 and s1
	{
	  arb_sub(tmp3,t1,tmp1,prec);
	  arb_mul_2exp_si(tmp3,tmp3,OP_ACC);
	  arb_sub_ui(tmp3,tmp3,1,prec);
	  if(arb_is_negative(tmp3))
	    {
	      arb_set(res,tmp1);
	      return(true);
	    }
	  arb_set(t0,tmp1);
	  arb_set(f0,tmp2);
	}
    }
}
      
	  
bool find_zeros(L_comp_t *Lc, L_family_t *Lf, L_upsample_t *Ls, uint64_t prec)
{
  bool first=true;
  arb_t tmp,tmp1,tmp2,t0;
  arb_init(tmp);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(t0);

  uint64_t n=1;
  uint8_t last_sign=sign(acb_realref(Lc->res[0]));
  uint8_t this_sign;
  while(true)
    {
      this_sign=sign(acb_realref(Lc->res[n]));
      while(this_sign==last_sign)
	{
	  n++;
	  if(n==Lc->NN/8)
	    return(true);
	  this_sign=sign(acb_realref(Lc->res[n]));
	}
      printf("Zero found between %lu and %lu\n",n-1,n);
      arb_mul_ui(tmp1,Lc->one_over_A,n-1,prec);
      arb_mul_ui(tmp2,Lc->one_over_A,n,prec);
      if(!isolate_zero(tmp,tmp1,tmp2,acb_realref(Lc->res[n-1]),acb_realref(Lc->res[n]),last_sign,this_sign,Lc,Lf,Ls,prec))
	return(false);
      printf("Zero at ");
      arb_printd(tmp,100);
      printf("\n");
      if(first)
	first=false;
      else
	exit(0);
      last_sign=this_sign;
      n++;
    }

  return(true);
}


  

int main(int argc, char **argv)
{
  printf("Command line:- %s",argv[0]);
  for(uint64_t i=0;i<argc;i++)
    printf(" %s",argv[i]);
  printf("\n");
  if(argc!=7)
    {
      printf("Usage:- %s <prec> <comp file> <family file> <func file names>  <out file> <M0>.\n",argv[0]);
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
  FILE *fnamefile=fopen(argv[4],"r");
  if(!fnamefile)
    {
      printf("Failed to open L function file %s. Exiting.\n",argv[4]);
      return(0);
    }
  L_comp_t L_comp;
  L_comp.M0=atol(argv[6]);
  FILE *cfile;
  L_func_t L_func;
  L_error_t L_errors;
  L_upsample_t L_upsample;

  if(!read_comp(&L_comp,cfile,arb_prec))
    {
      printf("Error reading computation file. Exiting.\n");
      return(0);
    }
  
  
  if(!read_family(&L_family,&L_comp,&L_errors,ffile,arb_prec))
    {
      printf("Error reading L function family file. Exiting.\n");
      return(0);
    }
  
  //print_family(&L_family);


  if(!init_func(&L_func))
    {
      printf("Error initialising L-func structure. Exiting.\n");
      return(0);
    }

  if(!set_errors(&L_errors))
    {
      printf("Failed to set error terms. Exiting.\n");
      return(0);
    }

  if(!setup_upsampling(&L_upsample,&L_comp,&L_family,arb_prec))
    {
      printf("Error setting up for upsampling. Exiting.\n");
      return(0);
    }

  // iterate over a set of L functions in the same family
  char fname[1024];
  uint64_t wt;
      
  while(fscanf(fnamefile,"%lu %lu %s\n",&wt,&L_func.N,fname)==3)
    {
      printf("Doing L function described in %s.\n",fname);
      FILE *fufile=fopen(fname,"r");
      if(!fufile)
	{
	  printf("Failed to open L Function file at %s. Exiting.\n",fname);
	  return(0);
	}

      double normalisation=((double) wt-1.0)*0.5;
      printf("Normalisation for coefficients = n^(%10.8e)\n",-wt);

      if(!read_func(&L_func,&L_comp,&L_family,fufile,normalisation,arb_prec))
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
      if(!finish_convolves(&L_comp,&L_family,&L_func,arb_prec))
	{
	  printf("Error finishing convolutions. Exiting.\n");
	  return(0);
	}
      printf("Finished convolutions...\n");

      if(!do_pre_iFFT_errors(&L_comp,&L_family,&L_func,&L_errors,arb_prec))
	{
	  printf("Error doing pre-iFFT error terms. Exiting.\n");
	  return(0);
	}
      
      if(!final_ifft(&L_comp,&L_func,&L_family,&L_errors,arb_prec))
	{
	  printf("Error doing final IFFT. Exiting.\n");
	  return(0);
	}

      
      normalise(&L_comp,&L_family,&L_upsample,arb_prec);

      if(!find_zeros(&L_comp,&L_family,&L_upsample,arb_prec))
	{
	  printf("Error finding zeros. Exiting.\n");
	  return(0);
	}
      
    }
  // just do this for the last L function as a check
  
  for(int64_t n=L_comp.NN/2+1,nn=n-L_comp.NN;n<L_comp.NN;n++,nn++)
    {
      printf("Z( %10.8e )= ",(double) nn/L_comp.A);
      arb_printd(acb_realref(L_comp.res[n]),20);
      printf("\n");
    }   
  
  
  for(uint64_t n=0;n<=L_comp.NN/2;n++)
    {
      printf("Z( %10.8e )= ",(double) n/L_comp.A);
      arb_printd(acb_realref(L_comp.res[n]),20);
      printf("\n");
    }   
  
  return(0);
}
