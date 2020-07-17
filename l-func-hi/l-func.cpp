//
// l-func.cpp
//
// hi precision, rigorous computation of Dirichlet L-functions
//
// computes \Lambda_\chi(1/2_it)=\omega (q/\pi)^{it}\Gamma((1/2+it+a_\chi)/2)\exp(Pi t/4)L_\chi(1/2+it)
// where a_chi=(1-\chi(-1))/2
// omega=
//
#include <iostream>
#include <cstdlib>
#include "characters.h"

using namespace std;

#define PREC (300)

mpfi_c_t ctemp,ctemp1;
mpfi_t itemp,itemp2,mpfi_2_pi;

// read a complex from ifile
// stored as mpz_t left right exponent for real and imaginary parts
void read_complex(mpfi_c_t res, FILE *infile, mpz_ptr l, mpz_ptr r)
{
  if(!mpz_inp_raw(l,infile))
    {
      printf("Error reading from hurwitz file. Exiting.\n");
      exit(0);
    }
  if(!mpz_inp_raw(r,infile))
    {
      printf("Error reading from hurwitz file. Exiting.\n");
      exit(0);
    }
  mpfi_interv_z(res->re,l,r);
  if(!mpz_inp_raw(l,infile))
    {
      printf("Error reading from hurwitz file. Exiting.\n");
      exit(0);
    }
  int64_t exp=mpz_get_si(l);
  mpfi_mul_2si(res->re,res->re,exp);
  
  if(!mpz_inp_raw(l,infile))
    {
      printf("Error reading from hurwitz file. Exiting.\n");
      exit(0);
    }
  if(!mpz_inp_raw(r,infile))
    {
      printf("Error reading from hurwitz file. Exiting.\n");
      exit(0);
    }
  mpfi_interv_z(res->im,l,r);
  if(!mpz_inp_raw(l,infile))
    {
      printf("Error reading from hurwitz file. Exiting.\n");
      exit(0);
    }
  exp=mpz_get_si(l);
  mpfi_mul_2si(res->im,res->im,exp);
}

void write_real(mpfi_ptr x, FILE *outfile, mpz_ptr l, mpfr_ptr m)
{
  mpfr_exp_t ex;
  mpfi_get_left(m,x);
  ex=mpfr_get_z_2exp(l,m);
  if(!mpz_out_raw(outfile,l))
    {
      printf("Failed to write mantissa in write-real. Exiting.\n");
      exit(0);
    }
  if(fwrite(&ex,sizeof(mpfr_exp_t),1,outfile)!=1)
    {
      printf("Failed to write exponent in write-real. Exiting.\n");
      exit(0);
    }
  mpfi_get_right(m,x);
  ex=mpfr_get_z_2exp(l,m);
  if(!mpz_out_raw(outfile,l))
    {
      printf("Failed to write mantissa in write-real. Exiting.\n");
      exit(0);
    }
  if(fwrite(&ex,sizeof(mpfr_exp_t),1,outfile)!=1)
    {
      printf("Failed to write exponent in write-real. Exiting.\n");
      exit(0);
    }
}


// set s_vec[i]=s(s+1)...(s+i)/(i+1)!
// so s_vec[0]=s/1!
//    s_vec[1]=s(s+1)/2!
//    s_vec[2]=s(s+1)(s+2)/3!
void init_s_vec(double t, double gap, uint64_t N_COLUMNS, mpfi_c_t *s_vec)
{
  mpfi_set_d(s_vec[0]->re,0.5);
  mpfi_set_d(s_vec[0]->im,t); // temp=s=1/2+it
  mpfi_c_set(ctemp,s_vec[0]);
  for(uint64_t i=1;i<N_COLUMNS;i++)
    {
      mpfi_add_ui(ctemp->re,ctemp->re,1);
      mpfi_c_mul(s_vec[i],s_vec[i-1],ctemp);
      mpfi_c_div_ui(s_vec[i],s_vec[i],i+1);
    }
}

// pretend s_vec[-1]=1
// zeta(s,M+alpha)=sum_{i=0}^N_COLUMNS zeta(s+i,M+alpha_0)s_vec[i-1](alpha_0-alpha)^i
void estimate_hurwitz(mpfi_c_t res, uint64_t a, uint64_t q, uint64_t N_COLUMNS, uint64_t N_ROWS, double gap, mpfi_c_t *hur_vec, mpfi_c_t *s_vec, mpfi_ptr alpha)
{
  double drow=((double) a/(double(q)))/gap+0.5;
  uint64_t row=drow;
  mpfi_c_t *this_row=hur_vec+N_COLUMNS*row;
  //mpfi_c_print_str("this_row[0]=",this_row[0]);
  mpfi_sub_d(itemp,alpha,row*gap); // row*gap is where zeta is evaluated at
  mpfi_neg(itemp,itemp);
  //printf("a=%lu q=%lu row=%lu ",a,q,row);mpfi_print_str("delta=",itemp);
  uint64_t i=N_COLUMNS-2;
  mpfi_c_mul(res,this_row[i+1],s_vec[i]);
  for(;i>0;i--)
    {
      mpfi_c_mul_i(res,res,itemp);
      mpfi_c_mul(ctemp,this_row[i],s_vec[i-1]);
      mpfi_c_add(res,res,ctemp);
    }
  mpfi_c_mul_i(res,res,itemp);
  mpfi_c_add(res,res,this_row[0]);
}

void do_omegas(mpfi_c_t *omegas, mpfi_c_t *in_vec, uint64_t q, DirichletGroup *G)
{
  mpfi_div_ui(itemp,mpfi_2_pi,q);
  //printf("In do_omegas.\n");fflush(stdout);
  for(uint64_t i=1;i<q;i++)
    if(G->is_coprime_to_q(i))
      {
	mpfi_mul_ui(itemp2,itemp,i);
	mpfi_sin(in_vec[i]->im,itemp2);
	mpfi_cos(in_vec[i]->re,itemp2);
      }
  G->DFTsum(omegas,in_vec);
  mpfi_set_ui(itemp,1);
  mpfi_div_ui(itemp2,itemp,q);
  mpfi_sqrt(itemp,itemp2);
  for(uint64_t i=2;i<q;i++)
    if(G->character(i).is_primitive())
      {
	mpfi_c_mul_i(omegas[i],omegas[i],itemp); // div by sqrt(q)
	if(G->character(i).is_even())
	  mpfi_neg(omegas[i]->im,omegas[i]->im); // conj()
	else
	  mpfi_swap(omegas[i]->re,omegas[i]->im); // 
	mpfi_c_sqrt(omegas[i],omegas[i]);
      }
  //printf("Exiting do_omegas.\n");fflush(stdout);
  
}

int main(int argc, char** argv)
{
  if(argc!=4)
    {
      printf("Usage:- %s <q> <hurwitz file> <outfile>. Exiting.\n",argv[0]);
      exit(0);
    }
 
  uint64_t q=atol(argv[1]);
  if((q&3)==2)
    {
      printf("moduli congruent to 2 mod 4 have no primitive characters. Exiting.\n");
      exit(0);
    }

  FILE *infile=fopen(argv[2],"rb");
  if(!infile)
    {
      printf("Failed to open %s for binary input. Exiting.\n",argv[2]);
      exit(0);
    }
  FILE *outfile=fopen(argv[3],"wb");
  if(!outfile)
    {
      printf("Failed to open %s for binary output. Exiting.\n",argv[3]);
      exit(0);
    }

  if(fwrite(&q,sizeof(uint64_t),1,outfile)!=1)
    {
      printf("Error writing q. Exiting.\n");
      exit(0);
    }

  uint64_t N_COLUMNS,N_ROWS,M;
  double t_start,t_end,del_t;
  if(fread(&N_COLUMNS,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading number of columns. Exiting.\n");
      exit(0);
    }
  if(fread(&N_ROWS,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading number of rows. Exiting.\n");
      exit(0);
    }
  if(fread(&M,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Error reading M. Exiting.\n");
      exit(0);
    }
  if(fread(&t_start,sizeof(double),1,infile)!=1)
    {
      printf("Error reading t_start. Exiting.\n");
      exit(0);
    }
  if(fwrite(&t_start,sizeof(double),1,outfile)!=1)
    {
      printf("Error writing t_start. Exiting.\n");
      exit(0);
    }
  if(fread(&t_end,sizeof(double),1,infile)!=1)
    {
      printf("Error reading t_end. Exiting.\n");
      exit(0);
    }
  if(fwrite(&t_end,sizeof(double),1,outfile)!=1)
    {
      printf("Error writing t_end. Exiting.\n");
      exit(0);
    }
  
  // check parameters and set up Taylor error
  if((t_end>220.0)||(N_ROWS<512)||(N_COLUMNS<50)||(M!=2))
    {
      printf("Taylor error is hard wired for t<=220, delta<=1/1024, N>=50 and M=2. Exiting.\n");
      exit(0);
    }

  if(fread(&del_t,sizeof(double),1,infile)!=1)
    {
      printf("Error reading del_t. Exiting.\n");
      exit(0);
    }
  if(fwrite(&del_t,sizeof(double),1,outfile)!=1)
    {
      printf("Error writing del_t. Exiting.\n");
      exit(0);
    }





  // how many individual values of s to process
  uint64_t num_s=(t_end-t_start)/del_t;

  // initialise bits
  dft_init(PREC);
  mpfi_c_init(ctemp);
  mpfi_c_init(ctemp1);
  mpfi_init(itemp);
  mpfi_init(itemp2);
  mpfi_init(mpfi_2_pi);
  mpfi_const_pi(mpfi_2_pi);
  mpfi_mul_ui(mpfi_2_pi,mpfi_2_pi,2);
  mpfi_c_t tay_err;
  mpfi_c_init(tay_err);
  mpfi_set_ui(tay_err->re,1);
  mpfi_put_si(tay_err->re,-1);
  mpfi_div_2ui(tay_err->re,tay_err->re,300); // error is <2^-300
  mpfi_set(tay_err->im,tay_err->re);

  // allocate memory
  // hur_vec holds the lattice of Hurwitz zeta values for a single s
  mpfi_c_t *hur_vec= (mpfi_c_t *)malloc(sizeof(mpfi_c_t)*N_COLUMNS*(N_ROWS+1));
  if(!hur_vec)
    {
      printf("Error allocating memory for hur_vec. Exiting.\n");
      exit(0);
    }
  // 
  mpfi_c_t *in_vec=(mpfi_c_t *)malloc(q*sizeof(mpfi_c_t));
  mpfi_c_t *out_vec=(mpfi_c_t *)malloc(q*sizeof(mpfi_c_t));
  mpfi_c_t *omegas=(mpfi_c_t *)malloc(q*sizeof(mpfi_c_t));

  if(!in_vec||!out_vec||!omegas)
    {
      printf("Error allocating memory for in_vec/out_vec/omegas. Exiting.\n");
      exit(0);
    }
  mpfi_c_t *s_vec=(mpfi_c_t *)malloc(N_COLUMNS*sizeof(mpfi_c_t));
  if(!s_vec)
    {
      printf("Error allocating memory for s_vec. Exiting.\n");
      exit(0);
    }

  for(uint64_t i=0;i<q;i++)
    {
      mpfi_c_init(in_vec[i]);
      mpfi_c_init(out_vec[i]);
      mpfi_c_init(omegas[i]);
    }
  for(uint64_t i=0;i<N_COLUMNS;i++)
    mpfi_c_init(s_vec[i]);
  for(uint64_t row=0,ptr=0;row<=N_ROWS;row++)
    for(uint64_t column=0;column<N_COLUMNS;column++)
      mpfi_c_init(hur_vec[ptr++]);

  mpz_t l,r;
  mpz_init(l);mpz_init(r);

  DirichletGroup G(q);
  mpfi_t **results;
  results=(mpfi_t **)malloc(sizeof(mpfi_t *)*q);
  if(!results)
    {
      printf("Error allocating memory for results. Exiting.\n");
      exit(0);
    }
  bool *is_prim,*is_ev;
  is_prim=(bool *)malloc(sizeof(bool)*q);
  is_ev=(bool *)malloc(sizeof(bool)*q);
  if(!is_prim||!is_ev)
    {
      printf("Error allocating memory for is_prim/is_ev. Exiting.\n");
      exit(0);
    }
  for(uint64_t i=1;i<q;i++)
    {
      if(G.character(i).is_primitive())
	{
	  is_prim[i]=true;
	  results[i]=(mpfi_t *)malloc(sizeof(mpfi_t)*num_s);
	  if(!(results[i]))
	    {
	      printf("Error allocating memory for results[%lu]. Exiting.\n",i);
	      exit(0);
	    }
	  for(uint64_t j=0;j<num_s;j++)
	    mpfi_init(results[i][j]);
	}
      else
	is_prim[i]=false;
      if(G.character(i).is_even())
	is_ev[i]=true;
      else
	is_ev[i]=false;
    }
  // set up the root numbers, used to compute Lambda from L
  do_omegas(omegas,in_vec,q,&G);

  mpfi_t alpha,lnqpi;
  mpfi_c_t minus_s,lng,lnga;
  mpfi_init(alpha);
  mpfi_init(lnqpi);
  mpfi_c_init(minus_s);
  mpfi_c_init(lng);
  mpfi_c_init(lnga);
  mpfi_set_d(minus_s->re,-0.5);
  mpfi_set_ui(lnqpi,q*2);
  mpfi_div(lnqpi,lnqpi,mpfi_2_pi);
  mpfi_log(lnqpi,lnqpi);
  mpfr_t diam;
  mpfr_init(diam);


  mpfr_t mrtemp;
  mpfr_init(mrtemp);

  double t=t_start,gap=(double)1.0/(double)N_ROWS;
  for(uint64_t s=0,ptr=0;s<num_s;s++,ptr=0,t+=del_t)
    {
      read_complex(lng,infile,l,r); // loggamma((1/2+it)/2)+Pi*t/4;
      read_complex(lnga,infile,l,r); // loggamma((3/2+it)/2)+Pi*t/4;
      //mpfi_c_print_str("lngamma=",lng);
      //mpfi_c_print_str("lngamma_a=",lnga);

      mpfi_mul_d(itemp,lnqpi,t/2.0); // t/2 log(q/pi)
      mpfi_add(lng->im,lng->im,itemp);
      mpfi_add(lnga->im,lnga->im,itemp);
      mpfi_c_exp(lng,lng);
      mpfi_c_exp(lnga,lnga);

      mpfi_set_d(minus_s->im,-t);
      init_s_vec(t,gap,N_COLUMNS,s_vec);
      for(uint64_t row=0;row<=N_ROWS;row++)
	for(uint64_t column=0;column<N_COLUMNS;column++)
	  {
	    //mpfi_c_print_str("",hur_vec[ptr]);
	    read_complex(hur_vec[ptr],infile,l,r);
	    ptr++;
	  }
      for(uint64_t i=1;i<q;i++)
	{
	  //printf("Checking coprimality.\n");fflush(stdout);
	  if(G.is_coprime_to_q(i))
	    {
	      //printf("In coprime case.\n");fflush(stdout);
	      mpfi_set_ui(itemp,i);
	      mpfi_div_ui(alpha,itemp,q);
	      estimate_hurwitz(in_vec[i],
			       i,q,N_COLUMNS,N_ROWS,gap,hur_vec,s_vec,alpha);
	      mpfi_c_add(in_vec[i],in_vec[i],tay_err);
	      for(uint64_t m=0;m<M;m++)
		{
		  mpfi_c_pow_i_c(ctemp,alpha,minus_s);
		  mpfi_c_add(in_vec[i],in_vec[i],ctemp);
		  mpfi_add_ui(alpha,alpha,1);
		}
	    }
	}
      //printf("About to DFT.\n");fflush(stdout);
      G.DFTsum(out_vec,in_vec);
      //printf("Completed DFT.\n");fflush(stdout);

      mpfi_set_ui(itemp,q);
      mpfi_c_pow_i_c(ctemp,itemp,minus_s);
      for(uint64_t i=1;i<q;i++)
	if(is_prim[i])
	  {
	    mpfi_c_mul(out_vec[i],ctemp,out_vec[i]);
	    //printf("L_{%lu,%lu}(1/2+i%f)=",i,q,t);mpfi_c_print_str("",out_vec[i]);
	    mpfi_c_mul(ctemp1,out_vec[i],omegas[i]);
	    if(is_ev[i])
	      mpfi_c_mul(out_vec[i],ctemp1,lng);
	    else
	      mpfi_c_mul(out_vec[i],ctemp1,lnga);
	    if(!mpfi_contains_zero(out_vec[i]->im))
	      {
		printf("Lambda was not real. Exiting.\n");
		printf("Lambda_{%lu,%lu}(1/2+i%f)=",i,q,t);
		mpfi_c_print_str("",out_vec[i]);
		exit(0);
	      }
	    mpfi_set(results[i][s],out_vec[i]->re);
	  }
    }
  for(uint64_t i=1;i<q;i++)
    if(is_prim[i])
      {
	uint64_t inverse=InvMod(i,q);
	//printf("inverse of %lu is %lu\n",i,inverse);
	if(inverse<i) continue;
	fwrite(&i,sizeof(uint64_t),1,outfile);
	for(uint64_t s=0;s<num_s;s++)
	  write_real(results[i][s],outfile,l,mrtemp);
	if(inverse!=i) // not a real character
	  {
	    fwrite(&inverse,sizeof(uint64_t),1,outfile);
	    for(uint64_t s=0;s<num_s;s++)
	      write_real(results[inverse][s],outfile,l,mrtemp);
	  }
      }
  return(0);
}

