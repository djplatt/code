#include "../includes/int_double12.0.h"

#define N_START (70368744177664L)
#define PI_START (2280998753949L)
#define N_STEP (2097152L)// 4194304L)
#define STEPS_PER_FILE (3145728L)// 1572864L)
#define N_FILES (160L)

#define LI_TERMS (127)
#define LOG_TERMS (7)
int_double li_terms[LI_TERMS];
int_double logns[7];

bool li_init=true;

void init_li()
{
  li_terms[0]=int_double(1.0);
  for(int i=1;i<LI_TERMS;i++)
    li_terms[i]=li_terms[i-1]/(i+1);
  for(int i=1;i<LI_TERMS;i++)
    li_terms[i]/=(i+1);
}

int_double comp_bin(long unsigned int n)
{
  int_double res=1.0;
  //printf("comp_bin called with %lu\n",n);
  for(int i=0;i<LOG_TERMS;i++)
    {
      if(n&1)
	res*=logns[i];
      n>>=1;
    }
  //print_int_double_str("comp_bin returning ",res);
  return(res);
}

long unsigned int li_count=0;

//
// first missing term < log(2^50)^129/129!/129
// common ratio < log(2^50)/130
// => error in [0,1e-21]
int_double li_err=int_double(0,1e-21);

int_double li(unsigned long int n)
{
  if(li_init)
    {
      init_li();
      li_init=false;
      //for(int i =0;i<10;i++) print_int_double_str("li_terms=",li_terms[i]);
    }
  //printf("Li called with %lu\n",n);
  li_count++;
  int_double nn=n,res=int_double(0.5772,0.5773);
  logns[0]=log(nn);
  res+=log(logns[0]);
  for(int i=1;i<LOG_TERMS;i++)
    logns[i]=sqr(logns[i-1]);

  for(int i=0;i<LI_TERMS;i++)
    res+=comp_bin(i+1)*li_terms[i];

  //print_int_double_str("Li returning (put truncation error term in Dave)",res);exit(0);
  return(res+li_err);
}

int main(int argc,char **argv)
{
  _fpu_rndd();

  int_double tmp,li_step,li_n,diff,logn,n_step;
  double low_diff=1e100,this_diff;
  n_step=int_double(N_STEP);
  char fname[256];
  FILE *infile;
  long unsigned int from,to,primes,n=N_START,pi_n=PI_START;
  li_n=li(n);
  logn=log(int_double(n));
  for(long unsigned int file_no=0;file_no<N_FILES;file_no++)
    {
      sprintf(fname,"sieve_%lu_%lu",file_no*STEPS_PER_FILE,STEPS_PER_FILE);
      printf("Processing file %lu of %lu named %s\n",file_no+1,N_FILES,fname);
      infile=fopen(fname,"r");
      if(!infile)
	{
	  printf("Failed to open file for input. Exiting.\n");
	  exit(0);
	}
      for(long unsigned int line=0;line<STEPS_PER_FILE;line++)
	{
	  diff=(li_n-pi_n)*logn;
	  //print_int_double_str("Diff=",diff);
	  if(diff.left<=N_STEP) // recompute Li
	    {
	      li_n=li(n);
	      diff=(li_n-pi_n)*logn;
	      if(diff.left<=N_STEP) // still fails
		{
		  /*
		  this_diff=(diff.left-diff.right)/2.0;
		  if(this_diff<low_diff)
		    low_diff=this_diff;
		  */
		  printf("Li test failed\n");
		  printf("n=%lu\npi(n)=%lu\n",n,pi_n);
		  print_int_double_str("Li(n)=",li_n);
		  print_int_double_str("diff*log(n)=",diff);
		  exit(0);
		  
		}
	    }
	  fscanf(infile,"Primes in %lu to %lu %lu\n",&from,&to,&primes);
	  if(from!=(n+1))
	    {
	      printf("Data mismatch. exiting\n");
	      exit(0);
	    }
	  tmp=n_step/logn;
	  li_step.right=tmp.right;
	  n+=N_STEP;
	  logn=log(int_double(n));
	  pi_n+=primes;
	  tmp=n_step/logn;
	  li_step.left=tmp.left;
	  li_n+=li_step;
	}
      fclose(infile);
    }
  printf("Final n was %lu\nFinal pi(n) was %lu\n",n,pi_n);
  print_int_double_str("Final li(n) was ",li_n);
  print_int_double_str("Final difference was ",diff);
  printf("Li was called %lu times\n",li_count);

  //printf("Lowest Difference found was %20.18e\n",low_diff);
  return(0);
}
