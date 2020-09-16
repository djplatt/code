#include "quad.h"

#define OP_ACC (101)

#define true (0==0)
#define false (1==0)


// read a 13 byte number from file
// structured 8,4,1
// read as if its exact
void in_bytes(arb_ptr t, FILE *infile, int64_t prec)
{
  uint64_t a;
  uint32_t b;
  uint8_t c;
  int res;

  if(fread(&a,sizeof(uint64_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (a). Exiting.\n");
      exit(0);
    }
  //printf("a=%lu\n",a);
  if(fread(&b,sizeof(uint32_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (b). Exiting.\n");
      exit(0);
    }
  //printf("b=%u\n",b);
  if(fread(&c,sizeof(uint8_t),1,infile)!=1)
    {
      printf("Fatal error in in_bytes (c). Exiting.\n");
      exit(0);
    }
  arb_set_ui(t,c);
  arb_mul_2exp_si(t,t,32);
  arb_add_ui(t,t,b,prec);
  arb_mul_2exp_si(t,t,64);
  arb_add_ui(t,t,a,prec);
  arb_mul_2exp_si(t,t,-OP_ACC);
}

#define MAXD_STEPS (40)
#define MOLIN_N (40)

void acb_theta(acb_t res, const acb_t t, int64_t prec)
{
  static int init=false;
  static arb_t quarter;
  static acb_t z,tmp,tmp1;
  static arb_t logpi;
  if(!init)
    {
      init=true;
      arb_init(quarter);
      arb_set_d(quarter,0.25);
      acb_init(z);
      acb_init(tmp);
      acb_init(tmp1);
      arb_set_d(acb_realref(z),0.25);
      arb_init(logpi);
      arb_const_pi(logpi,prec);
      arb_log(logpi,logpi,prec);
    }
  acb_mul_onei(z,t);
  acb_mul_2exp_si(z,z,-1); // it/2
  arb_add(acb_realref(z),acb_realref(z),quarter,prec); // z=1/4+it/2
  acb_lgamma(tmp,z,prec);
  acb_mul_arb(tmp1,t,logpi,prec);
  acb_mul_2exp_si(tmp1,tmp1,-1); // t/2 log pi
  acb_sub_arb(res,tmp1,acb_imagref(tmp),prec);
  acb_neg(res,res);
}


void arb_theta(arb_t res, const arb_t t, int64_t prec)
{
  static int init=false;
  static acb_t z,tmp;
  static arb_t logpi,tmp1;
  if(!init)
    {
      init=true;
      acb_init(z);
      acb_init(tmp);
      arb_init(tmp1);
      arb_set_d(acb_realref(z),0.25);
      arb_init(logpi);
      arb_const_pi(logpi,prec);
      arb_log(logpi,logpi,prec);
    }
  arb_mul_2exp_si(acb_imagref(z),t,-1); // z=1/4+it/2
  acb_lgamma(tmp,z,prec);
  arb_mul(tmp1,t,logpi,prec);
  arb_mul_2exp_si(tmp1,tmp1,-1); // t/2 log pi
  arb_sub(res,acb_imagref(tmp),tmp1,prec);
}

// 1/pi int_T^{T_eps} theta(t) dt
// need T, T_eps close enough so that when rescaled to -1,1,
// the pole of lngamma(1/4+it*/2) at t=i/2 is outside |2| 
void do_th_eps(arb_t res, arb_t T, arb_t T_eps, int64_t prec)
{
  static int init=false;
  static arb_t pi,md,tmp;
  if(!init)
    {
      init=true;
      arb_init(pi);
      arb_const_pi(pi,prec);
      arb_init(md);
      arb_init(tmp);
    }
  arb_maxd(md,acb_theta,T,T_eps,prec,MAXD_STEPS);
  molin_int(tmp,MOLIN_N,arb_theta,md,T,T_eps,prec);
  arb_div(res,tmp,pi,prec);
  //printf("Int ");arb_printd(T,10);printf(" to ");arb_printd(T_eps,10);printf(" returning ");arb_printd(res,20);printf("\n");
}
/*
#define T_MAX (600)
#define N (341)
*/

#define T_MAX (530) // bigger than 168 pi
#define N (291)

/*
#define T_MAX (50)
#define N (10)
*/

uint64_t do_N1_eps(arb_t res, arb_t N1, arb_t T_eps, arb_t *zeros, arb_t eps, uint64_t n, int64_t prec)
{
  static int init=false;
  static arb_t tmp1;
  if(!init)
    {
      init=true;
      arb_init(tmp1);
    }
  //printf("Doing N1 with n=%lu\nN1=",n);
  //arb_printd(N1,20);printf("\nT+eps=");
  //arb_printd(T_eps,20);
  //printf("\n");

  arb_mul_ui(tmp1,eps,n,prec);
  arb_add(res,tmp1,N1,prec); // N1 assuming no zeros in [T,T+eps]
  while(arb_lt(zeros[n],T_eps))
    {
      arb_sub(tmp1,T_eps,zeros[n],prec);
      printf("adding ");arb_printd(tmp1,20);printf(" to N1.\n");
      arb_add(res,res,tmp1,prec);
      n++;
      if(n==N)
	return n;
    }
  if(!arb_gt(zeros[n],T_eps))
    {
      printf("n=%d\nzero at ",n);
      arb_printd(zeros[n],20);printf("\n");
      printf("Zero straddled integration point. Giving up.\n");
      exit(0);
    }
  //printf("N1 returning ");arb_printd(res,20);printf("\n");
  return n;
}

// 2.067 + 0.059 log t  
void a_bit(arb_t res, arb_t t, int64_t prec)
{
  static int init=false;
  static arb_t A0,A1,tmp;
  if(!init)
    {
      init=true;
      arb_init(A0);
      arb_init(A1);
      arb_init(tmp);
      arb_set_ui(A0,2067);
      arb_div_ui(A0,A0,1000,prec); //2.067
      arb_set_ui(A1,59);
      arb_div_ui(A1,A1,1000,prec); // 0.059
    }
  arb_log(res,t,prec);
  arb_mul(tmp,res,A1,prec);
  arb_add(res,tmp,A0,prec);
}

// check |S1(t)|<A0+A1 log t for t in [6,T_MAX]
void do_s1(arb_t *zeros, int64_t prec)
{
  arb_t eps; // step size
  arb_init(eps);
  arb_set_d(eps,1.0/256.0); // small enough for T_MAX=600
  arb_t T,N1,th_int,T_eps,th_eps,N1_eps,tmp1,tmp2,tmp3,tmp4,zero,ab;
  arb_init(T);
  arb_set_ui(T,0);
  arb_init(N1); // int_0^T N(t) dt
  arb_init(th_int); // 1/pi int_6^T theta(t) dt
  arb_init(T_eps);
  arb_init(th_eps);
  arb_init(N1_eps);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  arb_init(tmp4);
  arb_init(zero);
  arb_init(ab);
  uint64_t n=0; // N(T)
  int first_fail=true;
  arb_const_pi(tmp1,prec);
  arb_mul_ui(tmp4,tmp1,168,prec); // 168 pi
  n=do_N1_eps(N1,zero,tmp4,zeros,tmp4,n,prec);
  uint64_t i;
  arb_t tmp5;
  arb_init(tmp5);
  arb_mul_2exp_si(tmp3,tmp4,-13);
  for(i=0;i<1024*8;i++)
    {
      arb_mul_ui(tmp1,tmp3,i,prec);
      arb_mul_ui(tmp2,tmp3,i+1,prec);
      do_th_eps(tmp5,tmp1,tmp2,prec);
      arb_add(th_int,th_int,tmp5,prec);
    }
  printf("Finished with T = ");arb_printd(tmp4,10);
  printf("\n1/Pi Int_0^T theta(t) dt = ");arb_printd(th_int,10);
  printf("\nN1(T) = ");arb_printd(N1,10);
  arb_sub(tmp1,N1,tmp4,prec);
  arb_sub(tmp2,tmp1,th_int,prec);
  printf("\nS1(T) = ");arb_printd(tmp2,20);printf("\n");
  return;
}


int main(int argc, char **argv)
{
  printf("Command line was:- ");
  for(unsigned long int i=0;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");

  if(argc!=2)
    {
      printf("Usage:- %s <prec>\n",argv[0]);
      exit(0);
    }
  FILE *lfile=fopen("S1_zeros.lst","rb");
  if(!lfile)
    {
      printf("Failed to open file S1_zeros.lst for binary input. Exiting.\n");
      exit(0);
    }
  arb_t T;
  arb_init(T);
  int64_t prec=atol(argv[1]);
  uint64_t z;
  arb_t zeros[N+1];
  for(z=0;z<=N;z++)
    arb_init(zeros[z]);
  arb_set_ui(zeros[N],T_MAX*2);
  char fname[1024];
  arb_t del_t,t;
  arb_init(t);arb_init(del_t);
  arb_t z_err;
  arb_set_ui(z_err,1);
  arb_mul_2exp_si(z_err,z_err,-OP_ACC-1);
  arb_t res,tmp,tmp1;
  arb_init(res);arb_init(tmp);arb_init(tmp1);
  uint64_t zz=0;  
  while((zz<N)&&(fscanf(lfile,"%s\n",fname)==1))
    {
      printf("Processing file %s\n",fname);
      FILE *zfile=fopen(fname,"rb");
      if(zfile == NULL)
	{
	  perror("Exiting: ");
	  exit(0);
	}
      long int num_its;
      fread(&num_its,sizeof(long int),1,zfile);
      double st[2];
      long int zs[2];
      unsigned long int ptr=0;
      
      for(long int it=0;it<num_its;it++)
	{
	  fread(st,sizeof(double),2,zfile);
	  fread(&zs[0],sizeof(long int),1,zfile);
	  if(st[0]==0.0)
	    continue;
	  arb_set_d(t,st[0]);
	  fread(&zs[1],sizeof(long int),1,zfile);
	  for(long int z=zs[0]+1;z<=zs[1];z++)
	    {
	      zz++;
	      if(zz!=z)
		{
		  printf("zz=%lu z=%lu\n",zz,z);
		  exit(0);
		}
	      if(zz>N)
		break;
	      in_bytes(del_t,zfile,prec);
	      if(arb_contains_zero(del_t))
		{
		  printf("Two zeros 0 apart. Exiting.\n");
		  exit(0);
		}
	      arb_add(t,t,del_t,prec);
	      arb_set(zeros[zz-1],t);
	      arb_add_error(zeros[zz-1],z_err);
	    }
	  if(zz>N)
	    break;
	}
      fclose(zfile);
    }

  
  if(zz<N)
    {
      printf("Ran out of zeros after %lu\n",zz);
      exit(0);
    }

  do_s1(zeros,prec);

  return 0;
}
