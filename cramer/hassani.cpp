#include "../includes/int_double15.0.h"

// read our 13 byte structure representing a zero gap
// into an int_double
int_double in_bytes(FILE *infile)
{
  uint64_t a;uint32_t b;uint8_t c;
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
  double r1,r2;
  //printf("Read a=%lX b=%x c=%x\n",a,b,c);
  r1=(double) c/32.0 + (double) b/(32.0*256.0*256.0*256.0*256.0)+(double) (a&0xfff8000000000000)/(32.0*65536.0*65536.0*65536.0*65536.0*65536.0*65536.0); // /2^3, /2^37, /2^101
  r2=(double) (a&0x0007ffffffffffff)/(32.0*65536.0*65536.0*65536.0*65536.0*65536.0*65536.0); // /2^101
  //printf("r1=%50.48e\nr2=%50.48e\n",r1,r2);

  return(int_double(r1)+r2);
}

int_double recip_sum(int_double *gammas, uint64_t N)
{
  int_double res=0.0;
  for(int64_t n=N-1;n>=0;n--)
    res+=1.0/gammas[n];
  return res;
}


int main(int argc, char ** argv)
{
  _fpu_rndd();
  std::cout << "Command line was:- ";
  for(unsigned long int i=0;i<argc;i++)
    std::cout << argv[i] << " ";
  std::cout << std::endl;

  if(argc!=5)
    {
      printf("Usage:- %s <zeros file list> <N> <T num> <T den>\n",argv[0]);
      exit(0);
    }
  FILE *lfile=fopen(argv[1],"rb");
  if(!lfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  uint64_t N=atol(argv[2]);
  int_double *zeros=(int_double *)malloc(sizeof(int_double)*N);
  if(zeros == NULL)
    {
      printf("Failed to allocate memeory for zeros. Exiting.\n");
      exit(0);
    }
  int_double T=atol(argv[3]);
  T/=atol(argv[4]);
  print_int_double_str("T (not checked against N) set to ",T);
  uint64_t zz=0;
  double OP_ACC=1.0;
  for(long int i=0;i<102;i++)
    OP_ACC/=2.0;
  int_double err=int_double(-OP_ACC,OP_ACC);
  char fname[1024];
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
      int_double del_t,t;
      //std::cout << "OP_ACC set to " << OP_ACC << std::endl;
      for(long int it=0;it<num_its;it++)
	{
	  //if((it%100)==0)
	  //printf("Starting block %lu/%lu\n",it+1,num_its);
	  fread(st,sizeof(double),2,zfile);
	  fread(&zs[0],sizeof(long int),1,zfile);
	  if(st[0]==0.0)
	    continue;
	  t=st[0];
	  fread(&zs[1],sizeof(long int),1,zfile);
	  printf("Processing zero %ld to %ld=%ld in total.\n",zs[0]+1,zs[1],zs[1]-zs[0]);
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
	      del_t=in_bytes(zfile);
	      if(contains_zero(del_t))
		{
		  printf("Two zeros 0 apart. Exiting.\n");
		  exit(0);
		}
	      t+=del_t;
	      zeros[zz-1]=t+err;
	    }
	  if(zz>N)
	    break;
	}
      fclose(zfile);
    }

  //printf("zz=%lu\n",zz);
  
  if(zz<N)
    {
      printf("Ran out of zeros after %lu\n",zz);
      exit(0);
    }
  
  printf("The %lu'th zero is ",N);
  print_int_double_str("",zeros[N-1]);
  int_double res=recip_sum(zeros,N),c1,four_pi;
  c1=d_two_pi*exp(int_double(1.0));
  four_pi=2.0*d_two_pi;
  print_int_double_str("sum 1/gamma = ",res);
  res=res-N/T-(sqr(log(T/c1))+1)/four_pi+(7.0/8.0)/T;
  print_int_double_str("H estimate = ",res);
  return 0;
}

