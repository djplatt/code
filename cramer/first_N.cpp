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

int_double f(int_double g1, int_double g2, int_double lg2)
{
  int_double s,c;
  sin_cos(lg2*(g1-g2),&s,&c);
  return 4*sqrt((sqr(4*c-1)+16*sqr(s))/((4*sqr(g1)+1)*(4*sqr(g2)+1)*(4+sqr(g1-g2))));
  //return mod((int_double(4.0)*exp(int_complex(zero,lg2*(g1-g2)))-1)/(int_complex(half,g1)*int_complex(half,-g2)*int_complex(two,g1-g2)));
}

int_double cramer_sum(int_double *gammas, uint64_t N)
{
  int_double res=0.0,lg2=log(int_double(2.0));
  for(uint64_t n=0;n<N;n++)
    {
      int_double g1=gammas[n];
      res+=2*(f(g1,g1,lg2)+f(g1,-g1,lg2));
      for(uint64_t m=n+1;m<N;m++)
	{
	  int_double g2=gammas[m];
	  res+=4*(f(g1,g2,lg2)+f(g1,-g2,lg2));
	  //res+=2.0*f(g1,g2,lg2)+2*f(g1,-g2,lg2);//+f(-g1,g2,lg2);//+f(-g1,-g2,lg2);
	}
    }
  return res;
}


int main(int argc, char ** argv)
{
  _fpu_rndd();
  std::cout << "Command line was:- ";
  for(unsigned long int i=0;i<argc;i++)
    std::cout << argv[i] << " ";
  std::cout << std::endl;

  if(argc!=3)
    {
      printf("Usage:- %s <zeros file list> <N>\n",argv[0]);
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
  uint64_t zz=0;
  double OP_ACC=1.0;
  for(long int i=0;i<102;i++)
    OP_ACC/=2.0;
  int_double err=int_double(-OP_ACC,OP_ACC);
  char fname[1024];
  while(fscanf(lfile,"%s\n",fname)==1)
    {
      FILE *zfile=fopen(fname,"rb");
      if(!zfile)
	{
	  printf("Failed to open file %s for binary input. Exiting.\n",fname);
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
	  //printf("Processing zero %ld to %ld=%ld in total.\n",zs[0]+1,zs[1],zs[1]-zs[0]);
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
      if(zz>N)
	break;
    }
  if(zz!=N+1)
    {
      printf("Ran out of zeros at %lu\n",zz);
      exit(0);
    }
  printf("The %lu'th zero is ",N);
  print_int_double_str("",zeros[N-1]);
  int_double res=cramer_sum(zeros,N);
  print_int_double_str("Sum = ",res);
  return 0;
}

