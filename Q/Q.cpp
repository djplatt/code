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

int_double biggest=0.0;

void do_zero(int_double gam, uint64_t N)
{
  int_double T_2pi = gam/(2*d_pi);
  int_double q0=int_double(N-1)-T_2pi*(log(T_2pi)-1)-7.0/8.0;
  int_double q1=q0+1.0;
  if(q0.left<0)
    q0=-q0;
  if(q1.left<0)
    q1=-q1;
  int_double lg=log(gam);
  int_double r=q0/lg;
  //print_int_double_str("r = ",r);
  if(r>biggest)
    {
      printf("Just before zero %li\n",N);
      print_int_double_str("New biggest ",r);
      print_int_double_str("At gamma = ",gam);
      biggest=r;
    }
  r=q1/lg;
  if(r>biggest)
    {
      printf("Just after zero %lu\n",N);
      print_int_double_str("New biggest ",r);
      print_int_double_str("At gamma = ",gam);
      biggest=r;
    }
}


int main(int argc, char ** argv)
{
  _fpu_rndd();
  std::cout << "Command line was:- ";
  for(unsigned long int i=0;i<argc;i++)
    std::cout << argv[i] << " ";
  std::cout << std::endl;

  if(argc!=2)
    {
      printf("Usage:- %s <zeros file> \n", argv[0]);
      exit(0);
    }
  FILE*zfile=fopen(argv[1],"rb");
  if(!zfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }

  long int num_its;
  fread(&num_its,sizeof(long int),1,zfile);
  double st[2];
  long int zs[2];
  unsigned long int ptr=0;
  int_double del_t,t;
  double OP_ACC=1.0;
  for(long int i=0;i<102;i++)
    OP_ACC/=2.0;
  std::cout << "OP_ACC set to " << OP_ACC << std::endl;
  int_double err=int_double(-OP_ACC,OP_ACC);
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
	  del_t=in_bytes(zfile);
	  if(contains_zero(del_t))
	    {
	      printf("Two zeros 0 apart. Exiting.\n");
	      exit(0);
	    }
	  t+=del_t;
	  int_double tt=t+err;
	  do_zero(tt,z);
	}
    }
  return(0);
}

