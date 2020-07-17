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

int_double al,max_ghat;

// g1(t)=sin(al*t)/(al*t)
int_double g1(int_double gam)
{
  static bool init=false;
  if(!init)
    {
      init=true;
      al=0.5*log(int_double(2));
      max_ghat=1.0/log(int_double(2));
      print_int_double_str("alpha = ",al);
      print_int_double_str("max ghat = ",max_ghat);
    }
  int_double s[1],c[1];
  sin_cos((al*gam),s,c);
  return s[0]/(al*gam);
}

// g2(t)=sin^2(al*t)/(al*t)^2
int_double g2a(int_double gam)
{
  static bool init=false;
  if(!init)
    {
      init=true;
      al=0.25*log(int_double(2));
      max_ghat=2.0/log(int_double(2));
      print_int_double_str("alpha = ",al);
      print_int_double_str("max ghat = ",max_ghat);
    }
  int_double s[1],c[1];
  sin_cos((al*gam),s,c);
  return sqr(s[0]/(al*gam));
}

//g3(t)=sin^3(al*t)/(al*t)^3
int_double g3a(int_double gam)
{
  static bool init=false;
  if(!init)
    {
      init=true;
      al=log(int_double(2))/6.0;
      max_ghat=9.0/(4.0*log(int_double(2)));
      print_int_double_str("alpha = ",al);
      print_int_double_str("max ghat = ",max_ghat);
    }
  int_double s[1],c[1];
  sin_cos((al*gam),s,c);
  int_double tmp=s[0]/(al*gam);
  return tmp*tmp*tmp;
}

double B,C;
//g4(t)=sin^4(al*t)/(al*t)^4
int_double g4a(int_double gam)
{
  static bool init=false;
  if(!init)
    {
      init=true;
      int_double l2=log(int_double(2));
      al=l2/8.0;
      if(B==0.0)
	max_ghat=8.0/(3*l2);
      else
	{
	  if(C==0.0)
	    max_ghat=(8.0*B-8)/(3.0*B*l2);
	  else
	    max_ghat=8.0/(3*log(int_double(2)))-8.0/(3.0*B*l2)-8.0/(3.0*C*l2)+384.0/(3.0*B*C*l2*l2*l2);
	}
      print_int_double_str("alpha = ",al);
      print_int_double_str("max ghat = ",max_ghat);
    }
  int_double s[1],c[1];
  sin_cos((al*gam),s,c);
  int_double tmp=s[0]/(al*gam);
  return sqr(sqr(tmp));
}


//g3(t)=sin^3(al*t)/(al*t)^3*(1-gam/B)
int_double g3(int_double gam)
{
  if(B>0.0)
    return g3a(gam)*(1.0-gam/B);
  else
    return g3a(gam);
}

int_double g2(int_double gam)
{
  if(B>0.0)
    return g2a(gam)*(1.0-gam/B);
  else
    return g2a(gam);
}

int_double g4(int_double gam)
{
  if(B==0.0)
    return g4a(gam);
  if(C>0.0)
    return g4a(gam)*(1.0-gam/B)*(1.0-gam/C);
  return g4a(gam)*(1.0-gam/B);
}


// not quite right at max |ghat(t)| is not quite at t=0
int_double g5a(int_double gam)
{
  static bool init=false;
  if(!init)
    {
      init=true;
      int_double lg2=log(int_double(2));
      al=lg2/10.0;
      max_ghat=575.0/(192.0*lg2);
      print_int_double_str("alpha = ",al);
      print_int_double_str("max ghat = ",max_ghat);
    }
  int_double s[1],c[1];
  sin_cos((al*gam),s,c);
  int_double tmp=s[0]/(al*gam);
  return tmp*sqr(sqr(tmp));
}

int_double g5(int_double gam)
{
  if(B==0.0)
    return g5a(gam);
  if(C>0.0)
    return g5a(gam)*(1.0-gam/B)*(1.0-gam/C);
  else
    return g5a(gam)*(1.0-gam/B);
}

/*
int_double g6_2(int_double gam)
{
  static bool init=false;
  if(!init)
    {
      init=true;
      int_double lg2=log(int_double(2));
      al=lg2/12.0;
      max_ghat=33.0/(10.0*lg2)+216.0/(lg2*lg2*lg2*B*C);
      print_int_double_str("alpha = ",al);
      print_int_double_str("max ghat = ",max_ghat);
    }
  int_double s[1],c[1];
  sin_cos((al*gam),s,c);
  int_double tmp=s[0]/(al*gam);
  if(B==0.0)
    return sqr(sqr(tmp)*tmp);
  if(C>0.0)
    return sqr(sqr(tmp)*tmp)*(1.0-gam/B)*(1.0-gam/C);
return sqr(sqr(tmp)*tmp)*(1.0-gam/B);
}
*/


int_double g(int_double gam)
{
  return abs(g3(gam));
}

int main(int argc, char ** argv)
{
  _fpu_rndd();
  std::cout << "Command line was:- ";
  for(unsigned long int i=0;i<argc;i++)
    std::cout << argv[i] << " ";
  std::cout << std::endl;

  if(argc!=6)
    {
      printf("Usage:- rs <zeros file> <B num> <B den> <C num> <C den>\n");
      exit(0);
    }
  FILE*zfile=fopen(argv[1],"rb");
  if(!zfile)
    {
      printf("Failed to open file %s for binary input. Exiting.\n",argv[1]);
      exit(0);
    }
  B=atof(argv[2]);
  if(B>0)
    B/=atof(argv[3]);
  C=atof(argv[4]);
  if(C>0)
    C/=atof(argv[5]);

  long int num_its;
  fread(&num_its,sizeof(long int),1,zfile);
  double st[2];
  long int zs[2];
  unsigned long int ptr=0;
  int_double del_t,t,s1,s2,r,r1,t1;
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
	  r=sqrt(0.25+tt*tt);
	  if(z==1)
	    {
	      print_int_double_str("gamma_1 = ",tt);
	      t1=tt; // gamma_1
	      r1=sqrt(0.25+t1*t1); // |rho_1|
	    }
	  else
	    s1+=abs(g(tt-t1))/r; // add |g(gamma_n-gamma_1)|/|rho_n|
	  s2+=abs(g(-tt-t1))/r; // add |g(-gamma_n-gamma_1)|/|rho_n|
	}
    }
  print_int_double_str("1/|rho_1= ",1.0/r1);
  print_int_double_str("s1= ",s1);
  print_int_double_str("s2= ",s2);
  printf("B=%e C=%e ",B,C);
  print_int_double_str("J1_5000(X)/max_ghat>= ",(1.0/r1-s1-s2)/max_ghat);
  print_int_double_str("J1_5000(X)>= ",(1.0/r1-s1-s2));
  return(0);
}

