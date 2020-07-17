//
// zetalim200.cpp
//
// DJ Platt Jan 2014
//
// Find upb on |zeta(1/2+it)|
// for t in [0,200]
#include "inttypes.h"
#include "../includes/int_double12.0.h"

//
// simple version of Euler Mac for Zeta(s) with sum of length N
//
int_complex zeta(const int_complex &s, uint64_t N)
{
  int_complex res(1);
  for(uint64_t n=2;n<=N;n++)
    res+=pow(n,-s);
  res+=pow(N,-s+1)/(s-1)-pow(N,-s)/2;
  int_complex term=s*pow(N,-s-1)/12; // k=1 term
  res+=term;
  term=s*(s+1)*(s+2)*pow(N,-s-3)/720; // k=2 term
  res-=term;
  int_double factor=mod(s+3)/(s.real+3);
  int_double error=factor*mod(term);
  error.left=error.right;
  res.real+=error;
  res.imag+=error;
  return(res);
} 

int main(int argc, char **argv)
{
  printf("Command Line:- %s ",argv[0]);
  for(uint64_t i=1;i<argc;i++)
    printf("%s ",argv[i]);
  printf("\n");
  if(argc!=4)
    {
      printf("Usage:- %s <delta num> <delta den> <outfile>.\n",argv[0]);
      exit(0);
    }

  // set up the SSE system for rounding and initialise constants
  _fpu_rndd();

  // this is how much we step by each time
  double full_delta=atof(argv[1])/atof(argv[2]);
  double ts=0.0;
  double te=200.0;

  FILE *outfile=fopen(argv[3],"wb");
  if(!outfile)
    {
      printf("Failed to open file %s for binary output. Exiting.\n",argv[3]);
      exit(0);
    }

  double largest=0.0;
  while(ts<te)
    {
      int_double t,res;
      int_complex cres;
      t=int_double(ts,ts+full_delta);
      cres=zeta(int_complex(int_double(0.5),t),100+floor(t.left));
      res=mod(cres);
      if(-res.right>largest)
	{
	  print_int_double_str("zeta=",res);
	  largest=-res.right;
	  fwrite(&t.left,sizeof(double),1,outfile);
	  fwrite(&largest,sizeof(double),1,outfile);
	  
	  //print_int_double_str("",res);
	  printf("new largest at %60.58e = %60.58e\n",t.left,largest);
	}
      ts+=full_delta;
      //printf("ts=%30.28e\n",ts);
    }
  return(0);
}
