#include "pari.h"
#include "math.h"

#define bitp ((unsigned long int)(1<<22)-1)

void print_usage()
{
  exit(0);
}

int liouville(GEN n)
{
  int om=bigomega(n);
  if(om&1)
    return(-1);
  else
    return(1);
}

int main(int argc, char **argv)
{
  long int n0,n1,n;
  double res=0.0,high_point=0.0;
  GEN pn;
  if(argc!=3)
    print_usage();
  n0=atoi(argv[1]);
  n1=atoi(argv[2]);
  if((n1<n0)||(n0<=0))
    print_usage();
  pari_init(1<<20,1<<20);
  pn=stoi(n0);
  for(n=n0;n<=n1;n++)
    {
      res=res+liouville(pn)/sqrt(n);
      if(res>high_point)
	high_point=res;
      gaddgsz(pn,1,pn);
    }
  printf("for      n=%ld-%ld high point was %20.18e\n",n0,n1,high_point);
  printf("sum over n=%ld-%ld high point was %20.18e\n",n0,n1,res);
  return(0);
}
