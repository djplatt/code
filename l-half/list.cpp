#include "stdio.h"
#include "math.h"

int main(int argc, char **argv)
{
  if(argc!=2)
    return(0);
  FILE *infile=fopen(argv[1],"rb");
  if(!infile)
    return(0);
  unsigned long int q;
  fread(&q,sizeof(unsigned long int),1,infile);
  printf("q=%lu\n",q);
  double sq=sqrt(1.0/(double) q);
  double z[2];
  while(true)
    {
      if(fread(z,sizeof(double),2,infile)!=2)
	return(0);
      printf("z=%20.18e %20.18e\n",z[0]*sq,z[1]*sq);
    }
  return(0);
}
