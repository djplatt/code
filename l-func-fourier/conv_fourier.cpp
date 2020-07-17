#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
//#include <complex>


using namespace std;

//typedef complex<double> dcomplex;

void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: conv_fourier (ifile)\n");
	printf("  (ifile) - file with dcomplex vals of Lambda(s,chi)\n");
	exit(0);
};

void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(0);
};

ifstream infile;

int main(int argc, char **argv)
{

//  dcomplex z;
  double *im_s,x;
  int i,q,j,num_s,s;
  FILE *infile;

  if(argc!=2)
    print_usage();
  if(!(infile=fopen(argv[1],"rb")))
	  fatal_error("Couldn't open file. Exiting.");
  fread(&num_s,sizeof(int),1,infile);
  im_s=(double *) malloc(sizeof(double)*num_s);
  cout << "number of s values = " << num_s << endl;
  for(j=0;j<num_s;j++)
    fread(&im_s[j],sizeof(double),1,infile);
  while(!feof(infile))
    {
      fread(&q,sizeof(int),1,infile);
      if(feof(infile))
		  break;
      cout << "q= " << q << endl;
      fread(&j,sizeof(int),1,infile);
      if(feof(infile))
		  break;
      cout << "j = " << j << endl;
      for(s=0;s<num_s;s++)
	{
	  fread(&x,sizeof(double),1,infile);
	  if(feof(infile))
	    break;
	  printf("%10.4e: %20.18e\n",im_s[s],x);
	}
    }
  fclose(infile);
  return(0);
};
