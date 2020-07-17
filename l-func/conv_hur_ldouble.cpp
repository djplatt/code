#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <complex>

using namespace std;

typedef complex<long double> dcomplex;

void print_usage()
  /* called when wrong arguments passed via command line */
{
	printf("Usage: test (ifile)\n");
	printf("  (ifile) - file with dcomplex vals of Z(s,a/q)\n");
	exit(0);
};

void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(0);
};

ifstream infile;

int main(int argv, char **argc)
{

  dcomplex z;
  long double x;
  int i,j,k,num_fracs;

  if(argv!=2)
    print_usage();
  infile.open(argc[1],ios::in|ios::binary);
  if(!infile.is_open())
    fatal_error("Couldn't open file. Exiting.");
  infile.read((char *) &i,sizeof(int));
  cout << "q_start = " << i << endl;
  infile.read((char *) &i,sizeof(int));
  cout << "q_end = " << i << endl;
  infile.read((char *)&num_fracs,sizeof(int));
  cout << "num_fracs = " << num_fracs << endl;
  infile.read((char *)&i,sizeof(int));
  cout << "num s points = " << i << endl;
  for(k=0;k<i;k++)
    {
      infile.read((char *) &x,sizeof(long double));
      printf("s=0.5+%20.18ei\n",(double)x);
      for(j=0;j<num_fracs;j++)
	{
	  infile.read((char *) &z,sizeof(dcomplex));
	  printf("%20.18e+%20.18ei\n",(double)real(z),(double)imag(z));
	}
    }

  infile.close();
  return(0);
};
