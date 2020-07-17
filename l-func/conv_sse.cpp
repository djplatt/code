#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;


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

  double z[4];
  int i;

  if(argv!=2)
    print_usage();
  infile.open(argc[1],ios::in|ios::binary);
  if(!infile.is_open())
    fatal_error("Couldn't open file. Exiting.");
  infile.read((char *) &i,sizeof(int));
  cout << "q_start = " << i << endl;
  infile.read((char *) &i,sizeof(int));
  cout << "q_end = " << i << endl;
  infile.read((char *)&i,sizeof(int));
  cout << "num_fracs = " << i << endl;
  while(true)
    {
      infile.read((char *) &z,4*sizeof(double));
      if(infile.eof())
	break;
      printf("[%20.18e,%20.18e]+[%20.18e,%20.18e]i\n",
	     z[1],-z[0],z[3],-z[2]);
      printf("\n");
    };

  infile.close();
  return(0);
};
