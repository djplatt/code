#include <stdio.h>
#include <stdlib.h>
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

typedef struct{double a;double b;double c;double d;} int_complex;

void print_int_complex(int_complex &z)
{
	printf("[%20.18e,%20.18e]+i[%20.18e,%20.18e]",z.a,-z.b,z.c,-z.d);
}

int main(int argv, char **argc)
{

  int_complex z;
  double *im_s;
  int i,j,num_fracs,num_s;

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
  infile.read((char *) &num_s,sizeof(int));
  cout << "num_s = " << num_s << endl;
  if(!(im_s=(double *)malloc(num_s*sizeof(double))))
	  fatal_error("Couldn't allocate memory for im_s. Exiting.");
  for(i=0;i<num_s;i++)
	  infile.read((char *) &im_s[i],sizeof(double));
  for(i=0;i<num_s;i++)
  {
	  cout << "s=0.5+i" << im_s[i] << endl;
	  for(j=0;j<num_fracs;j++)
	  {
		  infile.read((char *) &z,sizeof(int_complex));
		  printf("%5d ",j);
		  print_int_complex(z);
		  printf("\n");
	  }
  }

  infile.close();
  return(0);
};
