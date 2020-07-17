// MS Stuff, does no harm elsewhere
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include "../includes/int_double10.0.h"
#include "../includes/im_s.h"

#define VERSION ""
// 6/5/9
// 1.0 Original



void print_usage()
/* called when wrong arguments passed via command line */
{
	printf("Usage: hur_fft%s (hur-file) (ofile)\n",VERSION);
	printf("  (hur-file)  - file with hurwitz values\n");
	printf("  (ofile)     - output file\n");
	exit(1);
}

void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
	std::cout << error_string << endl;
	exit(1);
}

inline int gcd (unsigned int a, unsigned int b)
  /* Euclid algorithm gcd */
{
  unsigned int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

// make a the smaller
inline int co_prime(unsigned int a, unsigned int b)
{
  return(gcd(a,b)==1);
};

inline unsigned int phi(unsigned int q)
{
	unsigned int res=0,a;
	for(a=1;a<q;a++)
		if(co_prime(a,q))
			res++;
	return(res);
}

int main(int argc, char **argv)
{
	unsigned int q_start,q_end,num_fracs,num_s,q,i,a,phi_q,phi_qs;
	FILE *hur_file,*out_file;
	im_s *im_s_vec;
	int_complex *out_vals;

	_fpu_rndd();

	clock_t no_clicks;

	no_clicks=clock(); // start timing

	if(argc!=3)
		print_usage();

	hur_file=fopen(argv[1],"rb");
	if(!hur_file)
		fatal_error("Couldn't open hurwitz data file. Exiting.\n");

	out_file=fopen(argv[2],"wb");
	if(!out_file)
		fatal_error("Couldn't open output file. Exiting.\n");

	fread(&q_start,sizeof(unsigned int),1,hur_file);  // lowest q in file
	fwrite(&q_start,sizeof(unsigned int),1,out_file);
	fread(&q_end,sizeof(unsigned int),1,hur_file);    // highest q in file
	fwrite(&q_end,sizeof(unsigned int),1,out_file);
	fread(&num_fracs,sizeof(unsigned int),1,hur_file);     // total a/q's in file
	fwrite(&num_fracs,sizeof(unsigned int),1,out_file);
	fread(&num_s,sizeof(unsigned int),1,hur_file);    // number of s values in file
	fwrite(&num_s,sizeof(unsigned int),1,out_file);

	if(!(im_s_vec=(im_s *) _aligned_malloc(num_s*sizeof(im_s),16)))
		fatal_error("Failed to allocate memory for im_s.\n");
	for(i=0;i<num_s;i++)
	{
		fread(&im_s_vec[i].im_s,sizeof(double),1,hur_file);
		fread(&im_s_vec[i].lambda_s,sizeof(int_complex),1,hur_file);
		fread(&im_s_vec[i].lambda_s_a,sizeof(int_complex),1,hur_file);
	}
	fwrite(im_s_vec,sizeof(im_s),num_s,out_file);
#ifndef LINUX
	_aligned_free(im_s_vec);
#endif

	printf("allocating memory\n");
	printf("we want %ld bytes\n",num_s*num_fracs*sizeof(int_complex));
	if(!(out_vals=(int_complex *) _aligned_malloc(num_s*num_fracs*sizeof(int_complex),16)))
		fatal_error("Failed to allocate memory for out_vals.\n");

	printf("reading out_vals\n");
	fread(out_vals,sizeof(int_complex),num_s*num_fracs,hur_file);

	printf("writing out_vals\n");
	phi_qs=0;
	for(q=q_start;q<=q_end;q++)
	{
		if((q&3)!=2)
		{
			phi_q=phi(q);
			for(i=0;i<num_s;i++)
				for(a=0;a<phi_q;a++)
					fwrite(&out_vals[i*num_fracs+phi_qs+a],sizeof(int_complex),1,out_file);
			phi_qs+=phi_q;
		}
	}

	return(0);
}


	



