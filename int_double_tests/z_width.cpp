// z_width.cpp
// look at widths of intervals in CONVERTED fft_output
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
#include "../includes/int_double11.0.h"
#include "../includes/im_s.h"

#define VERSION ""
// 29/7/9
// 1.0 Original



void print_usage()
/* called when wrong arguments passed via command line */
{
	printf("Usage: z_width%s (fft-file)\n",VERSION);
	printf("  (fft-file)  - file from CONVERTED int-l-func9.9 routine\n");
	exit(0);
}

void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
	std::cout << error_string << " Exiting." << endl;
	exit(0);
}

double get_width(int_double x)
{
	return(-(x.left+x.right));
};

int main(int argc, char **argv)
{
	unsigned int num_s,q;
	unsigned int s,num_prims,prim,index;
	int_complex omega;
	FILE *hur_file;
	im_s im_st;
	int_double z,z1;
	double width,max_width;
	bool neg_one;

	_fpu_rndd();

	clock_t no_clicks;

	no_clicks=clock(); // start timing

	if(argc!=2)
		print_usage();

	hur_file=fopen(argv[1],"rb");
	if(!hur_file)
		fatal_error("Couldn't open fft data file. Exiting.\n");
	fread(&num_s,sizeof(unsigned int),1,hur_file);
	printf("num_s=%d\n",num_s);
	for(s=0;s<num_s;s++)
	{
		fread(&im_st,sizeof(im_s),1,hur_file);
		//printf("%10.8e\n",im_st.im_s);
	}
	while(fread(&q,sizeof(unsigned int),1,hur_file))
	{
		printf("q=%d\n",q);
		fread(&num_prims,sizeof(unsigned int),1,hur_file);
		for(prim=0;prim<num_prims;prim++)
		{
			fread(&index,sizeof(unsigned int),1,hur_file);
			fread(&omega,sizeof(int_complex),1,hur_file);
			fread(&neg_one,sizeof(bool),1,hur_file);
			max_width=0.0;
			for(s=0;s<num_s;s++)
			{
				fread(&z,sizeof(int_double),1,hur_file);
				if(contains_zero(z))
					print_int_double_str("z spanning zero ",z);
				width=get_width(z);
				if(width>max_width)
				{
					max_width=width;
					z1=z;
				}
			}
		printf("q=%d:%d max width found was %10.8e for z=[%10.8e,%10.8e]\n",q,index,max_width,z1.left,-z1.right);
		}
	}


	return(0);
}


	



