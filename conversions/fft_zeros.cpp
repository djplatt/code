// fft_zeros.cpp
// convert file from output of fft to input for zeros.cpp
// MS Stuff, does no harm elsewhere
// edited 10/8/2009
//
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"

#define VERSION ""
// 7/5/9
// 1.0 Original



void print_usage()
/* called when wrong arguments passed via command line */
{
	printf("Usage: fft_zeros%s (hur-file) (ofile)\n",VERSION);
	printf("  (hur-file)  - file from fft routine\n");
	printf("  (ofile)     - output file for zeros.cpp\n");
	exit(1);
}

void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
	std::cout << error_string << " Exiting." << endl;
	exit(1);
}


unsigned int conj_j (unsigned int j, unsigned int num_chi, unsigned int q)
{
	if(q&7) // q not divisible by 8
	  return(num_chi-j-1);
	if(j<(num_chi>>1))
	  return((num_chi>>1)-j-1);
	return(num_chi+(num_chi>>1)-j-1);
}


int main(int argc, char **argv)
{
	unsigned int q_start,q_end,num_s,q;
	unsigned int s,num_prims,prim,*indices,q_file,prim1;
	FILE *hur_file,*out_file;
	//im_s *im_s_vec;
	int_complex *omegas;
	int_double *zs;
	bool *neg_ones;

	_fpu_rndd();

	clock_t no_clicks;

	no_clicks=clock(); // start timing

	if(argc!=3)
		print_usage();

	hur_file=fopen(argv[1],"rb");
	if(!hur_file)
		fatal_error("Couldn't open fft data file. Exiting.\n");

	out_file=fopen(argv[2],"wb");
	if(!out_file)
		fatal_error("Couldn't open output file. Exiting.\n");

	fread(&q_start,sizeof(unsigned int),1,hur_file);  // lowest q in file
//	fwrite(&q_start,sizeof(unsigned int),1,out_file);
	fread(&q_end,sizeof(unsigned int),1,hur_file);    // highest q in file
//	fwrite(&q_end,sizeof(unsigned int),1,out_file);
	fread(&num_prims,sizeof(unsigned int),1,hur_file); // max number of primitives we will encounter
	fread(&num_s,sizeof(unsigned int),1,hur_file);    // number of s values in file
	fwrite(&num_s,sizeof(unsigned int),1,out_file);

	//printf("q_start=%d\nq_end=%d\nnum_prims=%d\nnum_s=%d\n",q_start,q_end,num_prims,num_s);

	//if(!(im_s_vec=(im_s *) _aligned_malloc(num_s*sizeof(im_s),16)))
	//fatal_error("Failed to allocate memory for im_s.\n");

	//fread(im_s_vec,sizeof(im_s),num_s,hur_file);
	//fwrite(im_s_vec,sizeof(im_s),num_s,out_file);
	//#ifndef LINUX
	//_aligned_free(im_s_vec);
	//#endif

	//printf("allocating memory\n");

	if(!(indices=(unsigned int *) malloc(num_prims*sizeof(unsigned int))))
		fatal_error("Couldn't allocate memory for indices.\n");
	if(!(omegas=(int_complex *) _aligned_malloc(num_prims*sizeof(int_complex),16)))
		fatal_error("Couldn't allocate memory for omegas.\n");
	if(!(neg_ones=(bool *) malloc(num_prims*sizeof(bool))))
		fatal_error("Couldn't allocate memory for neg_ones.\n");
	if(!(zs=(int_double *) _aligned_malloc(num_prims*num_s*sizeof(int_double),16)))
		fatal_error("Couldn't allocate memory for zs.\n");

	for(q=q_start;q<=q_end;q++)
	{
		if((q&3)==2)
		  continue;
		fread(&q_file,sizeof(unsigned int),1,hur_file);
		if(q_file!=q)
		  {
		    printf("Mismatch between q=%d in file and internal counter=%d. Exiting.\n",q_file,q);
		    exit(1);
		  }
		fwrite(&q_file,sizeof(unsigned int),1,out_file);
		fread(&num_prims,sizeof(unsigned int),1,hur_file);
		fwrite(&num_prims,sizeof(unsigned int),1,out_file);

		for(prim=0;prim<num_prims;prim++)
		{
			fread(&indices[prim],sizeof(unsigned int),1,hur_file);
			fread(&omegas[prim],sizeof(int_complex),1,hur_file);
			fread(&neg_ones[prim],sizeof(bool),1,hur_file);
			fread(&zs[prim*num_s],sizeof(int_double),1,hur_file);
		}

		for(s=1;s<num_s;s++)
			for(prim=0;prim<num_prims;prim++)
				fread(&zs[prim*num_s+s],sizeof(int_double),1,hur_file);
/* Now out put the Z values
   in conjugate pairs
   */
		for(prim=0;prim<num_prims;prim++)
		  {
		    prim1=conj_j(prim,num_prims,q);
		    if(prim1<prim)
		      continue;
		    //printf("Saving q:%d ind:%d ",q,prim);
		    fwrite(&indices[prim],sizeof(unsigned int),1,out_file);
		    fwrite(&omegas[prim],sizeof(int_complex),1,out_file);
		    //if(indices[prim]==1699)print_int_complex_str("omega=",omegas[prim]);
		    fwrite(&neg_ones[prim],sizeof(bool),1,out_file);
		    fwrite(&zs[prim*num_s],sizeof(int_double),num_s,out_file);
		    if(prim1>prim) // it has a conjugate so save it
		      {
			//printf("ind:%d\n",prim1);
			fwrite(&indices[prim1],sizeof(unsigned int),1,out_file);
			//		    if(indices[prim1]==1699)print_int_complex_str("omega=",omegas[prim1]);
			fwrite(&omegas[prim1],sizeof(int_complex),1,out_file);
			fwrite(&neg_ones[prim1],sizeof(bool),1,out_file);
			fwrite(&zs[prim1*num_s],sizeof(int_double),num_s,out_file);
		      }
		    //else
		    //printf("\n");
		  }
/*
		for(prim=0;prim<num_prims;prim++)
		{
			fwrite(&indices[prim],sizeof(unsigned int),1,out_file);
			fwrite(&omegas[prim],sizeof(int_complex),1,out_file);
			fwrite(&neg_ones[prim],sizeof(bool),1,out_file);
			fwrite(&zs[prim*num_s],sizeof(int_double),num_s,out_file);
		}
*/
	}
	//printf("Tot Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
	return(0);
}


	



