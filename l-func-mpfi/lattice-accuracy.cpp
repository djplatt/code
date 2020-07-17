#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <complex>
#include "../includes/int_double12.0.h"

// created: 26/01/09
// last edit: 26/01/09

#define ERR (-14)

inline int gap (double x, double y)
{
  __int128_t *x1,*y1,res;
  x1=(__int128_t *) &x;
  y1=(__int128_t *) &y;
  res=x1[0]-y1[0];
  if(res<0)
    return(-res);
  else
    return(res);
}


int main(int argc, char **argv)
{
	FILE* infile;
	int num_s,n_out,rn_out,no_gaps,i,j,k,rel_err,errs=0;
	double im_s,z[4];
	int max_gap=0,new_gap,worst_k,worst_j,worst_ak,worst_aj;
	double worst_left,worst_right,worst_aleft,worst_aright,max_abs_gap=0.0,new_abs_gap;
	_fpu_rndd();
	//im_s=1.0;
	//printf("nextafter(1.0)=%A\n",nextafter(im_s));
	//printf("im_s          =%A\n",im_s);
	//exit(0);
	if(argc!=2)
	{
		printf("Usage:- lattice-accuracy <infile>.\n");
		exit(0);
	}
	if(!(infile=fopen(argv[1],"rb")))
	{
		printf("Failed to open %s for binary input. Exiting.\n",argv[1]);
		exit(0);
	}

	fread(&num_s,sizeof(int),1,infile);
	fread(&n_out,sizeof(int),1,infile);

	fread(&rn_out,sizeof(int),1,infile);
	fread(&no_gaps,sizeof(int),1,infile);

	for(i=0;i<num_s;i++)
	{
		fread(&im_s,sizeof(double),1,infile);
		fread(z,sizeof(double),4,infile);
		fread(z,sizeof(double),4,infile);
		fread(z,sizeof(double),4,infile);
		for(j=0;j<no_gaps+1;j++)
		{
			for(k=0;k<n_out;k++)
			{
			  fread(z,sizeof(double),4,infile);
			  //			  z[1]=-z[1];z[3]=-z[3];
			  
			  new_abs_gap=z[1]-z[0];//gap(z[0],z[1]);
			  if(new_abs_gap>max_abs_gap)
			    {
			      max_abs_gap=new_abs_gap;
			      worst_aleft=z[0];
			      worst_aright=z[1];
			      worst_aj=j;
			      worst_ak=k;
			    }
			  new_abs_gap=z[3]-z[2];//gap(z[2],z[3]);
			  if(new_abs_gap>max_abs_gap)
			    {
			      max_abs_gap=new_abs_gap;
			      worst_aleft=z[2];
			      worst_aright=z[3];
			      worst_aj=j;
			      worst_ak=k;
			    }
			  
			  new_gap=gap(z[0],z[1]);
			  if(new_gap>max_gap)
			    {
			      max_gap=new_gap;
			      worst_left=z[0];
			      worst_right=z[1];
			      worst_j=j;
			      worst_k=k;
			    }
			  new_gap=gap(z[2],z[3]);
			  if(new_gap>max_gap)
			    {
			      max_gap=new_gap;
			      worst_left=z[2];
			      worst_right=z[3];
			      worst_j=j;
			      worst_k=k;
			    }
			}
		}
	}

	//printf("Worst case gap was %30.28e [%30.28e,%30.28e].\n",max_abs_gap,worst_left,worst_right);
	printf("Worst case rel gap was %d [%30.28e,%30.28e] absolute error %30.28e.\n",max_gap,worst_left,worst_right,worst_right-worst_left);
	printf("Happened at j=%d k=%d\n",worst_j,worst_k);
	printf("Worst case abs gap was %30.28e [%30.28e,%30.28e].\n",max_abs_gap,worst_aleft,worst_aright);
	printf("Happened at j=%d k=%d\n",worst_aj,worst_ak);
	return(0);
}



