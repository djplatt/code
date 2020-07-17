#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <complex>
#include "../includes/int_double10.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft.h"
#define POS 0
#define NEG 1
#define UNK 2

factor factors[MAX_Q+1];
unsigned int N_POINTS;

inline char sign(const int_double &x)
{
	if(x.right>0)
		return(NEG);
	if(x.left>0)
		return(POS);
	return(UNK);
}

// calculate a z value
// z is sum chi(n)zeta(s,a/q) for a=1..q-1
int_double z_val(unsigned int q, int_complex &z, int_complex &s, int_complex &omega, bool neg_one)
{
	int_complex res;
	double delta;
	if(neg_one)
		delta=1.0;
	else
		delta=0.0;
	sin_cos(imag(lngamma((s+delta)/2)),&res.imag,&res.real);
	res*=pow(d_pi,-(s+delta)/2);
	res*=omega*z*pow(q,-s/2);

	if(!contains_zero(res.imag))
		printf("Error, Imaginary part of res does not contain zero.\n");
	return(res.real);
}


unsigned int q_n[MAX_Q],a_n[MAX_Q],offset_vec[MAX_Q];

bool find_zero(unsigned int q, unsigned int index, double &im_s0, double &im_s1, int_complex &omega, bool neg_one,
			   char this_sign,double *im_s)
{
	unsigned int phi_q=factors[q].phi,fac;
	int dims[MAX_FACS];
	unsigned int no_dims;
	bool power_2;
	unsigned int pr=factors[q].pr;
	unsigned int i,j,offset,s_done;
	int_complex s;
	int_double zv;
	char next_sign;

	printf("q: %ld index: %ld\n",q,index);

	for (i=0;i<N_POINTS-1;i++)
		im_s[i]=im_s0+(i+1)*(im_s1-im_s0)/N_POINTS;

	if(pr)
	{
		fill_q_ns(q_n,pr,phi_q,q);
		init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
		for(i=0;i<N_POINTS-1;i++)
		{
			s=int_complex(int_double(0.5),int_double(im_s[i]));
//			s=int_complex(int_double(0.5),int_double(im_s0));
			for(j=1;j<q;j++)
				if(co_prime(j,q))
					fft_vec[q_n[j]]=hurwitz(s,int_double(j)/q);
			bluestein_fft(1,fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare);
			next_sign=sign(z_val(q,fft_vec[index],s,omega,neg_one));
			if((next_sign!=UNK)&&(next_sign!=this_sign))
				break;
		}
		if(i==N_POINTS-1)
		{
			printf("No change of sign found.\n");
			return(false);
		}
	}
	else
	{
		j=0;
		for(i=1;i<q;i++)
			if(co_prime(i,q))
				a_n[j++]=i;
		no_dims=factors[q].num_facs;
		fac=factors[q].facs[0];        // the first p^n
		power_2=(factors[fac].pr==0);  // no conductor => p=2, n>=3
		if(power_2)
		{
			no_dims++;
			fill_q_ns_2s(q_n,fac);    // use the {-1,1}X{5} trick
			for(i=1;i<factors[q].num_facs;i++) // move all the dimensions up one
				dims[i+1]=factors[factors[q].facs[i]].phi;
			dims[1]=factors[q].facs[0]/4;
			dims[0]=2;                         // slip in a two
		}
		else
		{
			fill_q_ns(q_n,factors[fac].pr,factors[fac].phi,fac); // use the generator
			for(i=0;i<factors[q].num_facs;i++)
				dims[i]=factors[factors[q].facs[i]].phi;
		}

		offset=fac;
		for(i=1;i<factors[q].num_facs;i++)  // do the rest on the factors, all will have generators
		{
			fac=factors[q].facs[i];			
			pr=factors[fac].pr;
			fill_q_ns(&q_n[offset],pr,factors[fac].phi,fac);  // use the generator
			offset+=fac;
		}
		s_done=0;
		make_offsets(q,q_n,offset_vec,factors,a_n);    // reverse q_n so we know how to populate fft vector

		unsigned int w_ptr=0;
		for(i=0;i<no_dims;i++)
		{
			ws_ptr[i]=w_ptr;
			if((dims[i]==2)||dims[i]==4)
				continue;

			if(bin_power_2(dims[i]))
			{
				w_ptr+=dims[i];
				continue;
			}
			init_bluestein_fft(dims[i],&conv_sizes[i],&bs[w_ptr],&b_star_conjs[w_ptr]);
			w_ptr+=conv_sizes[i];
		}

		
		for(i=0;i<N_POINTS-1;i++)
		{
			s=int_complex(int_double(0.5),int_double(im_s[i]));
			offset=0;
			for(j=1;j<q;j++)
				if(co_prime(j,q))
					fft_vec[offset_vec[offset++]]=hurwitz(s,int_double(j)/q);
			do_nd_fft(fft_vec,no_dims,dims,phi_q);
			next_sign=sign(z_val(q,fft_vec[index],s,omega,neg_one));
			if((next_sign!=UNK)&&(next_sign!=this_sign))
				break;
		}
		if(i==N_POINTS-1)
		{
			printf("No change of sign found.\n");
			return(false);
		}

	}

	return(true);
}

// reads stationary point records from argv[1]
// each record consists of
// q: unsigned int
// index: unsigned int 
// s0: double - the false stationary point is between
// s1: double -     1/2+i*s0 and 1/2+i*s1
// omega: int_complex - the root number for this q,index
// neg_one: true if chi(-1)=-1
// this_sign: the sign of Z(0.5+i*s0)
//
// argv[2] is the facs_file
//
// argv[3] N_POINTS is an integer>=3. The gap between s0 and s1 is
// divided into this many gaps. Make successive ones coprime
// to avoid duplicating points.
//
// argv[4] an output file for records in same format as input
//         file for stationary points where we fail to find a
//         zero at this resolution.
int main(int argc, char **argv)
{
	FILE *infile,*outfile;
	ifstream facs_file;
	unsigned int q,index;
	double im_s0,im_s1,*im_s;
	char this_sign;
	int_complex omega;
	
	bool neg_one;

	_fpu_rndd();

	if(argc!=5)
	{
		printf("Incorrect command line. Exiting.\n");
		exit(0);
	}

	infile=fopen(argv[1],"rb");
	if(!infile)
	{
		printf("Failed to open %s for binary input. Exiting.\n",argv[1]);
		exit(0);
	}

	facs_file.open(argv[2]);
	if(!facs_file)
	{
		printf("Failed to open %s for input. Exiting.\n",argv[2]);
		exit(0);
	}

	if(!read_factors(factors,MAX_Q,facs_file))
	{
		printf("Fatal error reading factors file. exiting.\n");
		exit(0);
	}
	facs_file.close();

	N_POINTS=atoi(argv[3]);
	if(N_POINTS<3)
	{
		printf("N_POINTS must be at least 3. Exiting.\n");
		exit(0);
	}

	im_s=(double *) malloc(sizeof(double)*(N_POINTS-1));
	if(!im_s)
	{
		printf("Error allocating memory for im_s. Exiting.\n");
		exit(0);
	}

	outfile=fopen(argv[4],"wb");
	if(!outfile)
	{
		printf("Failed to open file %s for binary output. Exiting\n");
		exit(0);
	}
	init_ws(_ws);

	while(fread(&q,sizeof(unsigned int),1,infile))
	{
		fread(&index,sizeof(unsigned int),1,infile);
		fread(&im_s0,sizeof(double),1,infile);
		fread(&im_s1,sizeof(double),1,infile);
		fread(&omega,sizeof(int_complex),1,infile);
		fread(&neg_one,sizeof(bool),1,infile);
		fread(&this_sign,sizeof(char),1,infile);
		if(!find_zero(q,index,im_s0,im_s1,omega,neg_one,this_sign,im_s))
		{
			fwrite(&q,sizeof(unsigned int),1,outfile);
			fwrite(&index,sizeof(unsigned int),1,outfile);
			fwrite(&im_s0,sizeof(double),1,outfile);
			fwrite(&im_s1,sizeof(double),1,outfile);
			fwrite(&omega,sizeof(int_complex),1,outfile);
			fwrite(&neg_one,sizeof(bool),1,outfile);
			fwrite(&this_sign,sizeof(char),1,outfile);
			printf("q: %ld index: %ld Failed to find change of sign between Im(s)= %10.8e and %10.8e.\n",
				q,index,im_s0,im_s1);
		}
	}

	fclose(infile);
	fclose(outfile);
	return(0);
}
