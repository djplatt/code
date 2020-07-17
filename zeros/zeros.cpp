#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <complex>
#include "../includes/int_double11.0.h"
#include "../includes/im_s.h"


#define MAX_Q (100000)
#define POS 0
#define NEG 1
#define UNK 2
#define UP 0
#define DOWN 1

using namespace std;

char sign(int_double &x)
{
	if(x.left>0.0)
		return(POS);
	if(x.right>0.0)
		return(NEG);
	return(UNK);
}

void print_sign(char sign)
{
	if(sign==POS)
		printf("+");
	else
	{
		if(sign==NEG)
			printf("-");
		else
			printf("?");
	}
}

void print_signs(unsigned int n, int_double *re_zs, unsigned int step, im_s *im_s_vec)
{
	unsigned int i;

	printf("%6.5e\n",im_s_vec[0].im_s);
	print_sign(sign(re_zs[0]));

	for(i=step;i<n;i+=step)
	{
		if(sign(re_zs[i])!=sign(re_zs[i-step]))
		{
			printf("\n%6.5e\n",im_s_vec[i].im_s);
		}
		print_sign(sign(re_zs[i]));
	}
	printf("\n");
}

unsigned int do_stat_points(int_double *re_zs, unsigned int n, unsigned int q, unsigned int index, 
					unsigned int step, im_s *im_s_vec, FILE *out_file,int_complex &omega, bool neg_one,
					unsigned int MIN_N)
{
	unsigned int missed=0;
	unsigned int i;
	char this_sign;

	for(i=MIN_N;i+step+step<n;i+=step)
	{
		this_sign=sign(re_zs[i]);
		if(this_sign==UNK)
			continue;
		if(this_sign==POS)
		{
			if(sign(re_zs[i+step])==POS)
				if(sign(re_zs[i+step+step])==POS)
					if((re_zs[i]>re_zs[i+step])&&(re_zs[i+step]<re_zs[i+step+step]))
					{
						missed+=2;
//						printf("q: %d Index: %d ",q,index);
//						printf("Min found above zero in Im(s)=[%6.5e,%6.5e]\n",im_s_vec[i].im_s,im_s_vec[i+2*step].im_s);
						if(out_file)
						{
							fwrite(&q,sizeof(unsigned int),1,out_file);
							fwrite(&index,sizeof(unsigned int),1,out_file);
							fwrite(&im_s_vec[i].im_s,sizeof(double),1,out_file);
							fwrite(&im_s_vec[i+2*step].im_s,sizeof(double),1,out_file);
							fwrite(&omega,sizeof(int_complex),1,out_file);
							fwrite(&neg_one,sizeof(bool),1,out_file);
							fwrite(&this_sign,sizeof(char),1,out_file);
						}
						else
							printf("Stat point in Turing zone, Q: %ld index: %ld\n",q,index);
					}
		}
		else
		{
			if(sign(re_zs[i+step])==NEG)
				if(sign(re_zs[i+step+step])==NEG)
					if((re_zs[i]<re_zs[i+step])&&(re_zs[i+step]>re_zs[i+step+step]))
					{
						missed+=2;
//						printf("q: %d Index: %d ",q,index);
//						printf("Max found below zero in Im(s)=[%6.5e,%6.5e]\n",im_s_vec[i].im_s,im_s_vec[i+2*step].im_s);
						if(out_file)
						{
							fwrite(&q,sizeof(unsigned int),1,out_file);
							fwrite(&index,sizeof(unsigned int),1,out_file);
							fwrite(&im_s_vec[i].im_s,sizeof(double),1,out_file);
							fwrite(&im_s_vec[i+2*step].im_s,sizeof(double),1,out_file);
							fwrite(&omega,sizeof(int_complex),1,out_file);
							fwrite(&neg_one,sizeof(bool),1,out_file);
							fwrite(&this_sign,sizeof(char),1,out_file);
						}
						else
							printf("Stat point in Turing zone, Q: %ld index: %ld\n",q,index);

					}
		}
	}

	return(missed);
}

void read_check(unsigned int expect, unsigned int i, FILE **infiles)
{
	unsigned int tmp;
	if(fread(&tmp,sizeof(unsigned int),1,infiles[i]))
		if(tmp==expect)
			return;
	printf("Data mismatch between input files. Exiting.\n");
	exit(0);
}

#define little_step(a) (sizeof(int_double))
#define big_step(a) (little_step(a)+sizeof(int_complex)+sizeof(bool)+sizeof(unsigned int))

// integrate observed N(t) from a..b where N(a)=0
int_double n_twiddle(unsigned int a, unsigned int b, int_double *re_zs, im_s *im_s_vec)
{
	unsigned int ptr=a+1,zeros_so_far=0;
	char last_sign=sign(re_zs[a]),this_sign;
	int_double res=int_double(0.0),last_change=int_double(im_s_vec[a].im_s),this_change;
	int_double end_point=int_double(im_s_vec[b].im_s);

	// skip leading UNK's
	while(last_sign==UNK)
	{
		if(ptr>b)
			return(res);
		last_sign=sign(re_zs[ptr++]);
	}

	while(ptr<=b)
	{
		for(this_sign=sign(re_zs[ptr++]);this_sign==last_sign;this_sign=sign(re_zs[ptr++]))
			if(ptr>b)
				break;
		if(this_sign!=UNK)
		{
			res+=end_point-int_double(im_s_vec[ptr-2].im_s,im_s_vec[ptr-1].im_s);
			last_sign=this_sign;
		}
		else
		{
			if(ptr>b) // the last entry was an UNK
				return(res);
			this_sign=sign(re_zs[ptr++]);
			if(this_sign==UNK)
			{
				printf("Fatal error, two Uknowns in succession. Exiting.\n");
				exit(0);
			}
			if(this_sign!=last_sign) // isolated unknowns get ignored e.g. ++++?++++
			{
				last_sign=this_sign;
				res+=end_point-int_double(im_s_vec[ptr-3].im_s,im_s_vec[ptr-1].im_s);
			}
		}
	}
	return(res);
}

int_double s_of_t(unsigned int q, double t2) // Rumely Th. 2 p 429
{
	int_double tmp2=int_double (t2),tmp=int_double(1242), res=int_double(18397);
	tmp/=10000;res/=10000;
	tmp2*=q;
	tmp2/=2;
	tmp2/=d_pi;
	res=res+tmp*log(tmp2);
	res.left=res.right;
	return(res);
}

// returns (t1+t2)*ln(q/pi)/4
int_double ln_term (double t1, double t2, unsigned int q)
{
	int_double res=int_double(q);
	res/=d_pi;
	return(log(res)*(t1+t2)/4);
}

// count zeros between a and b
int num_zeros(int a, int b, int_double *re_zs)
{
	int zero_count=0,ptr=a+1;
	char this_sign,last_sign=sign(re_zs[a]);
	while(last_sign==UNK)
	{
		last_sign=sign(re_zs[ptr++]);
		if(ptr>b)
		{
			printf("Catastrophic failure in num_zeros. Exiting.\n");
			exit(0);
		}
	}

	while(ptr<=b)
	{
		this_sign=sign(re_zs[ptr++]);
		while((this_sign==UNK)&&(ptr<=b))
			this_sign=sign(re_zs[ptr++]);
		if(this_sign==UNK) // file finished with one or more unknowns
			return(zero_count);
		if(this_sign!=last_sign)
		{
			zero_count++;
			last_sign=this_sign;
		}
	}
	return(zero_count);
}	

// compare number of zeros observed with that calculate by Turing's method
// If Turing's method does not identify an integer, returns false
// If Turing's method identifies an integer, but there are more/less zeros
//    found, then return false
// Otherwise return true
int check_zeros(unsigned int t1_ptr, unsigned int t2_ptr,unsigned int q,
				 int_double &gamma_int, int_complex &omega, int_double *re_zs, im_s *im_s_vec,
				 unsigned int index)
{
	double t1=im_s_vec[t1_ptr].im_s;
	double t2=im_s_vec[t2_ptr].im_s;
	int_double ntwiddle=n_twiddle(t1_ptr,t2_ptr,re_zs,im_s_vec);

	int_double arg_omega;
	
	if(re_zs[0].right>=0.0) // negative
		omega=-omega;
	arg_omega=argument(omega);

	int_double s_t=s_of_t(q,t2);
	int_double lnt=ln_term(t1,t2,q);
	int_double calc_zeros;
	int calc_zeros_int;
	int count_zeros=num_zeros(0,t1_ptr,re_zs); // number of zeros found from t1 to t2

// temporary stuff
	int fft_ptr=0;
	if(index==821)
	{
		for(fft_ptr=0;fft_ptr+512<=t1_ptr;fft_ptr+=512)
			printf("Zeros %ld-%ld = %ld\n",fft_ptr,fft_ptr+511,num_zeros(fft_ptr,fft_ptr+511,re_zs));
		calc_zeros=(arg_omega+lnt+gamma_int/(t2-t1))/d_pi-ntwiddle/(t2-t1)+s_t/(t2-t1); // Turing's estimate
		print_int_double_str("Calculated # zeros = ",calc_zeros);
		FILE *temp;
		temp=fopen("z_821.dat","wb");
		fwrite(re_zs,sizeof(int_double),t1_ptr,temp);
		fclose(temp);


//		printf("Zeros %ld-%ld = %ld\n",fft_ptr,t1_ptr,num_zeros(fft_ptr,t1_ptr,re_zs));
//		printf("Zeros %ld-%ld = %ld\n",t1_ptr-511,t1_ptr,num_zeros(t1_ptr-511,t1_ptr,re_zs));
		exit(0);
	}

	// end of temporary stuff



	calc_zeros=(arg_omega+lnt+gamma_int/(t2-t1))/d_pi-ntwiddle/(t2-t1)+s_t/(t2-t1); // Turing's estimate
	
	calc_zeros_int=(int)ceil(calc_zeros.left);

	if(calc_zeros_int!=floor(-calc_zeros.right)) // Turing's estimate did not bracket a unique integer
	{
		printf("q: %ld index: %ld Problem with Turing's estimate:- ",q,index);
		print_int_double(calc_zeros);
		printf("\n");
		return(-1);
	}

//	printf("q: %ld Index: %ld Turing's estimate was %ld but I found %ld.\n",q,index,calc_zeros_int,count_zeros);
	return(calc_zeros_int-count_zeros);
}

unsigned int conj_j (unsigned int j, unsigned int num_chi, unsigned int q)
{
	if(q&1) // odd q
		return(num_chi-j-1);
	if(j<(num_chi>>1))
		return((num_chi>>1)-j-1);
	return(num_chi+(num_chi>>1)-j-1);
}

#define MAX_F_NAME (256)

int main(int argc, char **argv)
{
	FILE *spec_file,**infiles,*out_file;
	im_s *im_s_vec;
	unsigned int num_files,num_top_files,i,j,num_s,index,*num_ss,*start_pos;
	unsigned int q,num_chi,re_z_ptr;
	unsigned int missed=0,first_top=0,last_top;
	int_double gamma,gamma_a;
	int_double *re_zs;
	char fname[MAX_F_NAME];
	bool neg_one;
	int_complex omega;
	int diffs[MAX_Q];

	_fpu_rndd();

	if(argc!=3)
	{
		printf("Incorrect command line. Exiting.\n");
		exit(0);
	}

	spec_file=fopen(argv[1],"r");
	if(!spec_file)
	{
		printf("Failed to open %s for input. Exiting.\n",argv[1]);
		exit(0);
	}

	out_file=fopen(argv[2],"wb");
	if(!out_file)
	{
		printf("Failed to open %s for binary output. Exiting.\n",argv[2]);
		exit(0);
	}

	fscanf(spec_file,"%ld %ld\n",&num_files,&num_top_files);
	infiles=(FILE **) malloc(sizeof(spec_file)*num_files);
	if(!infiles)
	{
		printf("Fatal error allocating memory for infiles. Exiting.\n");
		exit(0);
	}

	num_ss=(unsigned int *)malloc(sizeof(unsigned int)*num_files);
	if(!num_ss)
	{
		printf("Fatal error allocating memory for num_ss. Exiting.\n");
		exit(0);
	}

	start_pos=(unsigned int *)malloc(sizeof(unsigned int)*num_files);
	if(!start_pos)
	{
		printf("Fatal error allocating memory for start_pos. Exiting.\n");
		exit(0);
	}

	for(i=0;i<num_files;i++)
	{
		fscanf(spec_file,"%s\n",fname);
		infiles[i]=fopen(fname,"rb");
		if(!infiles[i])
		{
			printf("Failed to open file %s for binary input. Exiting.\n",fname);
			exit(0);
		}
	}

	fscanf(spec_file,"%ld",&i);
	gamma.left=i;
	gamma.right=-(int)i-1;
	fscanf(spec_file,"%ld",&i);
	gamma_a.left=i;
	gamma_a.right=-(int)i-1;
	fscanf(spec_file,"%ld",&i);
	gamma/=i;
	gamma_a/=i;
	fclose(spec_file);

	num_s=0;
	for(i=0;i<num_files;i++)
	{
		fread(&num_ss[i],sizeof(unsigned int),1,infiles[i]);
		num_s+=num_ss[i];
	}

	for(i=0;i<num_files-num_top_files;i++)
		first_top+=num_ss[i];
	last_top=first_top-1;
	for(;i<num_files;i++)
		last_top+=num_ss[i];

	im_s_vec=(im_s *) _aligned_malloc(sizeof(im_s)*num_s,16);
	if(!im_s_vec)
	{
		printf("Fatal error allocating memory for im_s_vec. Exting.\n");
		exit(0);
	}

	re_zs=(int_double *) _aligned_malloc(sizeof(int_double)*num_s,16);
	if(!re_zs)
	{
		printf("Fatal error allocating memory for re_zs. Exting.\n");
		exit(0);
	}

	index=0;
	for(i=0;i<num_files;i++)
	{
		fread(&im_s_vec[index],sizeof(im_s),num_ss[i],infiles[i]);
		index+=num_ss[i];
	}

	while(fread(&q,sizeof(unsigned int),1,infiles[0]))
	{
		if(q>MAX_Q)
		{
			printf("q exceeds MAX_Q. Exiting.\n");
			exit(0);
		}
		fread(&num_chi,sizeof(unsigned int),1,infiles[0]);
		for(i=1;i<num_files;i++)
		{
			read_check(q,i,infiles);
			read_check(num_chi,i,infiles);
		}

		for(j=0;j<num_chi;j++)
		{
			re_z_ptr=0;
			for(i=0;i<num_files;i++)
			{
				fread(&index,sizeof(unsigned int),1,infiles[i]);
				fread(&omega,sizeof(int_complex),1,infiles[i]);
				fread(&neg_one,sizeof(bool),1,infiles[i]);
				fread(&re_zs[re_z_ptr],sizeof(int_double),num_ss[i],infiles[i]);
				re_z_ptr+=num_ss[i];
			}



// look at the Turing zone to see if there are any obvious missing zeros
			do_stat_points(&re_zs[first_top],last_top-first_top,q,j,1,&im_s_vec[first_top],NULL,omega,neg_one,0);

			diffs[j]=check_zeros(first_top,last_top,q,(neg_one ? gamma_a : gamma),omega,re_zs,im_s_vec,index);
			diffs[j]-=do_stat_points(re_zs,first_top,q,index,1,im_s_vec,out_file,omega,neg_one,10); // assume we'll find
			                                                                                        // these zeros
		}

		for(j=0;j<num_chi;j++)
		{
			if(diffs[j]!=0)
				if(diffs[j]+diffs[conj_j(j,num_chi,q)]!=0)
					printf("Problem at q: %ld j: %ld diffs[%ld]: %ld diffs[%ld]: %ld.\n",q,j,j,diffs[j],num_chi-j-1,diffs[num_chi-j-1]);
		}
	}
}
