#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <float.h>
#include <complex>
//#define LINUX
#ifndef LINUX
#include <windows.h>
#endif
#include "../includes/int_double8.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft.h"



#define VERSION "1.0"
// 8/3/09
// 1.0 Original



void print_usage()
/* called when wrong arguments passed via command line */
{
	printf("Usage: int-fft-new q_start q_end im_start num_s gap facs_file out_file\n");
	exit(0);
}

void fatal_error(const char *error_string)
/* print the error message to stdout and exit */
{
	std::cout << error_string << endl;
	exit(0);
}



void make_l(unsigned int q, 
			unsigned int num_s, 
			factor *factors,
			unsigned int *q_n,
			im_s *im_s_vec,
			unsigned int *a_n,
			unsigned int *offset_vec,
			int_complex *omegas,
			FILE *out_file)
{
	unsigned int phi_q=factors[q].phi,fac,coords[MAX_FACS],conv_sizes[MAX_DIMS];
	int_complex omega,omega_a,z,z1,q_minus_s_2;
	int dims[MAX_FACS];
	int no_dims;
	bool power_2,primitive,neg_one;
	unsigned int pr=factors[q].pr,n_prims;
	unsigned int i,j,offset,s_done;
	int_complex s;

	printf("Processing Q=%d\n",q);

	s.real=d_half;
	fwrite(&q,sizeof(unsigned int),1,out_file);
	n_prims=num_prims(q,factors); // no of prim characters mod q
	fwrite(&n_prims,sizeof(unsigned int),1,out_file);


	if(pr)  // q has a primitive root so nice and easy 1-d FFT
	{

		fill_q_ns(q_n,pr,phi_q,q);

		init_bluestein_fft(phi_q,conv_sizes,bs,b_star_conjs);
		prep_omegas(omegas,q,phi_q,q_n,a,bs,b_star_conjs,conv_sizes[0],b_spare);
		s_done=0;
		while(s_done<num_s)
		{
			s.imag=int_double(im_s_vec[s_done].im_s);
			q_minus_s_2=pow(q,-s/2);
			for(i=1;i<q;i++)
				if(co_prime(i,q))
					fft_vec[q_n[i]]=hurwitz(s,int_double(i)/q);

			bluestein_fft(1,fft_vec,phi_q,a,bs,b_star_conjs,conv_sizes[0],b_spare);
			
			// now fft_vec contains L(chi,s)*q^(s/2)
			for(i=1;i<phi_q;i++)  // i=0 corresponds to principal chi
				if(prim_p(i,factors[q].primes[0]))				
				{
					neg_one=(i&1);

					if(s_done==0) // first time through, so finish off omega values
					{
						omegas[i]=finish_omega(omegas[i],neg_one);
						fwrite(&i,sizeof(unsigned int),1,out_file);
						fwrite(&omegas[i],sizeof(int_complex),1,out_file);
						fwrite(&neg_one,sizeof(bool),1,out_file);
					}
					z=fft_vec[i];			
					if(neg_one) // chi(-1)=-1
						z=z*omegas[i]*im_s_vec[s_done].lambda_s_a*q_minus_s_2;
					else
						z=z*omegas[i]*im_s_vec[s_done].lambda_s*q_minus_s_2;

					if(!contains_zero(z.imag))
					{
						printf("Imaginary part of z does not contain zero.\n");
						exit(0);
					}

					fwrite(&z.real,sizeof(int_double),1,out_file);
				}
				s_done++;
		}
	}
	else // q doesn't have a primitive root
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
		prep_omegas_nd(offset_vec,q,omegas,no_dims,dims,phi_q);

		while(s_done<num_s)
		{
			s.imag=int_double(im_s_vec[s_done].im_s);
			q_minus_s_2=pow(q,-s/2);
			for(i=0;i<phi_q;i++)
				fft_vec[offset_vec[i]]=hurwitz(s,int_double(a_n[i])/q);

			do_nd_fft(fft_vec,no_dims,dims,phi_q);

			for(i=0;i<factors[q].num_facs;i++)
				coords[i]=0;

			for(i=0;i<phi_q;i++)
			{
				primitive=true;
				for(j=0;j<factors[q].num_facs;j++)
					if(coords[j]%factors[q].primes[j]==0)
					{
						primitive=false;
						break;
					}

					if(primitive)
					{
						z=fft_vec[i];
						neg_one=neg_one_p(coords,factors[q].num_facs);
						if(power_2&&(i<(phi_q>>1)))
							neg_one=!neg_one;

						if(s_done==0) // first time for this im_s and chi
						{
							omegas[i]=finish_omega(omegas[i],neg_one);
							fwrite(&i,sizeof(unsigned int),1,out_file);
							fwrite(&omegas[i],sizeof(int_complex),1,out_file);
							fwrite(&neg_one,sizeof(bool),1,out_file);
						}

						if(neg_one)
							z*=omegas[i]*im_s_vec[s_done].lambda_s_a*q_minus_s_2;
						else
							z*=omegas[i]*im_s_vec[s_done].lambda_s*q_minus_s_2;

						if(!contains_zero(z.imag))
						{
							printf("Non Zero Imag Part %d %d %d\n",q,i,j);
							exit(0);
						}

						fwrite(&z.real,sizeof(int_double),1,out_file);
					}

						j=factors[q].num_facs-1;
						coords[j]++;
						while(coords[j]==factors[factors[q].facs[j]].phi)
						{
							coords[j]=0;
							if(j==0)
								break;
							j--;
							coords[j]++;
					}
			}
			s_done++;
		}
	}


}


int main(int argc, char **argv)
{
	unsigned int q_start,q_end,num_s,q,i,*q_n,*offset_vec,*a_n;
	FILE *out_file;
	ifstream facs_file;
	factor *factors;
	int_complex *omegas,s,pi_minus_s_2;
	im_s *im_s_vec;
	double im_start,gap,im;


	clock_t no_clicks;

	no_clicks=clock(); // start timing

	_fpu_rndd;




	if(argc!=8)
		print_usage();
	q_start=atoi(argv[1]);
	if(q_start<3)
		fatal_error("q_start must be >=3.");
	q_end=atoi(argv[2]);
	if(q_end<q_start)
		fatal_error("q_end must be >= q_start.");
	im_start=atof(argv[3]);
	if(im_start<0.0)
		fatal_error("im_start must be >= 0.0.");
	num_s=atoi(argv[4]);
	if(num_s<=0)
		fatal_error("num_s must be > 0.");
	gap=atof(argv[5]);
	if(gap<=0.0)
		fatal_error("gap must be > 0.0.");
	facs_file.open(argv[6]);
	if(!facs_file.is_open())
		fatal_error("Couldn't open factors file. Exiting.\n");
	out_file=fopen(argv[7],"wb");
	if(!out_file)
		fatal_error("Couldn't open output file. Exiting.\n");
	if(q_end>MAX_Q)
	{
		printf("Fatal error, Q=%d exceeds MAX_Q=%d. Change defines and recompile.\n",q_end,MAX_Q);
		exit(0);
	}

	printf("Running from q= %d to %d.\n",q_start,q_end);

	if(!(im_s_vec=(im_s *) calloc(num_s,sizeof(im_s))))
		fatal_error("Error allocating memory for im_s_vec. Exiting.\n");

	s.real=d_half;
	for(i=0,im=im_start;i<num_s;i++,im+=gap)
	{
		s.imag=int_double(im);
		pi_minus_s_2=pow(d_pi,-s/2);
		im_s_vec[i].im_s=im;
		sin_cos(imag(lngamma(s/2)),&im_s_vec[i].lambda_s.imag,&im_s_vec[i].lambda_s.real);
		im_s_vec[i].lambda_s*=pi_minus_s_2;
		sin_cos(imag(lngamma((s+1)/2)),&im_s_vec[i].lambda_s_a.imag,&im_s_vec[i].lambda_s_a.real);
		im_s_vec[i].lambda_s_a*=pi_minus_s_2;
	}

	fwrite(&num_s,sizeof(unsigned int),1,out_file);
	fwrite(im_s_vec,sizeof(im_s),num_s,out_file);

	if(!(factors=(factor *) calloc(q_end+1,sizeof(factor))))
		fatal_error("Error allocating memory for factors. Exiting.\n");

	if(!read_factors(factors,q_end,facs_file))
		fatal_error("Error reading factor file. Exiting.\n");

	facs_file.close();

	if(!(omegas=(int_complex *) calloc(q_end-1,sizeof(int_complex))))
		fatal_error("Error allocating memory for omegas. Exiting.\n");

	if(!(q_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
		fatal_error("Error allocating memory for q_n. Exiting.\n");

	if(!(a_n=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
		fatal_error("Error allocating memory for a_n. Exiting.\n");

	if(!(offset_vec=(unsigned int *) calloc(q_end-1,sizeof(unsigned int))))
		fatal_error("Error allocating memory for offset_vec. Exiting.\n");

#ifndef LINUX
	manage_memory(factors,omegas,q_n,a_n,offset_vec);
#endif
	init_ws(_ws);

	for(q=q_start;q<=q_end;q++)
	{
		// hur_file is pointing to H(s_0,1/q_start)
		if((q&3)!=2)
			make_l(q,num_s,factors,q_n,
			im_s_vec,offset_vec,a_n,omegas,out_file);
	}


	fclose(out_file);

	printf("Time Elapsed = %8.2f secs.\n",((double) clock()-no_clicks)/((double) CLOCKS_PER_SEC));
	printf("%ld length 2 FFTs computed.\n",fft_2_count);
	printf("%ld length 4 FFTs computed.\n",fft_4_count);
	printf("%ld power of 2 length FFTs computed.\n",fft_power_2_count);
	return(0);
}


