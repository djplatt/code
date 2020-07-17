#ifndef HURWITZ
#define HURWITZ
#include "int_double8.0.h"

#define MAX_K (7)

int_double bernoulli[MAX_K];

void set_bernoulli()
{
	bernoulli[0]=int_double(1)/12;         // b(2)/2!
	bernoulli[1]=int_double(-1)/(30*24);   // b(4)/4!
	bernoulli[2]=int_double(1)/30240;      // b(6)/6!
	bernoulli[3]=int_double(-1)/1209600;    // b(8)/8!
	bernoulli[4]=int_double(1)/47900160;    // b(10)/10!
	bernoulli[5]=int_double(-691)/2730;
	bernoulli[5]/=479001600;                //b(12)/12!
	bernoulli[6]=int_double(7)/6;
	bernoulli[6]/=479001600;
	bernoulli[6]/=(13*14);                  //b(14)/14!

}

bool bernoulli_initialised=false;

#define DEFAULT_N (50)
int_complex hurwitz (int_complex s, int_double alpha)
{
	unsigned int i,k=MAX_K,n=DEFAULT_N;
	int_double err,n_alpha,s_mod=sqrt(norm(s));
	int_complex s1=s,res=c_zero,n_alpha_s,term;
	int_complex s_array[MAX_K];

	if(!bernoulli_initialised)
	{
		set_bernoulli();
		bernoulli_initialised=true;

	}

	if(n<-s_mod.right)
		n=(unsigned int) ceil(-s_mod.right);

	n_alpha=alpha+n;
	n_alpha_s=pow(n_alpha,-s+1);

	for(i=0;i<n;i++)
		res+=pow(alpha+i,-s);

	res+=n_alpha_s/(s-1);

	n_alpha_s/=n_alpha; // n_alpha_s=(n+alpha)^(-s)
		
	res+=n_alpha_s/2;

	s_array[0]=s*n_alpha_s/n_alpha;

	for(i=1;i<k;i++)
	{
		s1=s1+1;
		s_array[i]=s_array[i-1]*s1/n_alpha;
		s1=s1+1;
		s_array[i]*=s1/n_alpha;
	}

	for(i=0;i<k;i++)
	{
		term=s_array[i]*bernoulli[i];
		res+=term;
	}

	err=sqrt(norm(term))/(s.real+2*k-1)*(s_mod+2*k-1);
	if(err.left<=err.right)
		err.right=err.left;
	else
		err.left=err.right;

	res.real+=err;
	res.imag+=err;
	return(res);
}

#endif




