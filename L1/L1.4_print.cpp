/*
File: L1.4_print.cpp

Created: 30th Sept 2013

Implementation notes:

Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

By: DJ Platt
    Bristol University

Copyright 2011.

This work is funded by the UK ESPRC. */

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"
#define LINUX
#include "../includes/int-fft-half.h"

void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: L1.4_print (ifname)\n");
  printf("  (ifname)    - file with output from L1.4.\n");
  exit(1);
}

void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};

inline int phi(int p, int p_n)
{
  return(p_n-p_n/p);
}

long unsigned int pow_mod(long unsigned int a, long unsigned int b, long unsigned int m)
{
  long unsigned int a_pw=a,pw=b,res=1;
  while(true)
    {
      //printf("pw=%ld a_pw=%ld res=%ld\n",pw,a_pw,res);
      if(pw&1)
	res=(res*a_pw)%m;
      pw>>=1;
      if(pw==0)
	return(res);
      a_pw=(a_pw*a_pw)%m;
    }
}
      
unsigned long int pr(unsigned long int i, factor *factors)
{
  int phi=factors[i].phi;
  for(int p=2;p<i;p++)
    {
      if(gcd(p,i)!=1)
	continue;
      bool good=true;
      for(int j=0;j<factors[phi].num_facs;j++)
	{
	  if(pow_mod(p,phi/factors[phi].primes[j],i)!=1)
	    continue;
	  good=false;
	  break;
	}
      if(good)
	return(p);
    }
}

bool make_factors(factor *factors, int q_end)
{
  return(true);
  printf("Making factor database.\n");
  //printf("pow_mod(7,15,99)=%ld\n",pow_mod(7,19,99));
  int *primes,max_p=floor(sqrt(q_end));
  primes=(int *) malloc(sizeof(int)*(q_end+1));
  for(int i=2;i<=q_end;i++)
    primes[i]=0;
  for(int i=4;i<=q_end;i+=2)
    primes[i]=2;
  for(int i=3;i<=max_p;i+=2)
    if(primes[i]==0)
      for(int j=i*i;j<=q_end;j+=i)
	if(primes[j]==0)
	  primes[j]=i;
  printf("Prime sieve completed.\n");
  // now each entry primes[i] is 0 if i prime, else = smallest prime factor
  factors[3].num_facs=1;
  factors[3].primes[0]=3;
  factors[3].facs[0]=3;
  factors[4].num_facs=1;
  factors[4].primes[0]=2;
  factors[4].facs[0]=4;

  for(int f=5;f<=q_end;f++)
    if(primes[f]==0) // a prime
      {
	factors[f].num_facs=1;
	factors[f].primes[0]=f;
	factors[f].facs[0]=f;
      }
    else
      {
	factors[f].primes[0]=primes[f];
	int n=f/primes[f];
	if(factors[n].primes[0]==primes[f])
	  {
	    factors[f].num_facs=factors[n].num_facs;
	    factors[f].facs[0]=primes[f]*factors[n].facs[0];
	    for(int i=1;i<factors[n].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i];
		factors[f].facs[i]=factors[n].facs[i];
	      }
	  }
	else
	  {
	    factors[f].num_facs=factors[n].num_facs+1;
	    factors[f].facs[0]=primes[f];
	    factors[f].primes[0]=primes[f];
	    for(int i=1;i<factors[f].num_facs;i++)
	      {
		factors[f].primes[i]=factors[n].primes[i-1];
		factors[f].facs[i]=factors[n].facs[i-1];
	      }
	  }
      }
  free(primes);
  printf("Factors computed.\n");
  // now calculate phi(f)
  for(int i=3;i<=q_end;i++)
    {
      factors[i].phi=1;
      for(int j=0;j<factors[i].num_facs;j++)
	factors[i].phi*=phi(factors[i].primes[j],factors[i].facs[j]);
    }
  printf("phi computed.\n");

  //now do the prim roots
  factors[3].pr=2;
  factors[4].pr=3;
  //long unsigned int max_pr=3;
  for(int i=5;i<=q_end;i++)
    {
      if(((factors[i].num_facs==1)&&(factors[i].primes[0]!=2))|| // p^m, p an odd prime
	 ((factors[i].num_facs==2)&&(factors[i].facs[0]==2)))    // 2p^m
	{
	  factors[i].pr=pr(i,factors);
	  /*
	  if(factors[i].pr>max_pr)
	    {
	      max_pr=factors[i].pr;
	      printf("New largest pr=%lu for modulus %ld\n",max_pr,i);
	    }
	  */
	}
      else
	factors[i].pr=0;
    }
  printf("pr's computed.\n");
  /*  
  for(int i=3;i<50;i++)
    {
      printf("%ld has %ld factors, phi(%ld)=%ld, pr(%ld)=%ld, factors are",i,factors[i].num_facs,i,factors[i].phi,i,factors[i].pr);
      for(int j=0;j<factors[i].num_facs;j++)
	printf(" %ld %ld",factors[i].primes[j],factors[i].facs[j]);
      printf("\n");
    }
  */
  return(true);
}

factor factors[MAX_Q];

int main(int argc, char **argv)
{
  _fpu_rndd();
  if(argc!=2)
    print_usage();
  FILE *infile;
  infile=fopen(argv[1],"rb");
  if(!infile)
    fatal_error("Failed to open infile.");
  make_factors(factors,MAX_Q);
  int q_start,q_end;
  unsigned int q_read;
  int_double data[2];
  int_double bound=int_double(2012)/100000.0;
  print_int_double_str("Saad Eddin Bound set to ",bound);
  if(fread(&q_start,sizeof(int),1,infile)!=1)
    fatal_error("Data file corrupt.");
  if(fread(&q_end,sizeof(int),1,infile)!=1)
    fatal_error("Data file corrupt.");
  printf("Processing q from %d to %d\n",q_start,q_end);
  int_double nearest=-1000.0,furthest=0;
  int nearest_q,furthest_q;
  while(true)
    {
      if(fread(&q_read,sizeof(unsigned int),1,infile)!=1)
	break;
      if(fread(data,sizeof(int_double),2,infile)!=2)
	fatal_error("Data file corrupt.");
      //printf("q=%d ",q_read);print_int_double(data[0]);printf(" ");print_int_double(data[1]);printf("\n");
      int_double x=0.5*log(int_double(q_read))-bound;
      //print_int_double_str("|L(1)|=",data[0]);
      //printf("q=%u |L(1)|-0.5 log(q)=%f\n",q_read,x);
      if(-data[0].right>x.left)
	{
	  printf("Even character at q=%d with chi(2)=1 violates Saad-Eddin Bound\n",q_read);
	  print_int_double_str("|L(1,chi)|     =",data[0]);
	  print_int_double_str("1/2log(q)-0.212=",x);
	}
      int_double diff=data[0]-0.5*log(int_double(q_read));
      if(-diff.right>-nearest.right)
	{
	  nearest=diff;
	  nearest_q=q_read;
	}
      if(diff.left<furthest.right)
	{
	  furthest=diff;
	  furthest_q=q_read;
	}
    }
  printf("Nearest was at q=%d with |L(1)|-log(q)/2=",nearest_q);
  print_int_double_str("",nearest);
  printf("Furthest was at q=%d with |L(1)|-log(q)/2=",furthest_q);
  print_int_double_str("",furthest);

  return(0);
}
