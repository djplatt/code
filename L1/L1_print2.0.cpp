/*
File: L1_print2.0.cpp

Created: 4th October 2011

Version: <v> = 2.0
       : Analysis of L(1) as per Ramare email 5/11/11
       : Input from L1.1.cpp

Dialect: C++

Requires: -lrt

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
#include "../includes/int-fft-half.h"

void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: L1_print (ifname) (max_q)\n");
  printf("  (ifname)    - file containing list of files with output from L1.1.\n");
  printf("  (max_q)     - largest modulus in database.\n");
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
  //return(true);
  //printf("Making factor database.\n");
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
  //printf("Prime sieve completed.\n");
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
  //printf("Factors computed.\n");
  // now calculate phi(f)
  for(int i=3;i<=q_end;i++)
    {
      factors[i].phi=1;
      for(int j=0;j<factors[i].num_facs;j++)
	factors[i].phi*=phi(factors[i].primes[j],factors[i].facs[j]);
    }
  //printf("phi computed.\n");

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
  //printf("pr's computed.\n");
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


#define NO_C1 (41)
#define C1_GAP ((double) 0.1)
#define C1_RECIP ((double) 10.0)
#define C1_MAX ((double) 1.0)

typedef struct{
  unsigned long int odd_c1s[NO_C1];
  unsigned long int even_c1s[NO_C1];
} ramare;

int main(int argc, char **argv)
{
  if(argc!=3)
    print_usage();
  FILE *infile,*infile1;
  infile1=fopen(argv[1],"r");
  if(!infile1)
    fatal_error("Failed to open infile.");
  long int max_q=atol(argv[2]);
  if(max_q<3)
    fatal_error("max_q must be at least 3.");

  double c0s[7]={1.0/6.0,0.5,0.25,1.0/3.0,0.25,0.5,0.5};
  ramare ramares[7];
  for(long unsigned int r=0;r<7;r++)
    for(long unsigned int c1=0;c1<NO_C1;c1++)
      {
	ramares[r].odd_c1s[c1]=0;
	ramares[r].even_c1s[c1]=0;
      }

  factor *factors;
  factors=(factor *) malloc(sizeof(factor)*(max_q+1));
  if(!(factors&&ramares))
    fatal_error("Failed to allocate memory for factors/ramares.");
  //printf("Making factors.\n");
  make_factors(factors,max_q);
  //printf("Factor base made.\n");
  int q_start,q_end;
  unsigned int q,q_read;
  int_double data[4];

  char fname[1024];
  while(true)
    {
      fscanf(infile1,"%s",&fname);
      if(feof(infile1))
	break;
      //printf("Processing file %s\n",fname);
      infile=fopen(fname,"rb");
      if(!infile)
	fatal_error("Failed to open infile.\n");
      fread(&q_start,sizeof(int),1,infile);
      fread(&q_end,sizeof(int),1,infile);
      //printf("q_start=%lu q_end=%lu\n",q_start,q_end);
      while(true)
	{
	  if(fread(&q_read,sizeof(unsigned int),1,infile)!=1)
	    break;
	  if(fread(data,sizeof(int_double),4,infile)!=4)
	    fatal_error("Data file corrupt.");

	  unsigned long int m6;
	  if(factors[q_read].primes[0]==q_read) // prime q
	    m6=6;
	  else
	    m6=q_read%6;
	  double lnq=c0s[m6]*log((double) q_read);

	  if(q_read!=12) // no odd primitive characters
	    {
	      double x=-data[0].right-lnq; // largest odd
	      long unsigned int n=0;
	      for(double c1=C1_MAX;c1>=x;c1-=C1_GAP,n++);
	      for(;n<NO_C1;n++)
		ramares[m6].odd_c1s[n]=q_read;
	    }
	  
	  if(q_read!=4)
	    {
	      double x=-data[2].right-lnq; // largest even
	      long unsigned int n=0;
	      for(double c1=C1_MAX;c1>=x;c1-=C1_GAP,n++);
	      for(;n<NO_C1;n++)
		ramares[m6].even_c1s[n]=q_read;
	    }
	}
      fclose(infile);
    }
  /*
  printf("          ");
  for(long unsigned int c1=0;c1<NO_C1;c1++)
    printf("%7.1f ",(c1-(NO_C1-1)/2.0)*C1_GAP);
  for(long unsigned int r=0;r<7;r++)
    {
      printf("\nr=%1lu odd : ",r);
      for(long unsigned int c1=0;c1<NO_C1;c1++)
	printf("%7lu ",ramares[r].odd_c1s[c1]);
      printf("\n   even : ");
      for(long unsigned int c1=0;c1<NO_C1;c1++)
	printf("%7lu ",ramares[r].even_c1s[c1]);
    }
  */
  printf("\nLargest q<=%lu for which any primitive |L(1)|>C_0 log q + C_1\n",max_q);
  printf("Note that q which are both prime and =1,5 mod 6 appear in prime column.\n\n");
  printf("q mod 6  : ");
  for(long unsigned int r=0;r<6;r++)
    printf("|       %1lu       ",r);
  printf("|     Prime\n");
  printf("C_0      : ");
  for(long unsigned int r=0;r<7;r++)
    printf("|     %5.3f     ",c0s[r]);
  printf("\nParity   : ");
  for(long unsigned int r=0;r<7;r++)
    printf("|  Odd  |  Even ");
  for(long unsigned int c1=0;c1<NO_C1;c1++)
    {
      printf("\nC_1=%+3.1f : ",C1_MAX-c1*C1_GAP);
      for(long unsigned int r=0;r<7;r++)
	printf("|%7lu|%7lu",ramares[r].odd_c1s[c1],ramares[r].even_c1s[c1]);
    }
  printf("\n");
  //printf("Largest Odd |L(1)|-1/2log q at q=%u.",largest_q);print_int_double_str("",largest_x);
  return(0);
}
