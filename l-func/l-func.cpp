/*

File: l-func1.0.cpp

Created: 27th June 2008

Version: <v> = 1.0

Last Modified: 27th June 2008

Dialect: C++

Requires: No Libraries required

Implementation notes: Might get away with 1/2 of the fracs array
                      if its symmetric about 1/2.

Build instructions: g++ -ol-func<v> l-func<v>.cpp -O2

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "1.0"










#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define Q (10000)   /* maximum conductor */
#define NUM_FRACS (30397485)  /* number of distint Farey fractions
				 0 < p/q < Q with q <=Q           */
#define SUCCESS (1)
#define FAILURE (0)

using namespace std;


void print_usage()
{
	printf("Usage: l-func%s \n",VERSION);
	exit(FAILURE);
};

void fatal_error(char *error_string)
{
  cout << error_string << endl;
  exit(FAILURE);
};

typedef struct {double val; 
  unsigned short int num; unsigned short int den;} frac;

inline unsigned short int gcd (unsigned short int a, unsigned short int b)
{
  unsigned short int c;
  while(a!=0)
    {
      c=a;
      a=b%a;
      b=c;
    };
  return(b);
};

inline int co_prime(unsigned short int a, unsigned short int b)
{
  return(gcd(a,b)==1);
};

frac fracs[NUM_FRACS];

/* comparison function for qsort */
inline int _frac_comp (frac a, frac b)
{
  if(a.val-b.val>0)
    return(1);
  return(-1);
};

int frac_comp(const void *a,const void *b)
{
  return(_frac_comp(*(frac*)a,*(frac*)b));
};  

int make_fracs()
{
  unsigned short int num,den;
  unsigned int ptr,ptr2;

  cout << "Creating Farey fractions." << endl;
  ptr=0;
  for(den=1;den<=Q;den++)
    for(num=1;num<den;num++)
      if(co_prime(den,num))
	{
	  if(ptr==NUM_FRACS)
	    fatal_error(
	      "Generated more Farey fractions than expected. Exiting");
	  fracs[ptr].val=((double) num)/((double) den);
	  fracs[ptr].num=num;
	  fracs[ptr++].den=den;
	};

  cout << "Sorting." << endl;

  qsort(fracs,ptr,sizeof(frac),frac_comp);


  for(ptr2=1;ptr2<ptr;ptr2=ptr2<<1)
    cout << fracs[ptr2].val << endl;

  return(SUCCESS);
};

int main(int argc, char **argv)
{
  if(argc!=1)
    print_usage();

  if(!make_fracs())
    fatal_error("Error making Farey Fractions. Exiting.");
  
 return(SUCCESS);
};

