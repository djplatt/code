/*
File: theta.cpp

Created: 19th March 2015

Version: <v> = Original

Dialect: C++

Requires: -lrt

Build instructions: g++ -fomit-frame-pointer -frounding-math -finline-functions -lcrlibm

By: DJ Platt
    Heilbronn Institute
    Bristol University

Copyright 2015.


The author is a Heilbronn Research Fellow */

#define VERSION ""
#define MAX_X (40000000)

#include <stdio.h>
#include <stdlib.h>
//#include <iostream>
//#include <fstream>
//#include <time.h>
//#include <unistd.h>
#include <inttypes.h>

using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/factors.h"
//#include "../includes/hurwitz1.0.h"
//#include "../includes/im_s.h"
//#include "../includes/int-fft-half.h"
//#include "L1.h"

void print_usage()
  /* called when wrong arguments passed via command line */
{
  exit(1);
}


void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
};

int main(int argc, char **argv)
{

  _fpu_rndd();

  int64_t i;
  factor *factors;

  factors=make_factors(MAX_X);

  int_double theta[MAX_X+1];
  int_double logs[MAX_X+1];
  int_double sqrts[MAX_X+1];
  sqrts[1]=1;
  theta[1]=0.0;
  for(uint64_t i=2;i<=MAX_X;i++)
    {
      int_double id=i;
      logs[i]=log(id);
      sqrts[i]=sqrt(id);
      if(i==factors[i].primes[0]) // a prime
	theta[i]=logs[factors[i].primes[0]];
      else
	theta[i]=0.0;
    }

#define n_q (6)
  uint64_t qs[n_q]={3,5,7,11,13,17};
  for(uint64_t qq=0;qq<n_q;qq++) 
    {
      uint64_t q=qs[qq];
      uint64_t q2=q*q;
      uint64_t phi_q2=q*(q-1);
      int_double phi_q2d=1.0;
      phi_q2d/=phi_q2;
      double worst=0.0;
      uint64_t worst_a =0, worst_x=0;
      int_double worst_theta=0.0;
      for(uint64_t a=1;a<q2;a++)
	{
	  if(!co_prime(a,q2))
	    continue;
	  int_double theta_sum=0.0;
	  for(uint64_t y=a;y<=MAX_X;y+=q2)
	    if(theta[y].left!=0.0) // this is a prime
	      {
		int_double diff1=theta_sum-(y-1)*phi_q2d; // - y/phi(q^2)
		if(diff1.left<0.0) diff1=-diff1;
		diff1/=sqrts[y-1];
		if(-diff1.right>worst)
		  {
		    worst=-diff1.right;worst_a=a;worst_x=y-1;worst_theta=theta_sum;
		    printf("New record at q=%lu a=%lu just before y=%lu %20.18e\n",q2,a,y,worst);
		  }
		theta_sum+=logs[y];
		diff1=theta_sum-y*phi_q2d; // - y/phi(q^2)
		if(diff1.left<0.0) diff1=-diff1;
		diff1/=sqrts[y];
		if(-diff1.right>worst)
		  {
		    worst=-diff1.right;worst_a=a;worst_x=y-1;worst_theta=theta_sum;
		    printf("New record at q=%lu a=%lu just after  y=%lu %20.18e\n",q,a,y,worst);
		  }
	      }
	}
      printf("worst found for q=%lu was %20.18e at a=%lu x=%lu",q2,worst,worst_a,worst_x);
      print_int_double_str(" theta was ",worst_theta);

    }
  return(0);
}
