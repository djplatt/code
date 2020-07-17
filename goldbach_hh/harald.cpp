/*
File: int-l-func5.10.cpp

Created: 10th Jan 2009

Version: <v> = 5.10

Last Modified: 18th June 2012

5.0 changed summation code in calc_rn and calc_rn_neg
5.1 now multiply by (q/pi)^(it/2)
5.2 moved (q/pi)^(it/2) and q^(-s) into FFT
5.3 changed calc of q_comb for accuracy, introduced im_s_vec[].pi_minus_it_2
5.4 using semaphores to access lattice files
5.5 read lattice file in one go using read
5.6 improved calculation of hurwitz values
5.7 now uses crlibm 
5.8 makes its own factors, no longer uses semaphores
5.9 now includes int_fft1.2.h
    (this computes DFT of b once only for Bluestein convolution)
    thus no longer need b_spare
5.10 Exploit fact that phi_q is always even to halve FFT lengths.
     now includes int_fft1.3.h

Dialect: C++

Requires: -lrt

Implementation notes: Takes Output from l-func-mpfi-1.1
            num_s (int)
            N (int)
            rn (int) must be 1 in this implementation
            NO_GAPS (int)

            im_s (double)
            gam_s (dcomplex)
            gam_s_2 (dcomplex)
            N*(NO_GAPS+1) h_rn (dcomplex)

Build instructions: g++ -O1 -mfpmath=387 -fomit-frame-pointer -frounding-math -finline-functions

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "5.10"


#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <semaphore.h>
//#include <assert.h>
using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft1.3.h"
/*
extern "C" {
#include "../includes/nit.h"
}
*/
#include "../includes/upsamdefs.h"
#include "../includes/qit_struct.h"
//#include "../includes/qit.h"


//simple Euler Mac. calculation of zeta(s,alpha)
//cache s_mod,sigma s_array*bernoulli etc.
int_complex hurwitz_half (const double t, const int_double &alpha)
{
  unsigned int i,k=MAX_BERNOULLI_K,n=DEFAULT_HUR_N;
  int_double err,n_alpha;
  int_complex s=int_complex(0.5,t),s1=s,res=c_zero,n_alpha_s,term;
  int_complex s_array[MAX_BERNOULLI_K];

  if(!h_bernoulli_initialised)
    {
      set_h_bernoulli();
      h_bernoulli_initialised=true;
    }
//debug;
  if(n<t)
    n=t/2;
  double s_mod=((t<=0.5)?0.5:t);
  //print_int_double_str("alpha=",alpha);
  n_alpha=alpha+n;
  //debug;
  //print_int_double_str("n+alpha=",n_alpha);
  //print_int_complex_str("-s+1=",-s+1);
  n_alpha_s=pow(n_alpha,-s+1);
  //debug;
  for(i=0;i<n;i++)
    {
      res=res+pow(alpha+i,-s);
      //      if((i%1000)==0)
      //	printf("relative error in res is %d %d\n",rel_error(res.real),rel_error(res.imag));
    }

  //debug;
  res=res+n_alpha_s/(s-1);

  n_alpha_s=n_alpha_s/n_alpha; // n_alpha_s=(n+alpha)^(-s)
		
  res=res+n_alpha_s/2;

  s_array[0]=s*n_alpha_s/n_alpha;

  for(i=1;i<k;i++)
    {
      s1=s1+1;
      s_array[i]=s_array[i-1]*s1/n_alpha;
      s1=s1+1;
      s_array[i]=s_array[i]*s1/n_alpha;
    }
  
  for(i=0;i<k;i++)
    {
      term=s_array[i]*h_bernoulli[i];
      res=res+term;
    }

  err=sqrt(norm(term))/((s.real+2*k-1)*(s_mod+2*k-1));
  if(err.left<=err.right)
    err.right=err.left;
  else
    err.left=err.right;
  print_int_double_str("error term is ",err);
  res.real=res.real+err;
  res.imag=res.imag+err;
  return(res);
}

#define num_tests (10)
#define height ((double) 10000.0)

int main(int argc, char **argv)
{

  // initialsie the int_double package
  _fpu_rndd();

  int_complex res[num_tests];
  for(unsigned long int i=0;i<num_tests;i++)
    res[i]=hurwitz_half(height,int_double(i)/1000);
  return(0);
}
