/*
File: int-l-func5.5.cpp

Created: 10th Jan 2009

Version: <v> = 5.8

Last Modified: 16th September 2009

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


//int_complex a[32],bs[32],b_star_conjs[32];

int main()
{
  _fpu_rndd();
  unsigned int conv_sizes[1],i,lim=6;
  int_complex res[16];
  init_ws(_ws);
  printf("Initialising FFT\n");
  init_bluestein_fft(lim,conv_sizes,bs,b_star_conjs);
  printf("Conv size set to %d\n",conv_sizes[0]);
  for(i=0;i<lim;i++)
    a[i]=int_double(i*i);
  printf("Running FFT.\n");
  bluestein_fft(1,a,lim,a,bs,b_star_conjs,conv_sizes[0]);
  for(i=0;i<lim;i++)
    print_int_complex_str("a[]=",a[i]);
  lim=5;
  printf("Initialising FFT\n");
  init_bluestein_fft(lim,conv_sizes,bs,b_star_conjs);
  printf("Conv size set to %d\n",conv_sizes[0]);
  for(i=0;i<lim;i++)
    a[i]=int_double(i*i);
  printf("Running FFT.\n");
  bluestein_fft(1,a,lim,a,bs,b_star_conjs,conv_sizes[0]);
  for(i=0;i<lim;i++)
    print_int_complex_str("a[]=",a[i]);

  lim=6;
  printf("Initialising FFT\n");
  init_bluestein_fft(lim/2,conv_sizes,bs,b_star_conjs);
  printf("Conv size set to %d\n",conv_sizes[0]);
  for(i=0;i<lim;i++)
    res[i]=int_double(i*i);
  printf("Running FFT.\n");
  bluestein_fft(2,res,lim/2,a,bs,b_star_conjs,conv_sizes[0]);
  bluestein_fft(2,res+1,lim/2,a,bs,b_star_conjs,conv_sizes[0]);
  for(i=0;i<lim;i++)
    print_int_complex_str("res[]=",res[i]);
  printf("\n");
  int_complex ws[16];
  for(i=0;i<lim/2;i++)
    {
      sin_cospi(-d_one/(lim/2)*i,&ws[i].imag,&ws[i].real);
      print_int_complex_str("ws[]=",ws[i]);
    }
  printf("\n");
  for(i=0;i<lim;i++)
    {
      if(i<lim/2)
	print_int_complex_str("fft[]=",res[i+i]+ws[i]*res[i+i+1])
      else
	print_int_complex_str("fft[]=",res[i+i-lim]-ws[i-lim/2]*res[i+i-lim+1]);
    }

  return(0);
}
