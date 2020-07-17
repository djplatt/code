/*
File: L1_print.cpp

Created: 4th October 2011

Version: <v> = 1.0

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
  printf("Usage: L1.2_print (ifname)\n");
  printf("  (ifname)    - file with output from L1.2.\n");
  exit(1);
}


void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
}

int main(int argc, char **argv)
{
  if(argc!=2)
    print_usage();

  _fpu_rndd();
  FILE *infile;
  infile=fopen(argv[1],"rb");
  if(!infile)
    fatal_error("Failed to open infile.");
  //  make_factors(factors,MAX_Q);
  // int q_start,q_end;
  unsigned int q,q_read;
  int_double data[2];
  while(true)
    {
      if(fread(&q_read,sizeof(unsigned int),1,infile)!=1)
	break;
      if(fread(data,sizeof(int_double),2,infile)!=2)
	fatal_error("Data file corrupt.");
      printf("q=%u ",q_read);print_int_double(data[0]);print_int_double_str(" ",data[1]);
    }
  //printf("Largest Odd |L(1)|-1/2log q at q=%u.",largest_q);print_int_double_str("",largest_x);
  return(0);
}
