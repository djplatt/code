/*
File: L1.cpp

Created: 28th June 2011

Version: <v> = 1.0

Dialect: C++

Requires: -lrt

Implementation notes:
            num_s (int) = 1
            N (int) = 5 no of Taylor terms
            rn (int) must be 2
            NO_GAPS (int)

            N*(NO_GAPS+1) h_rn (dcomplex)

Build instructions: g++ -O1 -msse -march=nocona -fomit-frame-pointer -frounding-math -finline-functions

By: DJ Platt
    Bristol University

Copyright 2008.

This work is funded by the UK ESPRC. */

#define VERSION "1.0"
#define TAY_ERR ((double) 1.6e-20) // Take 5 taylor terms and RN=2

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unistd.h>
#include <semaphore.h>
#include <assert.h>
#include <algorithm>
using namespace std;



#include "../includes/int_double12.0.h"
#include "../includes/im_s.h"
#include "../includes/int-fft-half.h"
#include "L1.h"

void print_usage()
  /* called when wrong arguments passed via command line */
{
  printf("Usage: L1_merg ifile1 ifile2 ofile\n");
  exit(1);
}


void fatal_error(const char *error_string)
  /* print the error message to stdout and exit */
{
  cout << error_string << endl;
  exit(1);
}

void merge_heap(node *h1, node *h2, unsigned long int n, bool comp(node,node))
{
  unsigned long int i;
  for(i=0;i<n;i++)
    if(comp(h2[i],h1[0]))
      {
	/*
	printf("replacing\n");
	print_node(h1[0]);
	printf("with\n");
	print_node(h2[i]);
	*/
	pop_heap(h1,&h1[n],comp);
	h1[n-1]=h2[i];
	push_heap(h1,&h1[n],comp);
      }
}

int main(int argc, char **argv)
{

  _fpu_rndd();

  FILE *ifile1,*ifile2,*ofile;
  long unsigned int n,n1;
  node *h1,*h2;

  if(argc!=4)
    print_usage();
  
  ifile1=fopen(argv[1],"rb");
  if(!ifile1)
    {
      perror("Error opening in_file. ");
      printf("Failed to open file |%s|\n",argv[1]);
      fatal_error("Couldn't open in_file. Exiting.\n");
    }
  ifile2=fopen(argv[2],"rb");
  if(!ifile2)
    {
      perror("Error opening in_file. ");
      printf("Failed to open file |%s|\n",argv[2]);
      fatal_error("Couldn't open in_file. Exiting.\n");
    }

  ofile=fopen(argv[3],"wb");
  if(!ofile)
    {
      perror("Error opening out_file. ");
      printf("Failed to open file |%s|\n",argv[3]);
      fatal_error("Couldn't open in_file. Exiting.\n");
    }

  fread(&n,sizeof(long unsigned int),1,ifile1);
  fread(&n1,sizeof(long unsigned int),1,ifile2);
  if(n!=n1)
    {
      printf("Size of ifile1=%lu\nSize of ifile2=%lu\n",n,n1);
      fatal_error("File size mismatch");
    }
  fwrite(&n,sizeof(long unsigned int),1,ofile);
  h1=(node *)malloc(sizeof(node)*n);
  h2=(node *)malloc(sizeof(node)*n);
  fread(h1,sizeof(node),n,ifile1);
  fread(h2,sizeof(node),n,ifile2);
  printf("merging largest L(1)/log(q) based on left end point.\n"); 
  merge_heap(h1,h2,n,llcomp); // largest L(1)/log(q) based on left end point
  fwrite(h1,sizeof(node),n,ofile);
  fread(h1,sizeof(node),n,ifile1);
  fread(h2,sizeof(node),n,ifile2);
  printf("merging largest L(1)/log(q) based on right end point.\n"); 
  merge_heap(h1,h2,n,lrcomp); // largest L(1)/log(q) based on right end point
  fwrite(h1,sizeof(node),n,ofile);
  fread(h1,sizeof(node),n,ifile1);
  fread(h2,sizeof(node),n,ifile2);
  printf("merging smallest L(1)*log(q) based on left end point.\n"); 
  merge_heap(h1,h2,n,hlcomp); // smallest L(1)*log(q) based on left end point
  fwrite(h1,sizeof(node),n,ofile);
  fread(h1,sizeof(node),n,ifile1);
  fread(h2,sizeof(node),n,ifile2);
  printf("merging smallest L(1)*log(q) based on right end point.\n"); 
  merge_heap(h1,h2,n,hrcomp); // smallest L(1)*log(q) based on right end point
  fwrite(h1,sizeof(node),n,ofile);


  return(0);
}
