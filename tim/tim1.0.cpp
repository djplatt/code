//
// File: tmi.cpp
//
// Code to do Tim Trudgian's counting
//
// DJ Platt 
// Setup
//
// we are given a positive odd integer d>1 and 
// we have K variables, v_1 to v_k.
// We let the v_i range from 1 to ord_d(2) and count
// how many ways we can get each of the d possible values for N 
// N=2^{v_1}+2^{v_2}+...+2^{v_K} mod d
// We return the maximum found.
//
#include "stdlib.h"
#include "stdio.h"
#include "malloc.h"
#include "inttypes.h"

void print_W(uint64_t *W, uint64_t d)
{
  for(uint64_t i=0;i<d;i++)
    printf("W[%lu]=%lu\n",i,W[i]);
}

uint64_t sum_W(uint64_t *W, uint64_t d)
{
  uint64_t tot=0;
  for(uint64_t i=0;i<d;i++)
    tot+=W[i];
  return(tot);
}

// rotate W_spare right by r places and add it into W component wise
// r is in [0,d)
void rotate_and_add(uint64_t *W, uint64_t *W_spare, uint64_t r, uint64_t d)
{
  uint64_t ptr=d-r;
  for(uint64_t i=0;i<d;i++)
    {
      W[i]+=W_spare[ptr++];
      if(ptr==d)
	ptr=0;
    }
}

void doit(uint64_t k, uint64_t d, uint64_t ed)
{
  uint64_t *W,*Wk;
  W=(uint64_t *)malloc(sizeof(uint64_t)*d);
  Wk=(uint64_t *)malloc(sizeof(uint64_t)*d);
  for(uint64_t i=0;i<d;i++)
    W[i]=0;

  // set up W_1

  for(uint64_t i=1,ptr=2*k;i<=ed;i++)
    {
      ptr=(ptr+(1<<i))%d;
      W[ptr]++;
    }

  // now do W_2..W_K
  for(uint64_t i=2;i<=k;i++)
    {
      // copy W to Wk
      for(uint64_t j=0;j<d;j++)
	Wk[j]=W[j];
      //
      // now rotate Wk by 2,6,14,30... and add into W
      //
      for(uint64_t rl=2,drl=4,e=1;e<ed;e++)
	{
	  rotate_and_add(W,Wk,rl,d);
	  rl+=drl;
	  rl%=d;
	  drl<<=1;
	  drl%=d;
	}
    }

  // finished, so print some results
  uint64_t wc=0;
  wc=W[0];
  printf("W[0]=%lu\n",W[0]);
  for(uint64_t i=1;i<d;i++)
    {
      printf("W[%lu]=%lu\n",i,W[i]);
      if(W[i]<wc) wc=W[i];
    }
  printf("Min was %lu\n",wc);

}

// compute smallest n>0 s.t. 2^n = 1 mod d
// don't call this with even d!
uint64_t ord2(uint64_t d)
{
  uint64_t a=2,res=1;
  while(a!=1)
    {
      a<<=1;res++;
      if(a>d)
	a-=d;
    }
  return(res);
}

int main(int argc, char **argv)
{

  if(argc!=3)
    {
      printf("Usage:- %s <K> <d>\n",argv[0]);
      exit(0);
    }

  if(atol(argv[1])<1)
    {
      printf("K must be > 0. Exiting.\n");
      exit(0);
    }

  uint64_t k=atol(argv[1]);
  uint64_t d=atol(argv[2]);
  if(!(d&1))
    {
      printf("d must be odd. Exiting.\n");
      exit(0);
    }
  if(d<3)
    {
      printf("d must be >1. Exiting.\n");
      exit(0);
    }
  uint64_t ed=ord2(d);
  printf("Running with K=%lu d=%lu and Ord_%lu(2)=%lu.\n",k,d,d,ed);
  doit(k,d,ed);
  return(0);
}
