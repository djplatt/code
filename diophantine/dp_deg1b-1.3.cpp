//
// File:- dp_deg1b-1.3.cpp
//
// Author: DJ Platt
//
// Written: 23 January 2013
//
// Solutions to Del Pezzo surfaces of degree 1
//
// want w,x,y,z in Z>0
// (w,x,y,z)=1
// x,y<=H
// z<=H^2
// w<=H^3
//
// and
//
// w^2-ax^6=by^6-z^3
//
// we have two heaps, lnodes and rnodes
// we store w,x in H lnodes .a and .b
// initially a=0, b=0..H-1 (meaning 1 and 1..H)
// pq=w^2-ax^6
//
// we store y,z in H rnodes .a and .b
// initially a=0..H, b=0
// pq=by^6-z^3
//
// we compare the top of two heaps.
// if l.pq==r.pq we have a solution
// if l.pq<r.pq increase l.a
// if l.pq>r.pq increase r.b
//
// once l.a,r.b exceed H^3,H^2 respectively
// pop that node.
// continue until heap is empty

#include "inttypes.h"
#include "stdio.h"
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
//#include <assert.h>
#include <math.h>

long int H,w_H,z_H; // runs to height<=H

// 128 bit unsigned built in to GCC (efficiency?)
typedef __int128_t bigint;

// a node in our heap
// allign to 16 byte boundary
// so we can use fast SSE instructions
typedef struct{
  //unsigned int node_no;
  __attribute__ ((aligned(16)))
  long int a;
  long int b;
  bigint pq;
} node;

void print_bigint_pos(bigint i)
{
  if(i<10)
    printf("%1d",i);
  else
    {
      print_bigint_pos(i/10);
      printf("%1d",i%10);
    }
}

void print_bigint(bigint i)
{
  if(i<0)
    {
      printf("-");
      i=-i;
    }
  print_bigint_pos(i);
}

void print_node(node *n)
{
  printf("[%ld,%ld,",n->a+1,n->b+1);
  print_bigint(n->pq);
  printf("]\n");
}

void print_node1(node *n)
{
  printf("[%ld,%ld,",n->a+1,-(n->b+1));
  print_bigint(n->pq);
  printf("]\n");
}

void print_heap(node *heap, int last_one)
{
  for(int i=0;i<=last_one;i++)
    print_node(heap+i);
}

inline long int gcd (long int a, long int b)
/* Euclid algorithm gcd */
{
	int c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

inline long int gcd (long int a, long int b,
			      long int c)
{
  return(gcd(gcd(a,b),c));
}

inline long int gcd (long int a, long int b,
			      long int c, long int d)
{
  return(gcd(gcd(a,b),gcd(c,d)));
}

#define max(a,b) (a>b ? a : b)

long int solution_count=0;

long int sq(long int i) {return (i*i);}
long int cb(long int i) {return (i*sq(i));}
long int sx(long int i) {return (sq(cb(i)));}

double height(node *n1, node *n2)
{
  double res=(double) max(max(sq(n1->a+1),sx(n1->b+1)),max(sx(n2->a+1),cb(n2->b+1)));
  return(exp(log(res)/6.0));
}


void print_solution(node *n1, node *n2)
{
  if(gcd(n1[0].a+1,n1[0].b+1,n2[0].a+1,n2[0].b+1)==1)
    {
      printf("solution %ld found\n",++solution_count);
      print_node(n1);
      print_node1(n2);
      printf("Height= %8f.\n",height(n1,n2));
      //printf("Terminating search.\n");
      //exit(0);
    }
}

// res<-a*x^4
inline bigint p(long int x, long int a)
{
  bigint res=x;
  res*=res;
  res*=res;
  res*=a;
  return(res);
}

// returns 0 if top <left and top < right
// returns 1 if top > left < right
// returns -1 if top > right < left
inline int comp_nodes(node top, node left, node right)
{
  if(left.pq<=right.pq) // left <= right
    {
      if(top.pq<=left.pq) // top <= left <= right
	return(0);
      else // top > left and top <= right
	return(1); // so swap left
    }
  // right < left
  if(top.pq<=right.pq) // top <= right < left
    return(0);
  else
    return(-1);
} 

// use the 128 bit XMM registers to swap two nodes
#define swap_node(x,y)						\
  {								\
    __asm(							\
	  "movapd %0,%%XMM0\n\t"				\
	  "movapd %1,%%XMM1\n\t"				\
	  "movapd %%XMM0,%1\n\t"				\
	  "movapd %%XMM1,%0\n\t"				\
	  "movapd %2,%%XMM0\n\t"				\
	  "movapd %3,%%XMM1\n\t"				\
	  "movapd %%XMM0,%3\n\t"				\
	  "movapd %%XMM1,%2"					\
	  : "=m" (x), "=m" (y), "=m" (x.pq), "=m" (y.pq)	\
	  :							\
	  : "xmm0","xmm1");					\
  }

inline void balance_heap (node *heap, int heap_end)
{
  int this_ptr=0,left_ptr,right_ptr,cmp;
  while(true)
    {
      left_ptr=(this_ptr<<1)+1;
      right_ptr=(this_ptr+1)<<1;
      if(left_ptr<=heap_end) // it has a left child
	{
	  if(right_ptr<=heap_end) // it has a left and right children
	    {
	      cmp=comp_nodes(heap[this_ptr],heap[left_ptr],heap[right_ptr]);
	      //printf("comp of %d %d %d returned %d\n",this_ptr,left_ptr,right_ptr,cmp);
	      if(cmp>0) // swap with left
		{
		  swap_node(heap[this_ptr],heap[left_ptr]);
		  this_ptr=left_ptr;
		  continue;
		}
	      if(cmp<0) // swap with right
		{
		  swap_node(heap[this_ptr],heap[right_ptr]);
		  this_ptr=right_ptr;
		  continue;
		}
	      // parent node is in right place so stop
	      return;
	    }
	  else // this node only has a left child
	    {
	      if(heap[left_ptr].pq<heap[this_ptr].pq) // left<this
		swap_node(heap[this_ptr],heap[left_ptr]);
	      return; // tree is balanced
	    }
	}
      else // this node is at bottom
	return;
    }
}

#define pop_node(heap,heap_end){heap[0]=heap[heap_end];\
    heap_end--;\
balance_heap(heap,heap_end);}

void check_eqn(node *lnodes, node *rnodes,
	       bigint *squares, bigint *cubes, bigint *as, bigint *bs)
{
  int cmp;
  bigint temp,temp1;

  // set up heaps
  for(long int i=0;i<H;i++)
    {
      lnodes[i].a=0;
      lnodes[i].b=H-i-1;
      lnodes[i].pq=1-as[H-1-i];
    }
  int left_end=H-1;
  for(long int i=0;i<H;i++)
    {
      rnodes[i].a=i;
      rnodes[i].b=z_H-1;
      rnodes[i].pq=bs[i]-cubes[z_H-1];
    }
  int right_end=H-1;

  while(true)
    {
      //printf("Iterating with left=");print_node(lnodes);
      //printf("          and right=");print_node(rnodes);
      if(lnodes[0].pq<rnodes[0].pq) // left<right
	{
	  if(lnodes[0].a<w_H-1) // more left nodes to add
	    {
	      lnodes[0].a++;
	      lnodes[0].pq=squares[lnodes[0].a]-as[lnodes[0].b];
	      //printf("L Pushing ");print_node(lnodes[0]);
	      balance_heap(lnodes,left_end);
	      //print_heap(lnodes,H-1);
	    }
	  else // no more left nodes, so just pop
	    {
	      pop_node(lnodes,left_end);
	      if(left_end<0) // heap empty
		{
		  //print_node(rnodes[0]);
		  printf("Heap max reached.\n");
		  return;
		}
	    }
	  continue;
	}
      if(lnodes[0].pq>rnodes[0].pq) // right<left
	{
	  if(rnodes[0].b>0)
	    {
	      rnodes[0].b--;
	      rnodes[0].pq=bs[rnodes[0].a]-cubes[rnodes[0].b];
	      //printf("R Pushing ");print_node(rnodes[0]);
	      balance_heap(rnodes,right_end);
	      //print_heap(rnodes,H-1);
	    }
	  else
	    {
	      pop_node(rnodes,right_end);
	      if(right_end<0)
		{
		  //print_node(lnodes[0]);
		  printf("Heap max reached.\n");
		  return;
		}
	    }
	  continue;
	}
      // right=left
      //if(gcd(lnodes[0].a,lnodes[0].b,rnodes[0].a,rnodes[0].b)==1)
      print_solution(lnodes,rnodes);//exit(0);
      //printf("Heaps contain...\n");
      //print_heap(lnodes,left_end);
      //printf("\n");
      //print_heap(rnodes,right_end);
      bool left_duplicate=false;
      if(left_end>0)
	{
	  if(lnodes[0].pq==lnodes[1].pq)
	    left_duplicate=true;
	  else
	    {
	      if(left_end>1)
		left_duplicate=(lnodes[0].pq==lnodes[2].pq);
	    }
	}
      bool right_duplicate=false;
      if(right_end>0)
	{
	  if(rnodes[0].pq==rnodes[1].pq)
	    right_duplicate=true;
	  else
	    {
	      if(right_end>1)
		right_duplicate=(rnodes[0].pq==rnodes[2].pq);
	    }
	}
      if(left_duplicate)
	{
	  if(right_duplicate) // this has got to be rare?
	    {
	      printf("Fatal error ");
	      print_bigint(lnodes[0].pq);
	      printf(" can be formed 2 or more ways on both sides. Exiting.\n");
	      exit(0);
	    }
	  else
	    {
	      if(lnodes[0].a<w_H-1) // more left nodes to add
		{
		  lnodes[0].a++;
		  lnodes[0].pq=squares[lnodes[0].a]-as[lnodes[0].b];
		  //printf("L Pushing ");print_node(lnodes[0]);
		  balance_heap(lnodes,left_end);
		  //print_heap(lnodes,H-1);
		}
	      else // no more room on this lnode, so pop it
		pop_node(lnodes,left_end);
	    }
	}
      else
	{
	  if(right_duplicate)
	    {
	      if(rnodes[0].b>0) // more right nodes to add
		{
		  rnodes[0].b--;
		  rnodes[0].pq=bs[rnodes[0].a]-cubes[rnodes[0].b];
		  //printf("L Pushing ");print_node(lnodes[0]);
		  balance_heap(rnodes,right_end);
		  //print_heap(lnodes,H-1);
		}
	      else // no more room on this lnode, so pop it and move RHS on
		pop_node(rnodes,right_end);
	    }
	  else
	    {
	      if(lnodes[0].a<w_H-1) // more left nodes to add
		{
		  lnodes[0].a++;
		  lnodes[0].pq=squares[lnodes[0].a]-as[lnodes[0].b];
		  //printf("L Pushing ");print_node(lnodes[0]);
		  balance_heap(lnodes,left_end);
		  //print_heap(lnodes,H-1);
		}
	      else // no more room on this lnode, so pop it
		{
		  pop_node(lnodes,left_end);
		  if(left_end<0)
		    printf("Heap max reached.\n");
		  return;
		}
	    }
	}
    }
}

int main(int argc, char **argv)
{
  node *lnodes,*rnodes,left_node,right_node;
  bigint *squares,*cubes,*as,*bs;
  int cmp;
  long int As[2]; // this defines the equation

  if(argc!=4)
    {
      printf("Incorrect command line:- %s H a b.\n",argv[0]);
      exit(0);
    }

  H=atoi(argv[1]);
  z_H=H*H;
  w_H=H*z_H;
  for(int i=0;i<2;i++)
    As[i]=atoi(argv[i+2]);

  printf("RHS max H set to %ld.\n",H);

  lnodes=(node *) memalign(16,H*sizeof(node));
  rnodes=(node *) memalign(16,H*sizeof(node));
  squares=(bigint *) malloc(w_H*sizeof(bigint));
  cubes=(bigint *) malloc(z_H*sizeof(bigint));
  as=(bigint *) malloc(H*sizeof(bigint));
  bs=(bigint *) malloc(H*sizeof(bigint));
  if(!(lnodes&&rnodes&&squares&&cubes&&as&&bs))
    {
      printf("Problem allocating memory. Exiting.\n");
      exit(0);
    }

  printf("Looking for solutions to w^2-%ldx^6=%ldy^6-z^3 in (Z>0)^4 to height %ld\n",
	 As[0],As[1],H);
  for(int i=0;i<H;i++)
    {
      squares[i]=i+1;squares[i]*=i+1;
      cubes[i]=squares[i]*(i+1);
      bigint sixth=cubes[i]*cubes[i];
      as[i]=As[0]*sixth;
      bs[i]=As[1]*sixth;
    }
  for(int i=H;i<z_H;i++)
    {
      squares[i]=i+1;squares[i]*=i+1;
      cubes[i]=squares[i]*(i+1);
    }
  for(int i=z_H;i<w_H;i++)
    {
      squares[i]=i+1;squares[i]*=i+1;
    }
      
  /*
  printf("bigint is of length %d\n",sizeof(bigint));
  bigint foo=1;
  for (int i=1;i<=128;i++,foo<<=1) {printf("2^%d=",i);print_bigint(foo);printf("\n");} exit(0);
  print_bigint(as[H-1]);exit(0);
  */
  check_eqn(lnodes,rnodes,squares,cubes,as,bs);
}
