// version 2.1 use SSE instructions to move nodes
// stop after first solution
//
#include "inttypes.h"
#include "stdio.h"
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <assert.h>
#include "math.h"

// 128 bit unsigned built in to GCC (efficiency?)
typedef __uint128_t bigint;

// a node in our heap
typedef struct{
  //unsigned int node_no;
  __attribute__ ((aligned(16)))
  long unsigned int x;
  long unsigned int y;
  bigint n; // 4x^2+y^2
} node;


void print_bigint(bigint i)
{
  printf("%ju",(uintmax_t) i);
}

void print_node(node n)
{
  printf("[%lu,%lu,",n.x,n.y);
  print_bigint(n.n);
  printf("]\n");
}

void print_heap(node *heap, int last_one)
{
  for(int i=0;i<=last_one;i++)
    print_node(heap[i]);
  printf("\n");
}

inline long unsigned int gcd (long unsigned int a, long unsigned int b)
/* Euclid algorithm gcd */
{
	unsigned int c;
	while(a!=0)
	{
		c=a;
		a=b%a;
		b=c;
	};
	return(b);
};

inline long unsigned int gcd (long unsigned int a, long unsigned int b,
			      long unsigned int c)
{
  return(gcd(gcd(a,b),c));
}

inline long unsigned int gcd (long unsigned int a, long unsigned int b,
			      long unsigned int c, long unsigned int d)
{
  return(gcd(gcd(a,b),gcd(c,d)));
}


// returns 0 if top <left and top < right
// returns 1 if top > left < right
// returns -1 if top > right < left
inline int comp_nodes(node top, node left, node right)
{
  if(left.n<=right.n) // left <= right
    {
      if(top.n<=left.n) // top <= left <= right
	return(0);
      else // top > left and top <= right
	return(1); // so swap left
    }
  // right < left
  if(top.n<=right.n) // top <= right < left
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
	  "movapd %%XMM1,%2\n\t"					\
	  : "=m" (x), "=m" (y), "=m" (x.n), "=m" (y.n) \
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
	      if(heap[left_ptr].n<heap[this_ptr].n) // left<this
		swap_node(heap[this_ptr],heap[left_ptr]);
	      return; // tree is balanced
	    }
	}
      else // this node is at bottom
	return;
    }
}

#define pop_node(heap,heap_end){heap[0]=heap[heap_end];\
balance_heap(heap,heap_end-1);}

bigint START,END;

long int heap1_end;

node *make_heap1(bigint start, bigint end)
{
  double xlowd=(start-1)/4.0;
  long unsigned int xlow=ceil(sqrt(xlowd));
  double xhid=(end-1)>>2;
  long unsigned int xhi=floor(sqrt(xhid));
  heap1_end=xhi-xlow;
  printf("xlow=%ld xhi=%ld,len=%ld\n",xlow,xhi,heap1_end);
  long unsigned int num_nodes=heap1_end+1;
  node *nodes=(node *) memalign(16,sizeof(node)*num_nodes);
  if(!nodes)
    {
      printf("Failed to allocate memory for nodes1. Exiting.\n");
      exit(0);
    }
  long unsigned int i=0,x=xlow;
  for(;x<=xhi;x++,i++)
    {
      nodes[i].x=x;
      nodes[i].y=1;
      nodes[i].n=((x*x)<<2)+1;
    }
  balance_heap(nodes,heap1_end);
}

inline bigint next_n1 (node *nodes) // add two to y
{
  if(heap1_end<0)
    return(0);
  bigint n=nodes[0].n;
  nodes[0].n+=(nodes[0].y+1)<<2;
  if(nodes[0].n>END)
    {
      nodes[0]=nodes[heap1_end];
      heap1_end--;
      balance_heap(nodes,heap1_end);
      return(n);
    }
  nodes[0].y+=2;
  balance_heap(nodes,heap1_end);
  return(n);
}

int main(int argc, char **argv)
{
  node *heap1;
  bigint n,next_n;
  bigint next_square,next_square2;
  bool primep;
  long unsigned int count=0;

  START=10;
  END=START+1000000;
  heap1=make_heap1(START,END);


  /*
  n=heap1[0].n;
  while(n!=0)
    {
      printf("n=");print_bigint(n);printf("\n");
      n=next_n1(heap1);
    }
  exit(0);
  */


  print_heap(heap1,heap1_end);
  next_square=(heap1[0].x<<1)+1;
  next_square2=next_square*next_square;
  //printf("next_square is now ");print_bigint(next_square2);printf("\n");
  //printf("it is the square of ");print_bigint(next_square);printf("\n");
  n=next_n1(heap1);
  //print_heap(heap1,heap1_end);
  while(n!=0)
    {
      //printf("n=");print_bigint(n);printf("\n");
      while(n>next_square2)
	{
	  next_square2+=(next_square<<1)+1;
	  next_square++;
	  //printf("next_square is now ");print_bigint(next_square2);printf("\n");
	  //printf("it is the square of ");print_bigint(next_square);printf("\n");
	}
      if(n==next_square2)
	{
	  //printf("skipping a square\n");

	  n=next_n1(heap1);
	  while(n==next_square2)
	    n=next_n1(heap1);	  
	  continue;
	}
      if((n%5)==0)
	{
	  //printf("skipping 0 mod 5\n");
	  n=next_n1(heap1);
	  continue;
	}
      if(((n%12)!=1)&&((n%12)!=5))
	{
	  //printf("skipping not mod 1 or 5\n");
	  n=next_n1(heap1);
	  continue;
	}
      primep=true;
      next_n=next_n1(heap1);
      while(next_n==n)
	{
	  primep=!primep;
	  next_n=next_n1(heap1);
	}
      
      if(primep)
	{
	  //print_bigint(n);
	  //printf(" is prime.\n");
	  count++;
	}
      
      n=next_n;
    }
  printf("%ld prime found\n",count);
  return(0);
}
