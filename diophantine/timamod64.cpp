// version 2.1 use SSE instructions to move nodes
// stop after first solution
//
#include "inttypes.h"
#include "stdio.h"
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <assert.h>

unsigned long int H; // runs to height<=H

//#define BIT128

#ifdef BIT128
// 128 bit unsigned built in to GCC (efficiency?)
typedef __int128_t bigint;
#else
typedef int64_t bigint; 
#endif

// a node in our heap
typedef struct{
  //int node_no;
  __attribute__ ((aligned(16)))
  long int x;
  long int y;
  bigint p;
#ifndef BIT128 // pack out to 128 bits
  bigint dummy;
#endif
} node;

#ifdef BIT128
void print_bigint(bigint i)
{
  printf("%ju",(uintmax_t) i);
}
#else
void print_bigint(bigint i)
{
  printf("%lu",i);
}
#endif

void print_node(node n)
{
  printf("[%lu,%lu,",n.x,n.y);
  print_bigint(n.p);
  printf("]\n");
}

void print_heap(node *heap, int last_one)
{
  for(int i=0;i<=last_one;i++)
    print_node(heap[i]);
}

// returns 0 if top <left and top < right
// returns 1 if top > left < right
// returns -1 if top > right < left
inline int comp_nodes(node top, node left, node right)
{
  if(left.p<=right.p) // left <= right
    {
      if(top.p<=left.p) // top <= left <= right
	return(0);
      else // top > left and top <= right
	return(1); // so swap left
    }
  // right < left
  if(top.p<=right.p) // top <= right < left
    return(0);
  else
    return(-1);
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

void print_solution(node n1, node n2)
{
  printf("solution found\n");
  print_node(n1);
  print_node(n2);
  //printf("Terminating search.\n");
  //exit(0);
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
	  : "=m" (x), "=m" (y), "=m" (x.p), "=m" (y.p)	\
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
	      if(heap[left_ptr].p<heap[this_ptr].p) // left<this
		swap_node(heap[this_ptr],heap[left_ptr]);
	      return; // tree is balanced
	    }
	}
      else // this node is at bottom
	return;
    }
}


void check_eqn(node *heap, bigint *xcubes, bigint *ycubes, bigint *zcubes, long int H)
{

  long int heap_end=H-1;

  for(long int i=0;i<H;i++)
    {
      heap[i].x=i;
      heap[i].y=0;
      heap[i].p=xcubes[i]+ycubes[0];
    }
  balance_heap(heap,heap_end);
  //print_heap(heap,heap_end);

  long int z=0;
  bigint z3=zcubes[0];

  while(true)
    {
      while(z3<heap[0].p)
	{
	  if(z==H-1)
	    {
	      printf("Z reached H limit.\n");
	      exit(0);
	    }
	  z++;
	  z3=zcubes[z];
	}
      //printf("z3=");print_bigint(z3);printf(" x=%lu y=%lu p=",heap[0].x,heap[0].y);print_bigint(heap[0].p);printf("\n");
      if(z3==heap[0].p)
	{
	  printf("Solution found with x=%lu, y=%lu, z=%lu\n",heap[0].x,heap[0].y,z); 
	  exit(0);
	}
      if(heap[0].y==heap[0].x)
	{
	  if(heap_end==0)
	    {
	      printf("Heap empty.\n");
	      exit(0);
	    }
	  //printf("popping ");print_node(heap[0]);
	  heap[0]=heap[heap_end--];
	  balance_heap(heap,heap_end);
	  //printf("new head ");print_node(heap[0]);
	}
      else
	{
	  heap[0].y++;
	  heap[0].p=xcubes[heap[0].x]+ycubes[heap[0].y];
	  balance_heap(heap,heap_end);
	  //printf("new head ");print_node(heap[0]);
	}
    }
}

int main(int argc, char **argv)
{
  node *nodes;
  bigint *xcubes,*ycubes,*zcubes;
  long int W,H,N;

  if(argc!=7)
    {
      printf("Incorrect command line, W H N xmod ymod zmod #args=%d\n",argc);
      exit(0);
    }


  H=atol(argv[2]);
  W=atol(argv[1]);
  N=atol(argv[3]);


  printf("Searching for solutions to x^3+y^3=z^3+%lu in N^3 to height %lu\n",W,H*N);

  assert(nodes=(node *) memalign(16,H*sizeof(node)));
  assert(xcubes=(bigint *) malloc(H*sizeof(bigint)));
  assert(ycubes=(bigint *) malloc(H*sizeof(bigint)));
  assert(zcubes=(bigint *) malloc(H*sizeof(bigint)));
  long int x=atol(argv[4]),y=atol(argv[5]),z=atol(argv[6]);
  printf("Working modulo %lu with x=%lu y=%lu and z=%lu.\n",N,x,y,z);

  for(long int i=0;i<H;i++)
    {
      xcubes[i]=x;xcubes[i]*=x;xcubes[i]*=x;x+=N;
      ycubes[i]=y;ycubes[i]*=y;ycubes[i]*=y;y+=N;
      zcubes[i]=z;zcubes[i]*=z;zcubes[i]*=z;zcubes[i]+=W;z+=N;
    }
  check_eqn(nodes,xcubes,ycubes,zcubes,H);
}
