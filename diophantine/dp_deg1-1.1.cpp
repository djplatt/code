//
// Solutions to Del Pezzo surfaces of degree 1
//
// w^2=z^3+ax^6+bx^6
//
#include "inttypes.h"
#include "stdio.h"
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
//#include <assert.h>
#include <math.h>

#define my_assert(f,str) if(!(f)) {printf("Problem with %s. Exiting.\n",str);exit(0);}
long int H; // runs to height<=H

// 128 bit unsigned built in to GCC (efficiency?)
typedef __int128_t bigint;

// a node in our heap
typedef struct{
  //unsigned int node_no;
  __attribute__ ((aligned(16)))
  long int a;
  long int b;
  bigint pq;
} node;


void print_bigint(bigint i)
{
  printf("%jd",(intmax_t) i);
}

void print_node(node *n)
{
  printf("[%ld,%ld,",n->a+1,n->b+1);
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

void print_solution(node *n1, node *n2)
{
  if(gcd(n1[0].a+1,n1[0].b+1,n2[0].a+1,n2[0].b+1)==1)
    {
      printf("solution %ld found\n",++solution_count);
      print_node(n1);
      print_node(n2);
      printf("Height= %ld.\n",max(max(n1->a+1,n1->b+1),max(n2->a+1,n2->b+1)));
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
      lnodes[i].b=i;
      lnodes[i].pq=squares[0]+cubes[i];
      rnodes[i].a=0;
      rnodes[i].b=i;
      rnodes[i].pq=as[0]+bs[i];
    }

  int left_end=H-1;
  int right_end=H-1;

  /*
  bigint max_left=H;
  max_left*=H;
  max_left+=max_left*H;
  while(rnodes[right_end].pq>max_left) right_end--;
  */
  while(true)
    {
      //printf("Iterating with left=");print_node(lnodes[0]);
      //printf("          and right=");print_node(rnodes[0]);
      if(lnodes[0].pq<rnodes[0].pq) // left<right
	{
	  if(lnodes[0].a<H-1) // more left nodes to add
	    {
	      lnodes[0].a++;
	      lnodes[0].pq=squares[lnodes[0].a]+cubes[lnodes[0].b];
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
	  if(rnodes[0].a<H-1)
	    {
	      rnodes[0].a++;
	      rnodes[0].pq=as[rnodes[0].a]+bs[rnodes[0].b];
	      //printf("R Pushing ");print_node(rnodes[0]);
	      balance_heap(rnodes,right_end);
	      /*
	      if(max_left<rnodes[right_end].pq)
		{
		  printf("Getting rid of ");
		  print_node(rnodes+right_end);
		  right_end--;
		}
	      */
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

      if(lnodes[0].a<H-1) // more left nodes to add
	{
	  lnodes[0].a++;
	  lnodes[0].pq=squares[lnodes[0].a]+cubes[lnodes[0].b];
	  //printf("L Pushing ");print_node(lnodes[0]);
	  balance_heap(lnodes,left_end);
	  //print_heap(lnodes,H-1);
	}
      else // no more room on this lnode, so pop it and move RHS on
	{
	  pop_node(lnodes,left_end);
	  if(left_end<0) // heap empty
	    {
	      //print_node(rnodes[0]);
	      printf("Heap max reached (left).\n");
	      return;
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

  for(int i=0;i<2;i++)
    As[i]=atoi(argv[i+2]);

  printf("RHS max H set to %ld.\n",H);

  my_assert(lnodes=(node *) memalign(16,H*sizeof(node)),"allocating memory for lnodes");
  my_assert(rnodes=(node *) memalign(16,H*sizeof(node)),"allocating memory for rnodes");
  my_assert(squares=(bigint *) malloc(H*sizeof(bigint)),"allocating memory for squares");
  my_assert(cubes=(bigint *) malloc(H*sizeof(bigint)),"allocating memory for cubes");
  my_assert(as=(bigint *) malloc(H*sizeof(bigint)),"allocating memory for as");
  my_assert(bs=(bigint *) malloc(H*sizeof(bigint)),"allocating memory for bs");


  printf("Looking for solutions to w^2+z^3=%ldx^6+%ldy^6 in (Z>0)^4 to height %ld\n",
	 As[0],As[1],H);
  for(int i=0;i<H;i++)
    {
      squares[i]=i+1;squares[i]*=i+1;
      cubes[i]=squares[i]*(i+1);
      bigint sixth=cubes[i]*cubes[i];
      as[i]=As[0]*sixth;
      bs[i]=As[1]*sixth;
    }
  check_eqn(lnodes,rnodes,squares,cubes,as,bs);
}
