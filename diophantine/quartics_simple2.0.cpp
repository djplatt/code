#include "inttypes.h"
#include "stdio.h"
#include <stdlib.h>
#include <malloc.h>
#include <iostream>
#include <assert.h>

unsigned long int H; // runs to height<=H

// 128 bit unsigned built in to GCC (efficiency?)
typedef __uint128_t bigint;

// a node in our heap
typedef struct{
  //unsigned int node_no;
  long unsigned int a;
  long unsigned int b;
  bigint pq;
} node;


void print_bigint(bigint i)
{
  printf("%ju",(uintmax_t) i);
}

void print_node(node n)
{
  printf("[%lu,%lu,",n.a,n.b);
  print_bigint(n.pq);
  printf("]\n");
}

void print_heap(node *heap, int last_one)
{
  for(int i=0;i<=last_one;i++)
    print_node(heap[i]);
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

void print_solution(node n1, node n2)
{
  printf("solution found\n");
  print_node(n1);
  print_node(n2);
  //exit(0);
}

// res<-a*x^4
inline bigint p(long unsigned int x, long unsigned int a)
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
int comp_nodes(node top, node left, node right)
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

#define swap_node(x,y) {\
long unsigned int tmp; bigint tmp1; tmp1=x.pq;x.pq=y.pq;y.pq=tmp1;\
tmp=x.a;x.a=y.a;y.a=tmp;\
tmp=x.b;x.b=y.b;y.b=tmp;}


//#define swap_node(x,y) {node tmp; tmp=x;x=y;y=tmp;}

// called with a heap with the 0th element possibly in
// wrong place. The last element is at heap[heap_end]
inline void balance_heap (node *heap, int heap_end)
{
  int this_ptr=0,left_ptr,right_ptr,cmp;
  while(true)
    {
      left_ptr=(this_ptr<<1)+1;
      right_ptr=(this_ptr+1)<<1;
      if(left_ptr<=heap_end) // it has a left child
	{
	  if(right_ptr<=heap_end) // it has a right child
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
      else
	return;
    }
}

#define pop_node(heap,heap_end){heap[0]=heap[heap_end];\
balance_heap(heap,heap_end-1);}

void check_eqn(long unsigned int A1, long unsigned int A2, 
	       long unsigned int A3, long unsigned int A4, 
	       node *lnodes, node *rnodes,
	       bigint *ps, bigint *qs)
{
  int cmp;
  bigint temp,temp1;

  // set up left heap
  temp1=p(1,A1);
  for(long unsigned int i=0;i<H;i++)
    {
      lnodes[i].a=1;
      lnodes[i].b=i+1;
      lnodes[i].pq=temp1+p(lnodes[i].b,A2);
    }
  //printf("Left heap now\n");print_heap(lnodes,H-1);

  // set up right heap
  temp1=p(1,A3);
  for(long unsigned int i=0;i<H;i++)
    {
      rnodes[i].a=1;
      rnodes[i].b=i+1;
      rnodes[i].pq=temp1+p(rnodes[i].b,A4);
   }
  //printf("Right heap now\n");print_heap(rnodes,H-1);

  int left_end=H-1;
  int right_end=H-1;

  while(true)
    {
      //printf("Iterating with left=");print_node(lnodes[0]);
      //printf("          and right=");print_node(rnodes[0]);
      if(lnodes[0].pq<rnodes[0].pq) // left<right
	{
	  if(lnodes[0].a<H) // more left nodes to add
	    {
	      lnodes[0].pq+=ps[lnodes[0].a-1];
	      lnodes[0].a++;
	      //printf("L Pushing ");print_node(lnodes[0]);
	      balance_heap(lnodes,H-1);
	      //print_heap(lnodes,H-1);
	    }
	  else // no more left nodes, so just pop
	    pop_node(lnodes,left_end--);
	  if(left_end<0) // heap empty
	    return;
	  continue;
	}
      if(lnodes[0].pq>rnodes[0].pq) // right<left
	{
	  if(rnodes[0].a<H)
	    {
	      rnodes[0].pq+=qs[rnodes[0].a-1];
	      rnodes[0].a++;
	      //printf("R Pushing ");print_node(rnodes[0]);
	      balance_heap(rnodes,H-1);
	      //print_heap(rnodes,H-1);
	    }
	  else
	    pop_node(rnodes,right_end--);
	  if(right_end<0)
	    return;
	  continue;
	}
      // right=left
      if(gcd(lnodes[0].a,lnodes[0].b,rnodes[0].a,rnodes[0].b)==1)
	print_solution(lnodes[0],rnodes[0]);
      // remove top of left heap
      if(lnodes[0].a<H) // more left nodes to add
	{
	  lnodes[0].pq+=ps[lnodes[0].a-1];
	  lnodes[0].a++;
	  //printf("L Pushing ");print_node(lnodes[0]);
	  balance_heap(lnodes,H-1);
	  //print_heap(lnodes,H-1);
	}
      else // no more left nodes, so just pop
	pop_node(lnodes,left_end--);
      if(left_end<0) // heap empty
	return;
    }
}

int main(int argc, char **argv)
{
  node *lnodes,*rnodes,left_node,right_node;
  bigint *ps,*qs;
  int cmp;
  unsigned long int As[4]; // this defines the equation

  if(argc!=6)
    {
      printf("Incorrect command line, #args=%d\n",argc);
      exit(0);
    }

  H=atoi(argv[1]);

  for(int i=0;i<4;i++)
    As[i]=atoi(argv[i+2]);

  assert(lnodes=(node *) malloc(H*sizeof(node)));
  assert(rnodes=(node *) malloc(H*sizeof(node)));
  assert(ps=(bigint *) malloc(H*sizeof(bigint)));
  assert(qs=(bigint *) malloc(H*sizeof(bigint)));

  printf("Looking for solutions to %lux^4+%luy^4=%luu^4+%luv^4 to height %lu\n",
	 As[0],As[1],As[2],As[3],H);
  for(int i=0;i<H;i++)
    {
      ps[i]=p(i+1,As[0]);
      qs[i]=p(i+1,As[2]);
    }
  for(int i=0;i<H-1;i++)
    {
      ps[i]=ps[i+1]-ps[i];
      qs[i]=qs[i+1]-qs[i];
    }
  check_eqn(As[0],As[1],As[2],As[3],lnodes,rnodes,ps,qs);
}
