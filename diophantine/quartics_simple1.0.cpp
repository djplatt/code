#include "stdio.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <gmp.h>
#include <assert.h>

unsigned long int H; // runs to height<=H

// using GMP's arbitary length integer type MPZ
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
  mpz_out_str(NULL,10,i);
}

void print_node(node n)
{
  printf("[%lu,%lu,",n.a,n.b);
  print_bigint(n.pq);
  printf("]\n");
}

void print_solution(node n1, node n2)
{
  printf("solution found\n");
  print_node(n1);
  print_node(n2);
  //exit(0);
}

// res<-a*x^4
inline void p(bigint res, long unsigned int x, long unsigned int a)
{
  mpz_set_ui(res,x);
  mpz_mul(res,res,res); // x^2
  mpz_mul(res,res,res); // x^4
  mpz_mul_ui(res,res,a);// ax^4
}

#define mpz_le(a,b) (mpz_cmp(a,b)<=0)

// returns 0 if top <left and top < right
// returns 1 if top > left < right
// returns -1 if top > right < left
int comp_nodes(node top, node left, node right)
{
  if(mpz_le(left.pq,right.pq)) // left <= right
    {
      if(mpz_le(top.pq,left.pq)) // top <= left <= right
	return(0);
      else // top > left and top <= right
	return(1); // so swap left
    }
  // right < left
  if(mpz_le(top.pq,right.pq)) // top <= right < left
    return(0);
  else
    return(-1);
} 

#define swap_node(x,y) {int tmp; mpz_swap(x.pq,y.pq);tmp=x.a;x.a=y.a;y.a=tmp;\
    tmp=x.b;x.b=y.b;y.b=tmp;}

// called with a heap with the 0th element possibly in
// wrong place. The last element is at heap[heap_end]
inline void balance_heap (node *heap, int heap_end)
{
  int this_ptr=0,left_ptr,right_ptr,cmp;
  while(true)
    {
      left_ptr=this_ptr<<1+1;
      right_ptr=(this_ptr+1)<<1;
      if(left_ptr<=heap_end) // it has a left child
	{
	  if(right_ptr<=heap_end) // it has a right child
	    {
	      cmp=comp_nodes(heap[this_ptr],heap[left_ptr],heap[right_ptr]);
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
	      if(mpz_cmp(heap[left_ptr].pq,heap[this_ptr].pq)<0) // left<this
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
	       mpz_t temp, mpz_t temp1, node *lnodes, node *rnodes,
	       mpz_t *ps, mpz_t *qs)
{
  int cmp;

  // set up left heap
  p(temp1,1,A1);
  for(long unsigned int i=0;i<H;i++)
    {
      lnodes[i].a=1;
      lnodes[i].b=i+1;
      mpz_set(lnodes[i].pq,temp1);
      p(temp,lnodes[i].b,A2);
      mpz_add(lnodes[i].pq,lnodes[i].pq,temp);
    }

  // set up right heap
  p(temp1,1,A3);
  for(long unsigned int i=0;i<H;i++)
    {
      rnodes[i].a=1;
      rnodes[i].b=i+1;
      mpz_set(rnodes[i].pq,temp1);
      p(temp,rnodes[i].b,A4);
      mpz_add(rnodes[i].pq,rnodes[i].pq,temp);
   }

  int left_end=H-1;
  int right_end=H-1;

  while(true)
    {
      //printf("Iterating with left=");print_node(lnodes[0]);
      //printf("          and right=");print_node(rnodes[0]);
      cmp=mpz_cmp(lnodes[0].pq,rnodes[0].pq);
      if(cmp<0) // left<right
	{
	  if(lnodes[0].a<H) // more left nodes to add
	    {
	      mpz_add(lnodes[0].pq,lnodes[0].pq,ps[lnodes[0].a-1]);
	      lnodes[0].a++;
	      balance_heap(lnodes,H-1);
	    }
	  else // no more left nodes, so just pop
	    pop_node(lnodes,left_end--);
	  if(left_end<0) // heap empty
	    return;
	  continue;
	}
      if(cmp>0) // right<left
	{
	  if(rnodes[0].a<H)
	    {
	      mpz_add(rnodes[0].pq,rnodes[0].pq,qs[rnodes[0].a-1]);
	      rnodes[0].a++;
	      balance_heap(rnodes,H-1);
	    }
	  else
	    pop_node(rnodes,right_end--);
	  if(right_end<0)
	    return;
	  continue;
	}
      // right=left
      print_solution(lnodes[0],rnodes[0]);
      // now replace or pop left and right heaps
      if(lnodes[0].a<H)
	{
	  mpz_add(lnodes[0].pq,lnodes[0].pq,ps[lnodes[0].a-1]);
	  lnodes[0].a++;
	  balance_heap(lnodes,H-1);
	}
      else
	pop_node(lnodes,left_end--);
      if(left_end<0)
	return;
      if(rnodes[0].a<H)
	{
	  mpz_add(rnodes[0].pq,rnodes[0].pq,qs[rnodes[0].a-1]);
	  rnodes[0].a++;
	  balance_heap(rnodes,H-1);
	}
      else
	pop_node(rnodes,right_end--);
      if(right_end<0)
	return;
    }
}


int main(int argc, char **argv)
{
  node *lnodes,*rnodes,left_node,right_node;
  mpz_t temp,temp1,*ps,*qs;
  int cmp;
  unsigned long int As[4]; // this defines the equation

  if(argc!=6)
    exit(0);

  H=atoi(argv[1]);

  for(int i=0;i<4;i++)
    As[i]=atoi(argv[i+2]);

  assert(lnodes=(node *) malloc(H*sizeof(node)));
  assert(rnodes=(node *) malloc(H*sizeof(node)));
  assert(ps=(mpz_t *) malloc(H*sizeof(mpz_t)));
  assert(qs=(mpz_t *) malloc(H*sizeof(mpz_t)));

  printf("Looking for solutions to %lux^4+%luy^4=%luu^4+%luv^4 to height %lu\n",
	 As[0],As[1],As[2],As[3],H);
  mpz_init(temp);
  mpz_init(temp1);
  for(int i=0;i<H;i++)
    {
      //lnodes[i].node_no=i;
      //rnodes[i].node_no=i;
      mpz_init(lnodes[i].pq);
      mpz_init(rnodes[i].pq);
      mpz_init(ps[i]);
      mpz_init(qs[i]);
    }
  for(int i=0;i<H;i++)
    {
      p(ps[i],i+1,As[0]);
      p(qs[i],i+1,As[2]);
    }
  for(int i=0;i<H-1;i++)
    {
      mpz_sub(ps[i],ps[i+1],ps[i]);
      mpz_sub(qs[i],qs[i+1],qs[i]);
    }
  check_eqn(As[0],As[1],As[2],As[3],temp,temp1,lnodes,rnodes,ps,qs);
}
