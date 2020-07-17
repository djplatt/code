#include "stdio.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <gmp.h>
#include <assert.h>

#define H ((unsigned long int) 100000) // runs to height<=H

// using GMP's arbitary length integer type MPZ
typedef mpz_t bigint;

// a node in our heap
typedef struct{
  unsigned int node_no;
  long unsigned int a;
  long unsigned int b;
  bigint pq;
} node;

// 
inline int node_comp(node n1, node n2)
{
  return(mpz_cmp(n1.pq,n2.pq)>0);
}

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


void check_eqn(long unsigned int A1, long unsigned int A2, 
	       long unsigned int A3, long unsigned int A4, 
	       mpz_t temp, mpz_t temp1, node *lnodes, node *rnodes,
	       mpz_t *ps, mpz_t *qs)
{
  node left_node,right_node;
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

  std::vector<node> v1(lnodes,lnodes+H);
  std::vector<node>::iterator it;

  std::make_heap(v1.begin(),v1.end(),node_comp);

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

  std::vector<node> v2(rnodes,rnodes+H);
  std::make_heap(v2.begin(),v2.end(),node_comp);

  // get two lowest entries
  std::pop_heap(v1.begin(),v1.end(),node_comp);v1.pop_back();
  left_node=v1[v1.size()];
  std::pop_heap(v2.begin(),v2.end(),node_comp);v2.pop_back();
  right_node=v2[v2.size()];

  while(true)
    {
      //printf("Iterating with left=");print_node(left_node);
      //printf("          and right=");print_node(right_node);
      cmp=mpz_cmp(left_node.pq,right_node.pq);
      if(cmp<0) 
	{
	  if(left_node.a<H)
	    {
	      mpz_add(left_node.pq,left_node.pq,ps[left_node.a-1]);
	      left_node.a++;
	      v1.push_back(left_node);
	      std::push_heap(v1.begin(),v1.end(),node_comp);
	    }
	  if(v1.size())
	    {
	      std::pop_heap(v1.begin(),v1.end(),node_comp);v1.pop_back();
	      left_node=v1[v1.size()];
	      continue;
	    }
	  else
	    return;
	}
      if(cmp>0)
	{
	  if(right_node.a<H)
	    {
	      mpz_add(right_node.pq,right_node.pq,qs[right_node.a-1]);
	      right_node.a++;
	      v2.push_back(right_node);
	      std::push_heap(v2.begin(),v2.end(),node_comp);
	    }
	  if(v2.size())
	    {
	      std::pop_heap(v2.begin(),v2.end(),node_comp);v2.pop_back();
	      right_node=v2[v2.size()];
	      continue;
	    }
	  else
	    return;
	}
      print_solution(left_node,right_node);
      if(left_node.a<H)
	{
	  mpz_add(left_node.pq,left_node.pq,ps[left_node.a-1]);
	  left_node.a++;
	  v1.push_back(left_node);
	  std::push_heap(v1.begin(),v1.end(),node_comp);
	}
      if(v1.size())
	{
	  std::pop_heap(v1.begin(),v1.end(),node_comp);v1.pop_back();
	  left_node=v1[v1.size()];
	  continue;
	}
      else
	return;
      if(right_node.a<H)
	{
	  mpz_add(right_node.pq,right_node.pq,qs[right_node.a-1]);
	  right_node.a++;
	  v2.push_back(right_node);
	  std::push_heap(v2.begin(),v2.end(),node_comp);
	}
      if(v2.size())
	{
	  std::pop_heap(v2.begin(),v2.end(),node_comp);v2.pop_back();
	  right_node=v2[v2.size()];
	  continue;
	}
      else
	return;
    }
}


int main(int argc, char **argv)
{
  node *lnodes,*rnodes,left_node,right_node;
  mpz_t temp,temp1,*ps,*qs;
  int cmp;
  unsigned long int As[4]; // this defines the equation

  if(argc!=5)
    exit(0);

  for(int i=0;i<4;i++)
    As[i]=atoi(argv[i+1]);

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
      lnodes[i].node_no=i;
      rnodes[i].node_no=i;
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
