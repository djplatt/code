#include "stdio.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <gmp.h>

#define H (100000) // runs to height 16*H

// using GMP's arbitary length integer type MPZ
typedef mpz_t bigint;

// a node in our heap
typedef struct{
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
	       mpz_t temp, node *lnodes, node *rnodes, long unsigned int cong)
{
  node left_node,right_node;
  int cmp;
  long unsigned int congs[4];
  congs[0]=(cong&1);
  congs[1]=(cong&2)/2;
  congs[2]=(cong&4)/4;
  congs[3]=(cong&8)/8;

  // don't consider solutions containing zeros
  for(int i=0;i<4;i++)
    if(congs[i]==0)
      congs[i]=16;

  // set up left heap
  for(long unsigned int i=0;i<H;i++)
    {
      lnodes[i].a=congs[0];
      lnodes[i].b=i*16+congs[1];
      p(lnodes[i].pq,lnodes[i].a,A1);
      p(temp,lnodes[i].b,A2);
      mpz_add(lnodes[i].pq,lnodes[i].pq,temp);
    }

  std::vector<node> v1(lnodes,lnodes+H);
  std::vector<node>::iterator it;

  std::make_heap(v1.begin(),v1.end(),node_comp);

  // set up right heap
  for(long unsigned int i=0;i<H;i++)
    {
      rnodes[i].a=congs[2];
      rnodes[i].b=(i<<4)+congs[3];
      p(rnodes[i].pq,rnodes[i].a,A3);
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
	  if(left_node.a<(H<<4))
	    {
	      p(temp,left_node.a+16,A1);
	      mpz_add(left_node.pq,left_node.pq,temp);
	      p(temp,left_node.a,A1);
	      mpz_sub(left_node.pq,left_node.pq,temp);
	      //left_node.pq=p(left_node.a+1)-p(left_node.a)+left_node.pq;
	      left_node.a+=16;
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
	  if(right_node.a<(H<<4))
	    {
	      p(temp,right_node.a+16,A3);
	      mpz_add(right_node.pq,right_node.pq,temp);
	      p(temp,right_node.a,A3);
	      mpz_sub(right_node.pq,right_node.pq,temp);
	      //right_node.pq=r(right_node.a+1)-r(right_node.a)+right_node.pq;
	      right_node.a+=16;
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
      if(left_node.a<(H<<4))
	{
	  p(temp,left_node.a+16,A1);
	  mpz_add(left_node.pq,left_node.pq,temp);
	  p(temp,left_node.a,A1);
	  mpz_sub(left_node.pq,left_node.pq,temp);
	  //left_node.pq=p(left_node.a+1)-p(left_node.a)+left_node.pq;
	  left_node.a+=16;
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
      if(right_node.a<(H<<4))
	{
	  p(temp,right_node.a+16,A3);
	  mpz_add(right_node.pq,right_node.pq,temp);
	  p(temp,right_node.a,A3);
	  mpz_sub(right_node.pq,right_node.pq,temp);
	  //right_node.pq=r(right_node.a+1)-r(right_node.a)+right_node.pq;
	  right_node.a+=16;
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

// confirm that a set of congruences mod 16 work
inline int congs_ok(long unsigned int *As, long unsigned int congs)
{
  int res=As[0]*(congs&1)+As[1]*((congs&2)>>1)-As[2]*((congs&4)>>2)
    -As[3]*((congs&8)>>3);
  //printf("congs=%X -> res=%d\n",congs,res);
  return(!(res&15));
}


int main()
{
  node lnodes[H],rnodes[H],left_node,right_node;
  mpz_t temp;
  int cmp;
  long unsigned int congs,As[4]={1,2,1,4};
  mpz_init(temp);
  for(int i=0;i<H;i++)
    {
      mpz_init(lnodes[i].pq);
      mpz_init(rnodes[i].pq);
    }
  for(long unsigned int congs=0;congs<15;congs++)
    if(congs_ok(As,congs))
      {
	printf("Congruence %lu is ok \n",congs);
	check_eqn(As[0],As[1],As[2],As[3],temp,lnodes,rnodes,congs);
      }
}
