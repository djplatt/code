#include "stdio.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <gmp.h>

#define H (2000000/(13*17)) // runs to height<=H*modulus ??
#define N_PRIMES (2)
unsigned long int primes[N_PRIMES]={13,17};
//unsigned long int primes[1];
#define modulus ((long unsigned int) 13*17)

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
inline void p(bigint res, long unsigned int x, long unsigned int a, long unsigned int offset)
{
  mpz_set_ui(res,(x-1)*modulus+offset);
  mpz_mul(res,res,res); // x^2
  mpz_mul(res,res,res); // x^4
  mpz_mul_ui(res,res,a);// ax^4
}

bool offsets_ok(unsigned long int *offsets, unsigned long int *As, unsigned long int prime)
{
  long int x=offsets[0]%prime;
  long int y=offsets[1]%prime;
  long int u=offsets[2]%prime;
  long int v=offsets[3]%prime;
  if((x==0)&&(y==0)&&(u==0)&&(v==0))
    return(false);
  //printf("Checking %lu %lu %lu %lu (mod %lu)\n",x,y,u,v,prime);
  x*=x;x*=x;y*=y;y*=y;u*=u;u*=u;v*=v;v*=v;
  return(((As[0]*x+As[1]*y-As[2]*u-As[3]*v)%prime)==0);
}

inline bool offsets_ok(unsigned long int *offsets, unsigned long int *As)
{
  for(int i=0;i<N_PRIMES;i++)
    if(!(offsets_ok(offsets,As,primes[i])))
      return(false);
  return(true);
}


void check_eqn(long unsigned int A1, long unsigned int A2, 
	       long unsigned int A3, long unsigned int A4, 
	       mpz_t temp, mpz_t temp1, node *lnodes, node *rnodes,
	       mpz_t *ps, mpz_t *qs,
	       long unsigned int *offsets)
{
  node left_node,right_node;
  int cmp;

  // set up left heap
  p(temp1,1,A1,offsets[0]);
  for(long unsigned int i=0;i<H;i++)
    {
      lnodes[i].a=1;
      lnodes[i].b=i+1;
      mpz_set(lnodes[i].pq,temp1);
      p(temp,lnodes[i].b,A2,offsets[1]);
      mpz_add(lnodes[i].pq,lnodes[i].pq,temp);
    }

  std::vector<node> v1(lnodes,lnodes+H);
  std::vector<node>::iterator it;

  std::make_heap(v1.begin(),v1.end(),node_comp);

  // set up right heap
  p(temp1,1,A3,offsets[2]);
  for(long unsigned int i=0;i<H;i++)
    {
      rnodes[i].a=1;
      rnodes[i].b=i+1;
      mpz_set(rnodes[i].pq,temp1);
      p(temp,rnodes[i].b,A4,offsets[3]);
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


int main()
{
  node lnodes[H],rnodes[H],left_node,right_node;
  mpz_t temp,temp1,ps[H],qs[H];
  int cmp;
  unsigned long int As[4]={1,2,1,4}; // this defines the equation
  unsigned long int offsets[4]={1,1,1,1};

  printf("Looking for solutions to %lux^4+%luy^4=%luu^4+%luv^4 to height %lu\n",
	 As[0],As[1],As[2],As[3],H*modulus);
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
  while(true)
    {
      if(offsets_ok(offsets,As))
	{
	  printf("Running with offsets= %lu %lu %lu %lu\n",offsets[0],offsets[1],offsets[2],offsets[3]);
	  for(int i=0;i<H;i++)
	    {
	      p(ps[i],i+1,As[0],offsets[0]);
	      p(qs[i],i+1,As[2],offsets[2]);
	    }
	  for(int i=0;i<H-1;i++)
	    {
	      mpz_sub(ps[i],ps[i+1],ps[i]);
	      mpz_sub(qs[i],qs[i+1],qs[i]);
	    }
	  check_eqn(As[0],As[1],As[2],As[3],temp,temp1,lnodes,rnodes,ps,qs,offsets);
	}
      if(offsets[3]!=modulus)
	offsets[3]++;
      else
	{
	  offsets[3]=1;
	  if(offsets[2]!=modulus)
	    offsets[2]++;
	  else
	    {
	      offsets[2]=1;
	      if(offsets[1]!=modulus)
		offsets[1]++;
	      else
		{
		  offsets[1]=1;
		  printf("First offset increased to %lu\n",offsets[0]+1);
		  if(offsets[0]!=modulus)
		    offsets[0]++;
		  else
		    break;
		}
	    }
	}
    }
}
