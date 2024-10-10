#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_STR (10000000)

typedef struct{
 double k;
 double l;
 double k_l;} k_l_t;

void print_k_l(k_l_t k_l)
{
  printf("(k,l)=(%16.14f,%16.14f) : k+l-1/2 = %16.14f.\n",k_l.k,k_l.l,k_l.k_l);
}

k_l_t B(k_l_t inp)
{
  k_l_t op;
  op.k=inp.l-0.5;
  op.l=inp.k+0.5;
  op.k_l=inp.k_l;
  return op;
}

k_l_t A(k_l_t inp)
{
  k_l_t op;
  double num=2.0*(1.0+inp.k);
  op.k=inp.k/num;
  op.l=0.5+inp.l/num;
  op.k_l=op.k+op.l-0.5;
  return op;
}

k_l_t A_inv(k_l_t inp)
{
  k_l_t op;
  op.k=2*inp.k/(2.0*inp.k-1.0);
  op.l=(inp.l-0.5)*2.0*(1.0+op.k);
  op.k_l=op.k+op.l-0.5;
}
/*
k_l_t AB(k_l_t inp)
{
  return A(B(inp));
}

k_l_t ABA(k_l_t inp)
{
  return AB(A(inp));
}

k_l_t ABAA(k_l_t inp)
{
  return ABA(A(inp));
}
*/

int main(int argc, char **argv)
{
  srand(time(0));
  k_l_t init;
  init.k=0.0;
  init.l=1.0;
  init.k_l=init.k+init.l-0.5;

  double t0=atof(argv[1]);
  double t1=atof(argv[2]);
  double del=atof(argv[3]);
  int64_t steps=atol(argv[4]);

  k_l_t old_k_l,new_k_l;
  char *string;
  string=(char*)malloc(MAX_STR);
  string[0]=0;
  int64_t string_ptr=0;
  for(double t=t0;t>t1;t*=del)
    {
      for(int64_t step=0;step<steps;step++)
	{
	  if(string_ptr>=MAX_STR)
	    {
	      printf("string length exceeded.\n");
	      return 0;
	    }
	  int r=rand()%4;
	  switch(r)
	    {
	    case 1: new_k_l=A(old_k_l);break;
	    case 2:
	      if((string_ptr>0)&&(string[string_ptr-1]=='A'))
		new_k_l=A_inv(old_k_l);
	      break;
	    case 3: new_k_l=B(old_k_l);break;
	    default: new_k_l=B(old_k_l);
	    }
	  
	  double benefit=old_k_l.k_l-new_k_l.k_l; // +ve if improved
	  if((benefit>=0.0)||
	     (expl(benefit/t)>(double)rand() / RAND_MAX))
	    {
	      switch(r)
		{
		case 1: // do an A
		  if(((string_ptr>0)&&(string[string_ptr-1]=='B'))||
		     ((string_ptr>1)&&(string[string_ptr-2]=='B'))||
		     ((string_ptr>2)&&(string[string_ptr-3]=='B')))
		    {
		      string[string_ptr++]='A';
		      old_k_l=new_k_l;
		    }
		  break;
		case 2: // undo an A
		  if((string_ptr>0)&&(string[string_ptr-1]=='A'))
		    {
		      string_ptr--;
		      old_k_l=new_k_l;
		    }
		  break;
		default:
		  string[string_ptr++]='B';
		  old_k_l=new_k_l;
		}
	      if(string_ptr>=2)
		if((string[string_ptr-1]=='B')&&(string[string_ptr-2]=='B'))
		  string_ptr-=2;
	      string[string_ptr]=0;
	    }
	}
    }
  printf("New string = %s\n",string);
  print_k_l(old_k_l);
  
  return 0;
}
  
