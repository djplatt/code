#include <stdio.h>
#include <stdlib.h>
#include "int_double.h"
#define POS (0)
#define NEG (1)
#define TO_ZERO (2)
#define FROM_ZERO (3)
#define CROSS_ZERO (4)

#define my_print(str,val) {printf(str);printf(": ");print_int_double(val);printf("\n");}

double rand_double()
{
  return (double)rand()/RAND_MAX;
}

int_double rand_int_double(char type)
{
  int_double x;
  
  switch(type)
    {
    case POS:
      x.left=rand_double();
      x.right=x.left+1.0/RAND_MAX;
      return(x);
    case NEG:
      x.left=rand_double();
      x.right=x.left+1.0/RAND_MAX;
      return(-x);
    case TO_ZERO:
      x.left=-rand_double();
      x.right=0.0;
      return(x);
    case FROM_ZERO:
      x.left=0.0;
      x.right=rand_double();
      return(x);
    case CROSS_ZERO:
      x.left=-rand_double();
      x.right=rand_double();
      return(x);
    default:
      printf("Fatal error in rand_int_double. Exiting.\n");
      exit(0);
      return(x);
    }
}

// ((x*y)/x)/y)=1
bool testa (unsigned int n)
{
  int_double x,y,z;
  unsigned int i;
  char type1,type2;
  
  for(type1=POS;type1<=NEG;type1++)
    for(type2=POS;type2<=NEG;type2++)
      for(i=0;i<n;i++)
	{
	  x=rand_int_double(type1);
	  y=rand_int_double(type2);
	  z=x*y;
	  z=z/x;
	  z=z/y;
	  z=z-1;
	  if(!contains_zero(z))
	    {
	      my_print("z",z);
	      my_print("x",x);
	      my_print("y",y);
	      my_print("x*y",z=x*y);
	      my_print("x*y/x",z=z/x);
	      my_print("x/y/x/y",z=z/y);
	      return(false);
	    }
	}
  return(true);
}



//(x+y)(x-y)=x^2-y^2
bool testb (unsigned int n)
{
  int_double x,y;
  unsigned int i;
  char type1,type2,type3,type4;
  
  for(type1=POS;type1<=CROSS_ZERO;type1++)
    for(type2=POS;type2<=CROSS_ZERO;type2++)
      for(i=0;i<n;i++)
	{
	  x=rand_int_double(type1);
	  y=rand_int_double(type2);
	  if(!contains_zero(sqr(x)-sqr(y)-(x+y)*(x-y)))
	    {
	      printf("Test number %d.\n",i);
	      my_print("x",x);
	      my_print("y",y);
	      my_print("x*x",x*x);
	      my_print("y*y",y*y);
	      my_print("x*x-y*y",sqr(x)-sqr(y));
	      my_print("x+y",x+y);
	      my_print("x-y",x-y);
	      my_print("(x+y)*(x-y)",(x+y)*(x-y));
	      return false;
	    }
	}
  
  return(true);
}

//(x+y)(x-y)=x^2-y^2 with y a double
bool testc (unsigned int n)
{
  int_double x;
  double y;
  unsigned int i;
  char type1,type2,type3,type4;
  
  for(type1=POS;type1<=CROSS_ZERO;type1++)
    for(type2=POS;type2<=TO_ZERO;type2++)
      for(i=0;i<n;i++)
	{
	  x=rand_int_double(type1);
	  y=rand_double();
	  if(type2==NEG)
	    y=-y;
	  if(type2==TO_ZERO)
	    y=0.0;
	  int_double y2=sqr(y);
	  if(!contains_zero(sqr(x)-y2-(x+y)*(x-y)))
	    {
	      printf("Test number %d.\n",i);
	      my_print("x",x);
	      printf("y %20.18e\n",y);
	      my_print("x*x",sqr(x));
	      my_print("y*y",y2);
	      my_print("x*x-y*y",sqr(x)-y2);
	      my_print("x+y",x+y);
	      my_print("x-y",x-y);
	      my_print("(x+y)*(x-y)",(x+y)*(x-y));
	      return(false);
	    }
	}
  
  return(true);
}

//sqrt(x)^2=x
bool testd (unsigned int n)
{
  for(int i=0;i<n;i++)
    {
      int_double x=rand_int_double(POS);
      int_double sx=sqrt(x);
      if(!contains_zero(sqr(sx)-x))
	{
	  printf("Test number %d.\n",i);
	  my_print("x",x);
	  my_print("sqrt(x)",sx);
	  my_print("sqrt(x)^2-x",sqr(sx)-x);
	  return(false);
	}
      x=rand_int_double(FROM_ZERO);
      sx=sqrt(x);
      if(!contains_zero(sqr(sx)-x))
	{
	  printf("Test number %d.\n",i);
	  my_print("x",x);
	  my_print("sqrt(x)",sx);
	  my_print("sqrt(x)^2-x",sqr(sx)-x);
	  return(false);
	}
    }
  return(true);
}
	



//sin^2(x)+cos^2(x)=1
bool test3(unsigned int n)
{
  int_double x,cos_x,sin_x;
  unsigned int i;
  char type1;
  for(x.left=0.0;x.left<M_PI/2;x.left=x.left+1.0/(double) n)
    {
      x.right=x.left+0.00000001;
      sin_cos(x,&sin_x,&cos_x);
      if(!contains_zero(sin_x*sin_x+cos_x*cos_x-1))
	{
	  printf("x: ");
	  print_int_double(x);
	  printf("\nsin(x): ");
	  print_int_double(sin_x);
	  printf("\ncos(x): ");
	  print_int_double(cos_x);
	  printf("\n");
	  print_int_double_str("sin^2+cos^2-1=",sin_x*sin_x+cos_x*cos_x-1);
	  return(false);
	}
    }
  return(true);
}

//log(exp(x))=x
bool test4 (unsigned int n)
{
  int_double x;
  unsigned int i;
  char type1;
  for(type1=POS;type1<=CROSS_ZERO;type1++)
    for(i=0;i<n;i++)
      {
	x=rand_int_double(type1);
	if(!contains_zero(log(exp(x))-x))
	  return(false);
      }
  return(true);
}

//exp(log(x))=x
bool test5 (unsigned int n)
{
  int_double x,lnx,explnx;
  unsigned int i;
  for(i=0;i<n;i++)
    {
      x=rand_int_double(POS);
      lnx=log(x);
      explnx=exp(log(x));
      if(!contains_zero(explnx-x))
	{
	  printf("i: %d\nx: ",i);
	  print_int_double(x);
	  printf("\nexp(log(x))-x: ");
	  print_int_double(exp(log(x))-x);
	  return(false);
	}
    }
  return(true);
}

//atan2(sin(x),cos(x))=x
bool test6 (unsigned int n)
{
  double theta0;
  int_double theta,si,co,atn;
  for(theta0=-3.14/2.0;theta0<=3.13;theta0+=1.0/(double) n)
    {
      theta=int_double(theta0,theta0+1.0/100.0);
      sin_cos(theta,&si,&co);
      atn=atan2(si,co);
      if(!(contains_zero(atn-theta)))
	{
	  print_int_double_str("theta=",theta);
	  print_int_double_str("sin  =",si);
	  print_int_double_str("cos  =",co);
	  print_int_double_str("atan =",atn);
	  return(false);
	}
    }
  return(true);
}

// cos(2x)=2cos^2(x)-1
bool test7 (unsigned int n)
{
  double theta0;
  int_double theta,si,co1,co2;
  for(theta0=-3.14;theta0<=3.13;theta0+=1.0/(double) n)
    {
      theta=int_double(theta0,theta0+1.0/100.0);
      sin_cos(theta*2,&si,&co1);
      sin_cos(theta,&si,&co2);
      if(!(contains_zero(co1+1.0-2*sqr(co2))))
	{
	  print_int_double_str("theta=",theta);
	  print_int_double_str("co1  =",co1);
	  print_int_double_str("co2  =",co2);
	  return(false);
	}
    }
  return(true);
}

int main()
{
	printf("Testing....\n");
	print_int_double_str("Pi=",d_pi);
	print_int_double_str("Pi/2=",d_pi_2);
	print_int_double_str("2 Pi=",d_two_pi);
	print_int_double_str("sqrt(2)=",sqrt(int_double(2)));
	print_int_double_str("[-1,2]^2=",sqr(int_double(-1,2)));
	if(!testa(100000))
	    printf("testa failed.\n");
	else
	    printf("testa passed.\n");

	if(!testb(100000))
	    printf("testb failed.\n");
	else
	    printf("testb passed.\n");

	if(!testc(100000))
	    printf("testc failed.\n");
	else
	    printf("testc passed.\n");

	if(!testd(100000))
	    printf("testd failed.\n");
	else
	    printf("testd passed.\n");
	
	
	if(!test3(100000))
	    printf("test3 failed.\n");
	else
	    printf("test3 passed.\n");
	
	if(!test4(100000))
	    printf("test4 failed.\n");
	else
	    printf("test4 passed.\n");
	
	  if(!test5(100000))
	    printf("test5 failed.\n");
	else
	    printf("test5 passed.\n");
	  
	if(!test6(100000))
	    printf("test6 failed.\n");
	else
	    printf("test6 passed.\n");
	if(!test7(100000))
	    printf("test7 failed.\n");
	else
	    printf("test7 passed.\n");
	return(0);
}

